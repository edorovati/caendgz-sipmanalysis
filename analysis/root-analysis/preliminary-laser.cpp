#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TList.h>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

#include "root-library/utils.h"
#include "root-library/waveform.h"
#include "root-library/gauss_stuff.h"
#include "root-library/timestuff.h"
#include "root-library/graphics.h"
#include "root-library/analysis.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TPaletteAxis.h>
#include <TStyle.h>

void starter(const std::vector<std::tuple<TString, TString, TString>>& file_list,
             const TString& out_filename)
{
    // ========================
    // PARAMETRI GLOBALI
    // ========================

    // Range temporale laser
    const double T_MIN = 38.0;
    const double T_MAX = 42.0;

    // Parametri istogrammi di ampiezza (reference)
    const int    N_BINS_MAX_AMP       = 4000;
    const double LOWER_RANGE_MAX_AMP  = -200.0;
    const double UPPER_RANGE_MAX_AMP  = 200.0;

    // Parametri istogrammi 2D (MaxAmp vs DeltaT)
    const int    N_BINS_MAXAMP_2D     = 1000;
    const double RANGE_MAXAMP_MIN     = 0.0;
    const double RANGE_MAXAMP_MAX     = 200.0;
    const int    N_BINS_DELTAT_2D     = 40000.0;
    const double RANGE_DELTAT_MIN     = -1000.0;
    const double RANGE_DELTAT_MAX     = 1000.0;

    // Parametri calcolo soglie
    const double LASER_THRESHOLD      = 70.0;   // mV
    const int    N_SAMPLES_PER_WF     = 275;
    const double N_SIGMA_NOISE_THR    = 5.0;

    // Stile grafico per istogrammi / proiezioni
    const Color_t LINE_COLOR_REF = kGreen-2;
    const int      MARKER_STYLE_REF = 20;
    const int      LINE_WIDTH_REF = 2;
    const int      FILL_STYLE_NICE = 3004;
    const double   FILL_ALPHA = 0.35;

    // ========================
    // CREAZIONE FILE DI OUTPUT
    // ========================
    auto outFile = std::unique_ptr<TFile>(TFile::Open(out_filename, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: cannot create output file " << out_filename << std::endl;
        return;
    }

    TDirectory* dir_ref_amps     = outFile->mkdir("reference_amplitude_distributions");
    TDirectory* dir_2d_deltat    = outFile->mkdir("deltaT_vs_maxAmp");
    TDirectory* dir_deltat_proj  = outFile->mkdir("deltaT_projections");
    TDirectory* dir_summary      = outFile->mkdir("summary_vs_vbias");

    // Per salvare le soglie calcolate dalle reference per ogni canale
    std::map<TString, double> threshold_map;

    // Per memorizzare risultati vs vbias per ogni canale
    struct SummaryVecs {
        std::vector<double> vbias_v;      // vbias numeric
        std::vector<double> val;          // probability o dcr
        std::vector<double> val_err;
    };
    std::map<TString, SummaryVecs> prob_vs_vbias; // per ogni canale: probability
    std::map<TString, SummaryVecs> dcr_vs_vbias;  // per ogni canale: dcr (kHz)

    // ========================
    // LOOP SUI FILE DI INPUT
    // ========================
    for (const auto& [filename, file_type, vbias] : file_list) {

        std::cout << "Processing file: " << filename.Data()
                  << ", type: " << file_type.Data()
                  << ", vbias: " << vbias.Data() << std::endl;

        auto inFile = std::unique_ptr<TFile>(TFile::Open(filename));
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: cannot open input file " << filename << std::endl;
            continue;
        }

        // ========================
        // CARICAMENTO DATI LASER
        // ========================
        TTree* tree_laser = dynamic_cast<TTree*>(inFile->Get("laser"));
        std::vector<std::pair<std::vector<double>, std::vector<double>>> laser_waveforms;
        Double_t *laser_amp = nullptr, *laser_time = nullptr;

        if (tree_laser) {
            laser_amp  = Utils::setupBranch_dgz(tree_laser, "amplitude");
            laser_time = Utils::setupBranch_dgz(tree_laser, "time");
            Long64_t nEntries = tree_laser->GetEntries();
            laser_waveforms.reserve(nEntries);

            for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
                tree_laser->GetEntry(iEntry);
                std::vector<double> time_vec(1024), amp_vec(1024);
                for (int i = 0; i < 1024; ++i) {
                    time_vec[i] = laser_time[i];
                    amp_vec[i]  = laser_amp[i];
                }
                laser_waveforms.emplace_back(std::move(time_vec), std::move(amp_vec));
            }
            std::cout << "Loaded " << laser_waveforms.size() << " laser waveforms\n";
        } else {
            std::cerr << "No laser tree found in file " << filename << std::endl;
            continue;
        }

        // ========================
        // CARICAMENTO CANALI
        // ========================
        std::map<TString, TTree*> tree_map;
        std::map<TString, Double_t*> time_map;
        std::map<TString, Double_t*> amp_map_raw;

        if (auto keys = inFile->GetListOfKeys()) {
            for (int i = 0; i < keys->GetEntries(); ++i) {
                TString name = keys->At(i)->GetName();
                if (name.BeginsWith("ch")) {
                    if (auto tree = dynamic_cast<TTree*>(inFile->Get(name))) {
                        auto amp  = Utils::setupBranch_dgz(tree, "amplitude");
                        auto time = Utils::setupBranch_dgz(tree, "time");
                        if (amp && time) {
                            tree_map[name]    = tree;
                            amp_map_raw[name] = amp;
                            time_map[name]    = time;
                            std::cout << "Loaded tree and branches for " << name << std::endl;
                        }
                    }
                }
            }
        }

        // ========================
        // ALLINEAMENTO & BASELINE
        // ========================
        Waveform waveAligner(tree_laser, laser_time, laser_amp);
        std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_baseline_corrected_waveforms;

        for (const auto& [ch_name, ch_tree] : tree_map) {
            auto ch_time = time_map[ch_name];
            auto ch_amp  = amp_map_raw[ch_name];
            if (ch_tree && ch_time && ch_amp) {
                auto raw           = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, false);
                auto corrected_raw = waveAligner.correctBaseline(raw);
                all_baseline_corrected_waveforms[ch_name] = std::move(corrected_raw);
            }
        }

        // ========================
        // ANALISI REFERENCE
        // ========================
        if (file_type == "reference") {
            dir_ref_amps->cd();
            for (const auto& [ch_name, waveforms] : all_baseline_corrected_waveforms) {

                // Nome e titolo istogramma
                TString hist_name = Form("%s-%s-%s-baseline", file_type.Data(), vbias.Data(), ch_name.Data());
                TH1F h_ref(hist_name, hist_name + ";Amplitude (mV);Counts",
                           N_BINS_MAX_AMP, LOWER_RANGE_MAX_AMP, UPPER_RANGE_MAX_AMP);

                // Riempimento
                for (const auto& [times_vec, amps_vec] : waveforms) {
                    for (const auto& amp : amps_vec) {
                        h_ref.Fill(amp);
                    }
                }

                // ========= canvas + stile =========
                TString canvas_name = Form("c_ref_%s_%s", vbias.Data(), ch_name.Data());
                TCanvas* c_ref = Graphics::CreateCanvas(canvas_name, hist_name, 1200, 800, false, true, true);
                c_ref->cd();

                

                // Format assi
                Graphics::FormatAxis(h_ref.GetXaxis(), h_ref.GetYaxis(),
                                     "amplitude (mV)", "counts",
                                     0.03, 0.03, 62,
                                     true, -10, 10,
                                     false, 0, 0, 42);

                // Fit gaussiano (fit dopo il draw per avere stat box aggiornato)
                TF1 fgaus("fgaus", "gaus", -5, 5);
                fgaus.SetNpx(1000);
                fgaus.SetLineColor(kRed+2);
                fgaus.SetLineStyle(2);
                fgaus.SetLineWidth(2);
                h_ref.Fit(&fgaus, "IMREQ");
                h_ref.Draw("SAME");
                
                float mu    = fgaus.GetParameter(1);
                float sigma = fgaus.GetParameter(2);
                double noise_thr = mu + N_SIGMA_NOISE_THR * sigma;

                // Calcolo probabilità falsi crossing (solo info)
                float p_tail = 0.5f * std::erfc((noise_thr - mu) / (std::sqrt(2.0f) * sigma));
                float p_fake_wf = 1.0f - std::pow(1.0f - p_tail, N_SAMPLES_PER_WF);
                double expected_fake_crossings = p_fake_wf * waveforms.size();

                threshold_map[ch_name] = noise_thr;

                // Linea soglia (verticale)
                TLine* thr_line = Graphics::CreateLine(noise_thr, 0.0, noise_thr, h_ref.GetMaximum() * 1.05, kMagenta+2, 6, 2);
                thr_line->Draw("SAME");
                // Setup estetico e draw istogramma
                Graphics::SetupHistogramStyle(&h_ref,
                                             true,           // showStats
                                             1111,           // statOptions
                                             LINE_COLOR_REF, // lineColor
                                             LINE_COLOR_REF, // markerColor
                                             MARKER_STYLE_REF,
                                             LINE_WIDTH_REF,
                                             "HISTE1",         // drawOption
                                             true,           // enableFillAlpha
                                             LINE_COLOR_REF,
                                             FILL_STYLE_NICE,
                                             FILL_ALPHA);
                // Legenda
                TLegend* leg = Graphics::CreateLegend(0.1, 0.7, 0.3, 0.9);
                Graphics::AddLegendEntry(leg, &h_ref, "amplitude distribution", "l");
                Graphics::AddLegendEntry(leg, &fgaus, "gaussian fit", "l");
                Graphics::AddLegendEntry(leg, thr_line, Form("threshold = %.2f mV", noise_thr), "l");
                leg->Draw();

                // TPave con formula threshold e info aggiuntive
                TString thr_formula;
                thr_formula.Form("thr = #mu + %.0f#times#sigma = %.2f mV",
                                 N_SIGMA_NOISE_THR, noise_thr);

                TPaveText* pave_info = Graphics::CreateInfoBox(0.1, 0.1, 0.3, 0.2, thr_formula, kWhite, kBlack, 0.03, 1);
                pave_info->Draw();

                // Salvo la canvas nel file ROOT
                c_ref->Write();

                // Pulizia oggetti dinamici
                delete thr_line;
                delete leg;
                delete pave_info;
                delete c_ref;
            }
            outFile->cd();
        }

        // ========================
        // ANALISI DATA
        // ========================
        if (file_type == "data") {
            dir_2d_deltat->cd();

            for (const auto& [ch_name, waveforms] : all_baseline_corrected_waveforms) {
                if (threshold_map.find(ch_name) == threshold_map.end()) {
                    std::cerr << "No threshold found for channel " << ch_name << std::endl;
                    continue;
                }
                double threshold = threshold_map[ch_name];

                TString hist2d_name = Form("h2_deltaT_vs_maxAmp_%s_%s", vbias.Data(), ch_name.Data());
                TH2F h_fixed_th_diff(hist2d_name, hist2d_name + ";Max Amp;Delta t (ns)",
                                     N_BINS_MAXAMP_2D, RANGE_MAXAMP_MIN, RANGE_MAXAMP_MAX,
                                     N_BINS_DELTAT_2D, RANGE_DELTAT_MIN, RANGE_DELTAT_MAX);

                // Riempio l'istogramma 2D
                size_t n_fill = std::min(waveforms.size(), laser_waveforms.size());
                for (size_t i = 0; i < n_fill; ++i) {
                    const auto& [laser_times, laser_amps] = laser_waveforms[i];
                    const auto& [times_vec, amps_vec]     = waveforms[i];

                    double crossing_time_laser = TimeStuff::computeCrossingTime(laser_times.data(), laser_amps.data(),
                                                                                 (int)laser_times.size(), LASER_THRESHOLD);
                    if (crossing_time_laser < 0) continue;

                    double crossing_time_ch = TimeStuff::computeCrossingTime(times_vec.data(), amps_vec.data(),
                                                                             (int)times_vec.size(), threshold);
                    if (crossing_time_ch >= 0) {
                        int crossing_index = -1;
                        for (int j = 1; j < (int)times_vec.size(); ++j) {
                            if (times_vec[j-1] < crossing_time_ch && times_vec[j] >= crossing_time_ch) {
                                crossing_index = j;
                                break;
                            }
                        }
                        if (crossing_index != -1) {
                            int max_idx   = TimeStuff::findMaxAfterCrossing(amps_vec, crossing_index, threshold);
                            double max_amp = amps_vec[max_idx];
                            double delta_t = crossing_time_ch - crossing_time_laser;
                            h_fixed_th_diff.Fill(max_amp, delta_t);
                        }
                    }
                } // end fill

                // Salvo TH2F su file
                
                // ========== Canvas 2D ==============
                TString c2d_name = Form("c2d_%s_%s", vbias.Data(), ch_name.Data());
                TCanvas* c2d = Graphics::CreateCanvas(c2d_name, hist2d_name, 1000, 700, false, false, true);
                c2d->cd();

                // disegno con l'estetica definita
                Graphics::Format2DHisto(&h_fixed_th_diff, kViridis, "COLZ", 0.92, 0.025, 0.45, false);
                // linea soglia verticale (sul asse X = Max Amp)
                TLine* thr_line_2d = Graphics::CreateLine(threshold, 30, threshold, 50, kRed+1, 2, 2);
                thr_line_2d->Draw("SAME");
                Graphics::FormatAxis(h_fixed_th_diff.GetXaxis(), h_fixed_th_diff.GetYaxis(),
                                     "amplitude (mV)", "#Deltat (ns)");
                // piccolo TPave con info threshold
                TPaveText* info2d = Graphics::CreateInfoBox(0.7, 0.7, 0.9, 0.9,
                                                            Form("%s @v_{bias}=%sV\n @thr=%.2f mV", ch_name.Data(), vbias.Data(), threshold),
                                                            kWhite, kBlack, 0.03, 1);
                info2d->Draw();

                c2d->Write();
                delete thr_line_2d;
                delete info2d;
                delete c2d;

                // ========================
                // PROIEZIONE SU Y (Δt)
                // ========================
                dir_deltat_proj->cd();
                TH1D* h_projY = h_fixed_th_diff.ProjectionY(Form("projY_%s_%s", vbias.Data(), ch_name.Data()), 1, -1, "E"); // projection of all x
                h_projY->SetTitle(Form("Delta t projection for %s %s;Delta t (ns);Counts", vbias.Data(), ch_name.Data()));

                // canvas per la proiezione
                TString cproj_name = Form("c_projY_%s_%s", vbias.Data(), ch_name.Data());
                TCanvas* cproj = Graphics::CreateCanvas(cproj_name, h_projY->GetTitle(), 1000, 700, false, false, true);
                cproj->cd();

                // disegno proiezione con stile
                Graphics::SetupHistogramStyle(h_projY, false, 1111, kBlue, kBlue, 47, 1, "HISTE1", true, kBlue, FILL_STYLE_NICE, 0.35);
                Graphics::FormatAxis(h_projY->GetXaxis(), h_projY->GetYaxis(), "#Deltat (ns)", "counts", 0.03, 0.03, 62, true, 30, 50);

                // trova bin più alto
                // Trova bin più alto
                int b = h_projY->GetMaximumBin();
                Analysis::SignalBackground sb = Analysis::CalculateSignalAndBackground(h_projY, b);
                int n_wf = static_cast<int>(std::min(waveforms.size(), laser_waveforms.size()));
                auto [prob, prob_err] = Analysis::CalculateProbability(sb.diff, n_wf);

                // --- Dark counts: finestra temporale fino a (tempo_picco - 2 ns) ---
                double time_peak = h_projY->GetXaxis()->GetBinCenter(b);
                double time_limit = time_peak - 2.0; // 2 ns prima del picco

                // Trova il bin corrispondente a time_limit
                int bin_end_dark = h_projY->GetXaxis()->FindBin(time_limit);
                int bin_start_dark = 1;
                int nbins = h_projY->GetNbinsX();

                for (int i = 1; i <= nbins; ++i) {
                    if (h_projY->GetBinContent(i) > 0) {
                        bin_start_dark = i;
                        break;
                    }
                }

                double dark_counts = Utils::CountEvents(h_projY, bin_start_dark, bin_end_dark);

                // Finestra temporale effettiva (in ns) come somma dei moduli
                double time_min = h_projY->GetXaxis()->GetBinLowEdge(bin_start_dark);
                double time_max = h_projY->GetXaxis()->GetBinLowEdge(bin_end_dark + 1);
                double time_window = std::fabs(time_min) + std::fabs(time_max);

                // Calcolo DCR (time_window in secondi)
                auto [dcr, dcr_err] = Analysis::CalculateDCR(dark_counts, time_window * 1e-9, n_wf);
                
                // TPave con risultati
                TString res_text;
                res_text.Form("%s @v_{bias}=%sV\n p_{signal} = %.2f #pm %.2f %%\n DCR = %.0f #pm %.0f kHz\n",
                              ch_name.Data(), vbias.Data(), prob*100, prob_err*100, dcr, dcr_err);

                TPaveText* pave_res = Graphics::CreateInfoBox(0.6, 0.7, 0.9, 0.9, res_text, kWhite, kBlack, 0.03, 1);
                pave_res->Draw();

                // salva canvas e hist
                cproj->Write();

                // Memorizza risultati vs vbias per plotting successivo
                double vbias_val = vbias.Atof(); // converte la stringa vbias in double (es. "100" -> 100)
                // Probability
                prob_vs_vbias[ch_name].vbias_v.push_back(vbias_val);
                prob_vs_vbias[ch_name].val.push_back(prob*100);
                prob_vs_vbias[ch_name].val_err.push_back(prob_err*100);
                // DCR
                dcr_vs_vbias[ch_name].vbias_v.push_back(vbias_val);
                dcr_vs_vbias[ch_name].val.push_back(dcr);
                dcr_vs_vbias[ch_name].val_err.push_back(dcr_err);

                // pulizia
                delete pave_res;
                delete cproj;
                delete h_projY;

                dir_2d_deltat->cd(); // torniamo nella cartella principale 2d
            } // end per canale
            outFile->cd();
        } // end if data
    } // end loop files

    // ========================
    // CREAZIONE GRAFICI SUMMARY VS VBIAS (per canale)
    // ========================
    dir_summary->cd();
    for (const auto& kv : prob_vs_vbias) {
        const TString ch_name = kv.first;
        const SummaryVecs& pv = kv.second;

        if (pv.vbias_v.empty()) continue;

        // Grafico probability vs vbias
        TGraphErrors* g_prob = new TGraphErrors((int)pv.vbias_v.size());
        for (size_t i = 0; i < pv.vbias_v.size(); ++i) {
            g_prob->SetPoint((int)i, pv.vbias_v[i], pv.val[i]);
            g_prob->SetPointError((int)i, 0.0, pv.val_err[i]);
        }
        // Formatta assi
        Graphics::FormatAxis(g_prob->GetXaxis(), g_prob->GetYaxis(), "v_{bias} (V)", "probability (%)");
        TCanvas* cgp = Graphics::CreateCanvas(Form("c_prob_%s", ch_name.Data()), Form("Probability %s", ch_name.Data()), 900, 600, false, false, true);
        cgp->cd();
        Graphics::CustomizeGraph(g_prob, 47, kBlue, "APE", kBlue, 6, 1);
        g_prob->Write(Form("g_prob_%s", ch_name.Data()));

        // Grafico DCR vs vbias
        const SummaryVecs& dv = dcr_vs_vbias[ch_name];
        if (!dv.vbias_v.empty()) {
            TGraphErrors* g_dcr = new TGraphErrors((int)dv.vbias_v.size());
            for (size_t i = 0; i < dv.vbias_v.size(); ++i) {
                g_dcr->SetPoint((int)i, dv.vbias_v[i], dv.val[i]);
                g_dcr->SetPointError((int)i, 0.0, dv.val_err[i]);
            }
            Graphics::FormatAxis(g_dcr->GetXaxis(), g_dcr->GetYaxis(), "v_{bias} (V)", "DCR (kHz)");
            TCanvas* cgd = Graphics::CreateCanvas(Form("c_dcr_%s", ch_name.Data()), Form("DCR %s", ch_name.Data()), 900, 600, false, false, true);
            cgd->cd();
            Graphics::CustomizeGraph(g_dcr, 47, kRed, "APE", kRed, 6, 1);
            g_dcr->Write(Form("g_dcr_%s", ch_name.Data()));

            delete cgd;
            delete g_dcr;
        }

        delete cgp;
        delete g_prob;
    }

    // ========================
    // CHIUSURA FILE
    // ========================
    outFile->Write();
    outFile->Close();

    std::cout << "Processing finished. Output saved to " << out_filename.Data() << std::endl;
}
