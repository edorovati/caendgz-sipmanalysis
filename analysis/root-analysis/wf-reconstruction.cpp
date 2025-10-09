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
#include "root-library/fit.h"

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

void waveform(const std::vector<std::tuple<TString, TString, TString>>& file_list,
              const TString& out_root_name)
{
    const double t_min = 38.0;   // range minimo temporale laser
    const double t_max = 42.0;   // range massimo temporale laser

    const int    n_bins_max_amp      = 1000;   // numero di bin
    const double lower_range_max_amp = 0.0;    // valore minimo asse x
    const double upper_range_max_amp = 200.0;  // valore massimo asse x

    const int n_peaks_requested = 5; // numero picchi per fit multi-Gauss

    auto outFile = std::unique_ptr<TFile>(TFile::Open(out_root_name, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: cannot create output file " << out_root_name << std::endl;
        return;
    }

    TDirectory* dir_amp   = outFile->mkdir("amp_max-distribution");
    TDirectory* dir_cal   = outFile->mkdir("calibration");
    TDirectory* dir_avgWF = outFile->mkdir("average_waveforms");
    TDirectory* dir_histowf = outFile->mkdir("histo_waveforms");
    TDirectory* dir_base = outFile->mkdir("baseline");
    
    
    static TTree* info_tree = nullptr;
    static double tree_pe = 0.0;
    static double tree_slope = 0.0;
    static double tree_slope_err = 0.0;
    static double tree_baseline = 0.0;
    static double tree_baseline_err = 0.0;
    static double tree_jitter = 0.0;
    static double tree_jitter_err = 0.0;

    if (!info_tree) {
        info_tree = new TTree("info", "Waveform analysis info");
        info_tree->Branch("pe", &tree_pe);
        info_tree->Branch("slope", &tree_slope);
        info_tree->Branch("slope_err", &tree_slope_err);
        info_tree->Branch("baseline", &tree_baseline);
        info_tree->Branch("baseline_err", &tree_baseline_err);
        info_tree->Branch("jitter", &tree_jitter);
        info_tree->Branch("jitter_err", &tree_jitter_err);
    }

    
    std::map<TString, TF1*> linearFitsMap;
    struct IntegralPoint {
        double vbias;
        double value;
        double error;
    };

    // Aggiorni la mappa così:
    std::map<TString, std::vector<IntegralPoint>> integralsMap;
    std::map<TString, std::vector<IntegralPoint>> tauMap;
    std::map<TString, std::vector<IntegralPoint>> tauMap_single;  // chiave = canale

    // Per ogni file
    for (const auto& [filename, file_type, vbias] : file_list) {
        auto inFile = std::unique_ptr<TFile>(TFile::Open(filename));
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: cannot open input file " << filename << std::endl;
            continue;
        }
        
        TTree* tree_laser = dynamic_cast<TTree*>(inFile->Get("laser"));
        if (!tree_laser) {
            std::cerr << "No 'laser' tree in file " << filename << std::endl;
            continue;
        }
        Double_t* laser_amp  = Utils::setupBranch_dgz(tree_laser, "amplitude");
        Double_t* laser_time = Utils::setupBranch_dgz(tree_laser, "time");
        
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
                        }
                    }
                }
            }
        }
        
        Waveform waveAligner(tree_laser, laser_time, laser_amp);
        std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_aligned_waveforms;
        
        for (const auto& [ch_name, ch_tree] : tree_map) {
            auto ch_time = time_map[ch_name];
            auto ch_amp  = amp_map_raw[ch_name];
            if (ch_tree && ch_time && ch_amp) {
                auto raw = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, true);
                auto corrected_raw = waveAligner.correctBaseline(raw);
                all_aligned_waveforms[ch_name] = std::move(corrected_raw);
            }
        }
        
        // Calibrazione e fit
        for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
            TString unique_hist_name = TString::Format("%s-%s", vbias.Data(), ch_name.Data());
            
            auto [hist_noise, hist_signal, hist_diff] = waveAligner.calibration(unique_hist_name, waveforms, t_min, t_max, n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
            if (!hist_diff) continue;
            
            if (hist_signal->Integral() < 10) {
                std::cout << "Histogram too empty for " << vbias << " " << ch_name << ", skipping fit." << std::endl;
                continue;
            }
            
            dir_amp->cd();
            TString fit_name = TString::Format("%s-%s-multi-gauss", vbias.Data(), ch_name.Data());
            TF1* fit = GaussStuff::fitHistogram(hist_signal.get(), n_peaks_requested, fit_name);
            //hist_diff->Write();
            //hist_noise->Write();
            hist_signal->Write();
            if (fit) {
                struct PeakInfo {
                    double amp;
                    double mean;
                    double mean_err;
                };
                
                int n_par   = fit->GetNpar();
                int n_peaks = n_par / 3;
                std::vector<PeakInfo> peaks;
                
                for (int i = 0; i < n_peaks; ++i) {
                    peaks.push_back({
                        fit->GetParameter(i*3+0),
                        fit->GetParameter(i*3+1),
                        fit->GetParError(i*3+1)
                    });
                }
                
                std::sort(peaks.begin(), peaks.end(),
                          [](const PeakInfo& a, const PeakInfo& b) {
                    return a.mean < b.mean;
                });
                
                std::vector<double> x_vals, y_vals, y_errs;
                
                // Ogni picco corrisponde al numero di fotoelettroni
                for (size_t i = 0; i < peaks.size(); ++i) {
                    double x    = static_cast<double>(i);   // 0 p.e., 1 p.e., 2 p.e., ...
                    double y    = peaks[i].mean;
                    double yerr = peaks[i].mean_err;
                    
                    x_vals.push_back(x);
                    y_vals.push_back(y);
                    y_errs.push_back(yerr);
                    
                    std::cout << "Peak " << i
                    << " -> x = " << x
                    << ", y = " << y
                    << " ± " << yerr
                    << std::endl;
                }
                
                TString graph_name = TString::Format("%s-%s_mean_vs_pe", vbias.Data(), ch_name.Data());
                TGraphErrors* graph = new TGraphErrors(
                                                       x_vals.size(),
                                                       x_vals.data(),
                                                       y_vals.data(),
                                                       nullptr,
                                                       y_errs.data()
                                                       );
                graph->SetName(graph_name);
                graph->SetTitle(TString::Format(
                                                "Mean vs Photo-electron number for %s %s;PE;Mean (charge units)",
                                                vbias.Data(), ch_name.Data()
                                                ));
                
                TString linear_fit_name = TString::Format("%s-%s_linearFit", vbias.Data(), ch_name.Data());
                TF1* linearFit = new TF1(linear_fit_name, "pol1", 0, n_peaks);
                graph->Fit(linearFit, "IMREQ");
                
                linearFitsMap[vbias + "-" + ch_name] = linearFit;
                
                dir_cal->cd();
                graph->Write();
            }
            else {
                // Fit fallito → fit fittizio
                TString fake_fit_name = TString::Format("%s-%s_linearFit_fake", vbias.Data(), ch_name.Data());
                TF1* fakeFit = new TF1(fake_fit_name, "pol1", 0, 2);
                fakeFit->SetParameter(0, 0.0);
                fakeFit->SetParameter(1, 0.0);
                linearFitsMap[vbias + "-" + ch_name] = fakeFit;
                
                std::cout << "⚠️  No valid peak fit for " << vbias << " " << ch_name
                << " → using slope=0, intercept=0" << std::endl;
            }
        }
        // Selezione waveform e media con ottimizzazione tramite ricerca binaria
        // Loop su tutti i canali e waveform
        // Loop su tutti i canali e waveform
        // ============================
        // Ciclo sulle channel
        // ============================
        for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
            auto keyFit = vbias + "-" + ch_name;
            auto fitIt = linearFitsMap.find(keyFit);
            if (fitIt == linearFitsMap.end()) continue;

            double slope     = fitIt->second->GetParameter(1);
            double intercept = fitIt->second->GetParameter(0);

            // ============================
            // 1) Distribuzione ampiezze pre-35 ns
            // ============================
            TH1D* h_pre35 = new TH1D(Form("h_pre35_%s_%s", vbias.Data(), ch_name.Data()),
                                      Form("Amplitude distribution t<35ns %s %s", vbias.Data(), ch_name.Data()),
                                      200, -20, 20);

            for (const auto& wf : waveforms) {
                for (size_t i=0; i<wf.first.size(); i++) {
                    if (wf.first[i] < 35.0) h_pre35->Fill(wf.second[i]);
                }
            }

            auto [baseline_val, baseline_err, baseline_sigma, baseline_sigma_err, baseline_fit] =
                GaussStuff::gauss_fit_projection(h_pre35);
            dir_base->cd();
            h_pre35->Write();

            // ============================
            // 2) Analisi con più soglie
            // ============================
            std::vector<std::pair<double,double>> thr_ranges = {
                {0.5, 1.5}, {1.5, 2.5}, {2.5, 3.5}, {3.5, 4.5}, {4.5, 5.5}, {5.5, 6.5},  {6.5, 7.5}, {7.5, 8.5}
            };
            
            size_t pe_idx = 0;
            for (auto [low_mult, up_mult] : thr_ranges) {
                double lower_thr = low_mult * slope + intercept;
                double upper_thr = up_mult * slope + intercept;
                double thr = 0.5 * slope + intercept;

                // Calcola la soglia fissa direttamente qui
                double fixed_threshold = (low_mult - 0.0) * slope + intercept;  // o la logica che hai definito
                if (low_mult == 0.5) fixed_threshold = 0.5*slope + intercept;
                else if (low_mult == 1.5) fixed_threshold = 1.0*slope + intercept;
                else if (low_mult == 2.5) fixed_threshold = 1.5*slope + intercept;
                else if (low_mult == 3.5) fixed_threshold = 2.0*slope + intercept;

                std::cout << "th-fixed " << fixed_threshold << std::endl;
                
                std::cout << lower_thr << " SOGLIE " << upper_thr << "\n";

                double crossing_time_avg = 0.0;
                double crossing_sigma_avg = 0.0;

                // ----------------------------
                // Selection con crossing e massimo
                // ----------------------------
                auto selected = Waveform::selection(waveforms, lower_thr, upper_thr, thr,
                                                    crossing_time_avg, crossing_sigma_avg);
                if (selected.empty()) continue;

                

                // ----------------------------
                // Calcolo waveform media
                // ----------------------------
                auto avgWF = Waveform::computeAverageWaveform(selected);

                // ----------------------------
                // Trova massimo della media
                // ----------------------------
                auto maxIt = std::max_element(avgWF.amplitudes.begin(), avgWF.amplitudes.end());
                double max_amp = *maxIt;
                int peakIdx = std::distance(avgWF.amplitudes.begin(), maxIt);
                double max_time = avgWF.times[peakIdx];

                std::cout << Form("Bias %s, Ch %s, Thr %.1f-%.1f: MaxAmp=%.3f mV @ %.3f ns",
                                  vbias.Data(), ch_name.Data(), low_mult, up_mult, max_amp, max_time)
                          << std::endl;

                // ----------------------------
                // Calcolo rise time 20%-80%
                // ----------------------------
                const double level80 = 0.8 * max_amp;
                const double level20 = 0.2 * max_amp;

                double t80 = NAN, t20 = NAN;
                int i80 = -1, i20 = -1;

                for (int i = peakIdx; i > 0; --i) {
                    if (avgWF.amplitudes[i-1] < level80 && avgWF.amplitudes[i] >= level80) {
                        t80 = avgWF.times[i-1] + (level80 - avgWF.amplitudes[i-1]) *
                              (avgWF.times[i] - avgWF.times[i-1]) / (avgWF.amplitudes[i] - avgWF.amplitudes[i-1]);
                        i80 = i;
                        break;
                    }
                }
                if (i80 < 0) { std::cout << "WARNING: t80 non trovato!" << std::endl; continue; }

                for (int i = i80 - 1; i > 0; --i) {
                    if (avgWF.amplitudes[i-1] < level20 && avgWF.amplitudes[i] >= level20) {
                        t20 = avgWF.times[i-1] + (level20 - avgWF.amplitudes[i-1]) *
                              (avgWF.times[i] - avgWF.times[i-1]) / (avgWF.amplitudes[i] - avgWF.amplitudes[i-1]);
                        i20 = i;
                        break;
                    }
                }
                if (i20 < 0) { std::cout << "WARNING: t20 non trovato!" << std::endl; continue; }

                // ----------------------------
                // Salva waveform media
                // ----------------------------
                dir_avgWF->cd();
                TString avgWF_name = TString::Format("%s-%s_thr%.1f-%.1f_avgWF",
                                                     vbias.Data(), ch_name.Data(), low_mult, up_mult);

                TGraphErrors* g_avg = new TGraphErrors(avgWF.times.size(),
                                                       avgWF.times.data(),
                                                       avgWF.amplitudes.data(),
                                                       nullptr,
                                                       avgWF.errors.empty() ? nullptr : avgWF.errors.data());
                g_avg->SetName(avgWF_name);
                g_avg->SetTitle(TString::Format("Average WF %s %s thr %.1f-%.1f;time (ns);amplitude (mV)",
                                                vbias.Data(), ch_name.Data(), low_mult, up_mult));

                // ----------------------------
                // Fit lineare 20%-80%
                // ----------------------------
                TF1 fitFunc("fitFunc", "[0]*x + [1]", t20, t80);
                fitFunc.SetParameters((level80 - level20)/(t80 - t20), 0);

                g_avg->Fit(&fitFunc, "RQ");
                double fit_slope     = fitFunc.GetParameter(0);
                double fit_slope_err = fitFunc.GetParError(0);

                std::cout << "Slope 20%-80% WF media: " << fit_slope << " ± " << fit_slope_err << std::endl;

                double ratio = baseline_sigma / fit_slope;
                double ratio_err = ratio * std::sqrt(std::pow(fit_slope_err/fit_slope,2) +
                                                     std::pow(baseline_sigma_err/baseline_sigma,2));
                std::cout << Form("   -> Slope / baseline_sigma = %.3f ± %.3f", ratio, ratio_err) << std::endl;

                ++pe_idx;
                tree_pe = pe_idx;
                tree_slope = fit_slope;
                tree_slope_err = fit_slope_err;
                tree_baseline = baseline_sigma;
                tree_baseline_err = baseline_sigma_err;
                tree_jitter = ratio;
                tree_jitter_err = ratio_err;

                dir_avgWF->cd();
                g_avg->Write();
                info_tree->Fill();

                // ----------------------------
                // Istogramma 2D tempo vs ampiezza
                // ----------------------------
                // ----------------------------
                // Istogrammi 2D
                // ----------------------------
                
                if (!selected.empty()) {
                    
                    // --- Primo 2D: tempo vs ampiezza ---
                    TH2D* h_wf2D_time_amp = new TH2D(Form("h2D_time_amp_%s_%s_thr%.1f-%.1f",
                                                         vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                                     Form("WF density %s %s thr %.1f-%.1f;time (ns);amplitude (mV)",
                                                          vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                                     1000.0, 30, 80,
                                                     1000, -10, 90);

                    // --- Secondo 2D: crossing time vs max amp ---
                    TH2D* h_wf2D_cross_amp = new TH2D(Form("h2D_cross_amp_%s_%s_thr%.1f-%.1f",
                                                          vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                                      Form("Crossing vs MaxAmp %s %s thr %.1f-%.1f;crossing time (ns);max amplitude (mV)",
                                                           vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                                      1000.0, 0, 100,
                                                      4000.0, 0, 200);

                    // Primo 2D: tutti i punti dei waveform
                    for (const auto& wf : selected) {
                        for (size_t i = 0; i < wf.times.size(); ++i) {
                            h_wf2D_time_amp->Fill(wf.times[i], wf.amplitudes[i]);
                        }
                    }

                    // Secondo 2D: crossing time vs max amplitude
                    for (const auto& wf : selected) {
                        h_wf2D_cross_amp->Fill(wf.max_amplitude, wf.crossing_time);
                    }
                    
                        // --- Istogramma per l'integrale tra 37 e 57 ns ---
                        TH1D* h_wf_integral = new TH1D(
                            Form("h_wf_integral_%s_%s_thr%.1f-%.1f", vbias.Data(), ch_name.Data(), low_mult, up_mult),
                            Form("Integral WF 37-57 ns %s %s thr %.1f-%.1f;integral (mV*ns);counts",
                                 vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                                       1000.0, 0, 500  // Adatta i bin e il range a seconda dei tuoi segnali
                        );

                        for (const auto& wf : selected) {
                            double integral = 0.0;
                            const auto& times = wf.times;
                            const auto& amps  = wf.amplitudes;

                            for (size_t i = 1; i < times.size(); ++i) {
                                // Considera solo il range 37-57 ns
                                if (times[i] < 37.0) continue;
                                if (times[i] > 57.0) break;

                                // Trapezio tra due punti consecutivi
                                double dt = times[i] - times[i-1];
                                double avg_amp = 0.5*(amps[i] + amps[i-1]);
                                integral += avg_amp * dt;
                            }

                            h_wf_integral->Fill(integral);
                        }

                    // Proiezione asse X del secondo 2D (crossing vs max amp)
                    TH1D* h_projX_cross_amp = h_wf2D_cross_amp->ProjectionX(
                        Form("h_projX_cross_amp_%s_%s_thr%.1f-%.1f",
                             vbias.Data(), ch_name.Data(), low_mult, up_mult)
                    );

                    h_projX_cross_amp->SetTitle(
                        Form("Proiezione asse X: max amp %s %s thr %.1f-%.1f",
                             vbias.Data(), ch_name.Data(), low_mult, up_mult)
                    );
                    h_projX_cross_amp->GetXaxis()->SetTitle("amplitude (mV)");
                    h_projX_cross_amp->GetYaxis()->SetTitle("counts");
                    // --- CFD al 50% ---
                    // Creiamo un istogramma per i tempi CFD
                    TH1D* h_wf_cfd = new TH1D(
                        Form("h_wf_cfd_%s_%s_thr%.1f-%.1f", vbias.Data(), ch_name.Data(), low_mult, up_mult),
                        Form("CFD 50%% %s %s thr %.1f-%.1f;time (ns);counts",
                             vbias.Data(), ch_name.Data(), low_mult, up_mult),
                                              400.0, 30, 50 // range e bin da adattare a seconda del segnale
                    );

                    for (const auto& wf : selected) {
                        double thresh50 = 0.5 * wf.max_amplitude; // soglia 50%
                        const auto& times = wf.times;
                        const auto& amps  = wf.amplitudes;

                        // Trova l'indice del massimo entro 42 ns
                        size_t max_index = 0;
                        double local_max = -1e9;
                        for (size_t j = 0; j < times.size(); ++j) {
                            if (times[j] > 42.0) break;   // stop oltre 42 ns
                            if (amps[j] > local_max) {
                                local_max = amps[j];
                                max_index = j;
                            }
                        }

                        // Trova il punto in cui si passa sopra la soglia andando all'indietro
                        for (size_t i = max_index; i > 0; --i) {
                            if (amps[i] >= thresh50 && amps[i-1] < thresh50) {
                                double t_cfd = times[i-1] + (thresh50 - amps[i-1]) * (times[i] - times[i-1]) / (amps[i] - amps[i-1]);
                                h_wf_cfd->Fill(t_cfd);
                                break;
                            }
                        }
                    }

                    // --- Fit con Q-Gauss ---
                    TF1* fitFunc = GaussStuff::FitQGaussPlusConst(h_wf_cfd, Form("fit_%dpe", int(low_mult)));  // o un indice logico

                   

                    // Proiezione asse Y del secondo 2D (crossing time vs max amp)
                    TH1D* h_projY_cross_amp = h_wf2D_cross_amp->ProjectionY(
                        Form("h_projY_cross_amp_%s_%s_thr%.1f-%.1f",
                             vbias.Data(), ch_name.Data(), low_mult, up_mult)
                    );

                    h_projY_cross_amp->SetTitle(
                        Form("Proiezione asse Y: crossing time %s %s thr %.1f-%.1f",
                             vbias.Data(), ch_name.Data(), low_mult, up_mult)
                    );
                    h_projY_cross_amp->GetXaxis()->SetTitle("crossing time (ns)");
                    h_projY_cross_amp->GetYaxis()->SetTitle("counts");

                    // --- Fit con Q-Gauss ---
                    TF1* fitFuncY = GaussStuff::FitQGaussPlusConst(h_projY_cross_amp, Form("fitY_%dpe", int(low_mult)));

                    TH1D* h_cross_fixed = new TH1D(
                            Form("h_cross_fixed_%s_%s_thr%.1f-%.1f", vbias.Data(), ch_name.Data(), low_mult, up_mult),
                            Form("Crossing time fixed %s %s thr %.1f-%.1f;time (ns);counts", vbias.Data(), ch_name.Data(), low_mult, up_mult),
                            400, 30, 50 // range e bin da adattare
                        );

                    for (const auto& wf : selected) {
                        double t_cross = TimeStuff::computeCrossingTime(wf.times.data(), wf.amplitudes.data(), wf.times.size(), fixed_threshold);
                        if (t_cross > 0) h_cross_fixed->Fill(t_cross);
                    }

                        // Fit con Q-Gauss
                    TF1* fit_cross = GaussStuff::FitQGaussPlusConst(h_cross_fixed, Form("fit_cross_%dpe", int(low_mult)));


                    // Salvataggio
                    
                    dir_histowf->cd();
                    h_wf2D_time_amp->Write();
                    h_wf2D_cross_amp->Write();
                    h_projX_cross_amp->Write();
                    h_wf_integral->Write();
                    h_wf_cfd->Write();
                    h_projY_cross_amp->Write();
                    h_cross_fixed->Write();
                }
            }
        }

        outFile->cd();
        info_tree->Write();
        outFile->Close();
    }
}




/*
 std::map<TString, std::vector<double>> tau_means_map;  // media tau per ciascun Vbias
 std::map<TString, std::vector<double>> tau_errors_map; // errore tau (sigma gaussiana)
 std::map<TString, std::vector<double>> vbias_map;      // vbias corrispondente
 
 
 for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
     // selezione wf come prima
     std::vector<double> tau_selected, amp_selected;
     
     // Istogrammi esistenti
     TH2D* h2_amp_tau      = new TH2D(TString::Format("%s-%s_tau-voltage", vbias.Data(), ch_name.Data()),
                                      "Amplitude vs Tau;amplitude (mV);#tau (ns)", 250.0, 0, 50, 200, 0, 100);
     TH2D* h2_slope_amp    = new TH2D(TString::Format("%s-%s_slope-voltage", vbias.Data(), ch_name.Data()),
                                      "Slope vs Amp;Amplitude (mV);Slope (mV/ns)", 1000, 0, 500, 500, 0, 50);
     TH2D* h2_integral_amp = new TH2D(TString::Format("%s-%s_integral-voltage", vbias.Data(), ch_name.Data()),
                                      "Integral vs Amp;Amplitude (mV);Integral (mV×ns)", 1000, 0, 500, 2000, 0, 2000);
     TH1D* h_amp_dist      = new TH1D(TString::Format("%s-%s_amplitude", vbias.Data(), ch_name.Data()),
                                      "Amplitude Distribution <35ns", 100.0, -10, 10);

     TH2D* h2_tmax_tau     = new TH2D(TString::Format("%s-%s_tmax-tau", vbias.Data(), ch_name.Data()),
                                      "t_{max} vs Tau; t_{max} (ns); Tau (ns)", 4000.0, 0, 200, 200, 0, 100);
     TH2D* h2_tmax_slope   = new TH2D(TString::Format("%s-%s_tmax-slope", vbias.Data(), ch_name.Data()),
                                      "t_{max} vs Slope; t_{max} (ns); Slope (mV/ns)", 4000.0, 0, 200, 500, 0, 50);
     TH2D* h2_tmax_integral = new TH2D(TString::Format("%s-%s_tmax-integral", vbias.Data(), ch_name.Data()),
                                       "t_{max} vs Integral; t_{max} (ns); Integral (mV×ns)", 4000.0, 0, 200, 2000, 0, 2000);

     
     TH2D* h2_amp_t_max = new TH2D(
         TString::Format("%s-%s_t-voltage", vbias.Data(), ch_name.Data()),
         "Amplitude vs Tau;Amplitude (mV);t_{max}(ns)",
         4000, 0.0, 200.0, 500,   0.0, 100.0);


     
     for (const auto& wf : selected) {
         const auto& times = wf.first;
         const auto& amps  = wf.second;
         if (amps.empty()) continue;
         
         // Tau & idx_max
         auto [tau, tau_err, idx_max] = Fitting::estimateTau(times, amps);
         double maxAmp = amps[idx_max];
         double t_max  = times[idx_max];

         tau_selected.push_back(tau);
         amp_selected.push_back(maxAmp);
         h2_amp_tau->Fill(maxAmp, tau);
         h2_tmax_tau->Fill(t_max, tau);
         h2_amp_t_max->Fill(maxAmp, t_max);

         // 🔹 Metodo 1: slope t80-t20
         auto [slope, slope_err, t0] = Fitting::computeSlopeT80T20(times, amps, idx_max);
         h2_slope_amp->Fill(maxAmp, slope);
         h2_tmax_slope->Fill(t_max, slope);

         // 🔹 Metodo 2: integrale da t0 a 1% maxAmp
         double integral = Waveform::computeIntegralFromT0Minus5ToEnd(times, amps, t0);
         h2_integral_amp->Fill(maxAmp, integral);
         h2_tmax_integral->Fill(t_max, integral);

         // 🔹 Metodo 3: distribuzione ampiezze <35ns
         Waveform::fillAmplitudeDistributionBefore35ns(times, amps, h_amp_dist);
     }

     TH1D* hProj_tau = Utils::projectionY(h2_amp_tau, lower_thr, upper_thr, TString::Format("%s-%s_tau-histo", vbias.Data(), ch_name.Data()));
     
     auto [tau_val, tau_err, tau_sigma, tau_sigma_err, tau_fit] = GaussStuff::q_gauss_projection(hProj_tau);
     
     TH1D* hProj_slope = Utils::projectionY(h2_slope_amp, lower_thr, upper_thr, TString::Format("%s-%s_slope-histo", vbias.Data(), ch_name.Data()));
     auto [slope_val, slope_err, slope_sigma, slope_sigma_err, slope_fit] = GaussStuff::q_gauss_projection(hProj_slope);
     
     TH1D* hProj_integral = Utils::projectionY(h2_integral_amp, lower_thr, upper_thr, TString::Format("%s-%s_integral-histo", vbias.Data(), ch_name.Data()));
     auto [integral_val, integral_err, integral_sigma, integral_sigma_err, integral_fit] = GaussStuff::gauss_fit_projection(hProj_integral);
     
     TH1D* hProj_amp = Utils::projectionX(h2_amp_tau, 0, 1000, TString::Format("%s-%s_amp-histo", vbias.Data(), ch_name.Data()));
     auto [amp_val, amp_err, amp_sigma, amp_sigma_err, amp_fit] = GaussStuff::gauss_fit_projection(hProj_amp);
     
     
     auto [baseline_val, baseline_err, baseline_sigma, baseline_sigma_err, baseline_fit] = GaussStuff::gauss_fit_projection(h_amp_dist);
     
     
     // 1️⃣ Creazione della canvas
     TCanvas* c1 = new TCanvas("c1", "Average Waveform", 800, 600);
     c1->cd(); // attiva la canvas corrente

   
     
     auto avgWF = Waveform::computeAverageWaveform(selected);
     dir_avgWF->cd();
     TString avgWF_name = TString::Format("%s-%s_avgWF", vbias.Data(), ch_name.Data());
     TGraphErrors* g_avg = new TGraphErrors(avgWF.times.size(), avgWF.times.data(), avgWF.amplitudes.data(), nullptr, avgWF.errors.data());
     g_avg->SetName(avgWF_name);
     g_avg->SetTitle(TString::Format("Average Waveform corrected baseline %s %s;time (ns);amplitude (mV)", vbias.Data(), ch_name.Data()));
     if (!g_avg) return;
     g_avg->Write();
     // 2️⃣ Disegna il TGraphErrors
     g_avg->SetMarkerStyle(20);
     g_avg->SetMarkerColor(kMagenta);
     g_avg->Draw("AP"); // A = draw axes, P = points

     
     // Costanti
     constexpr double R_ohm = 1000.0; // Ohm
     constexpr double q_e = 1.602176634e-19; // C
     constexpr double gain_dB = 20.0; // <-- valore reale in dB
     
     // Conversione da dB a rapporto di tensione
     const double voltage_real_gain = std::pow(10.0, gain_dB / 20.0);

     
     // Gain calcolato (numero di elettroni)
     double gain_val = (integral_val * 1e-12) / (R_ohm * q_e * voltage_real_gain);
     double gain_err = (integral_err * 1e-12) / (R_ohm * q_e * voltage_real_gain); // solo errore dell'integrale

     
     // SNR
     double snr_val = amp_val / baseline_sigma;
     double snr_err = std::sqrt(
         std::pow(amp_err / baseline_sigma, 2) +
         std::pow((-amp_val * baseline_sigma_err) / (baseline_sigma * baseline_sigma), 2)
     );
     
     std::cout << lower_thr << " " << upper_thr << " " << amp_val << " " << amp_err << " " << baseline_sigma << " " << baseline_sigma_err << "\n";
     // Scomposizione gain in mantissa * 10^exponent
     int gain_exp = (gain_val != 0.0) ? static_cast<int>(std::floor(std::log10(std::fabs(gain_val)))) : 0;
     double gain_mantissa = gain_val / std::pow(10.0, gain_exp);

     // Stessa cosa per l'errore (scalato allo stesso esponente)
     double gain_err_mantissa = gain_err / std::pow(10.0, gain_exp);

     // Testo TPaveText
     TString infoText;
     // ====== Calcolo slope T20-T80 dalla wf media ======
     
     // Trova max della wf media
     auto it_max = std::max_element(avgWF.amplitudes.begin(), avgWF.amplitudes.end());
     int idx_max_avg = std::distance(avgWF.amplitudes.begin(), it_max);
     double amp_max = avgWF.amplitudes[idx_max_avg];
     double t_max   = avgWF.times[idx_max_avg];

     // ====== Calcolo slope T20-T80 dalla wf media ======
     auto [slope_avg, slope_avg_err, t0_avg] = Fitting::computeSlopeT80T20(
         avgWF.times, avgWF.amplitudes, idx_max_avg
     );

     auto [tau_avg, tau_err_avg, idx_max_avg_tau] = Fitting::estimateTau(avgWF.times, avgWF.amplitudes);
     
     
     // definisci il fit esponenziale, shiftato a partire da t_max
     TF1* expDecayFit = new TF1("expDecayFit",
                                "[0]*exp(-(x-[2])/[1])",
                                t_max, avgWF.times.back());
     // Dai nomi ai parametri
     expDecayFit->SetParNames("amp_max", "#tau", "t_{max}");

     // usa i parametri stimati dal fit
     expDecayFit->SetParameter(0, amp_max); // ampiezza al massimo
     expDecayFit->SetParameter(1, tau_avg);     // tau medio dal fit
     expDecayFit->SetParameter(2, t_max);   // tempo del massimo (shift)

     // un po’ di estetica
     expDecayFit->SetLineColor(kBlack);
     expDecayFit->SetLineStyle(1);
     expDecayFit->Draw("same");
     
     
     // Livelli 20% e 80% dell'ampiezza massima
     double level_20 = 0.2 * amp_max;
     double level_80 = 0.8 * amp_max;

     // Trova i tempi corrispondenti ai livelli 20% e 80%
     double t_20 = 0.0, t_80 = 0.0;
     for (size_t i = 0; i < avgWF.amplitudes.size(); ++i) {
         if (avgWF.amplitudes[i] >= level_20 && t_20 == 0.0) t_20 = avgWF.times[i];
         if (avgWF.amplitudes[i] >= level_80) { t_80 = avgWF.times[i]; break; }
     }

     // Fit lineare tra t_20 e t_80
     TF1* fit20_80 = new TF1("fit20_80", "[0]*(x-[1])", t_20, t_80);
     fit20_80->SetParameter(0, slope_avg); // punto di partenza con slope media
     fit20_80->SetParameter(1, t0_avg);    // t0 iniziale (origine linea)
     fit20_80->SetLineColor(kBlue);
     fit20_80->SetLineStyle(1);
     fit20_80->Draw("same");
     
     // 2️⃣ Linea estesa oltre t_20-t_80 fino alla fine della waveform
     TF1* slopeLineExtended = new TF1("slopeLineExtended", "[0]*(x-[1])", 0, 50);
     slopeLineExtended->SetParameter(0, slope_avg);
     slopeLineExtended->SetParameter(1, t0_avg);
     slopeLineExtended->SetLineColor(kBlue);
     slopeLineExtended->SetLineStyle(1); // tratteggiata per distinguere
     slopeLineExtended->Draw("same");


     // Disegna linee orizzontali 20% e 80% per riferimento
     TLine* line20 = new TLine(0, level_20, 40, level_20);
     line20->SetLineColor(kRed);
     line20->SetLineStyle(2);
     line20->Draw("same");

     TLine* line80 = new TLine(0, level_80, 40, level_80);
     line80->SetLineColor(kRed);
     line80->SetLineStyle(2);
     line80->Draw("same");

     TLegend* legend = new TLegend(0.55, 0.60, 0.88, 0.88); // (x1,y1,x2,y2) in NDC
     legend->SetBorderSize(0);
     legend->SetFillStyle(0);
     legend->SetTextSize(0.035);

     legend->AddEntry(expDecayFit, "A e^{-(x - t_{max}) / #tau}", "l");
     legend->AddEntry(fit20_80,    "m (x - t_{0})", "l");
     legend->AddEntry(line20,      "thr @20%", "l");
     legend->AddEntry(line80,      "thr @80%", "l");

     legend->Draw();
     
     // Aggiorna canvas
     c1->Update();

     // Jitter elettronica
     double jitter_val = baseline_sigma / slope_avg;
     double jitter_err = std::sqrt(
         std::pow(baseline_sigma_err / slope_avg, 2) +
         std::pow((-baseline_sigma * slope_avg_err) / (slope_avg * slope_avg), 2)
     );
     std::cout << jitter_val << " " << jitter_err << "\n";
     // Aggiungi questa slope al TPaveText
     /*infoText.Form(
         "#tau = %.2f #pm %.2f ns\n"
         //"S_{fit} = %.2f #pm %.2f mV/ns\n"
         "S_{avgWF} = %.2f #pm %.2f mV/ns\n"
         //"G = %.3f #pm %.3f #times 10^{%d}\n"
         "#sigma_{jitter} = %.2f #pm %.2f ps\n"
         "SNR = %.1f #pm %.1f",
         tau_val, tau_err,
         slope_val, slope_err,            // slope dal fit gaussiano
         slope_avg, slope_avg_err,        // slope dalla wf media
         //gain_mantissa, gain_err_mantissa, gain_exp,
         jitter_val*1000, jitter_err*1000,
         snr_val, snr_err
     );*/
     
     /*infoText.Form(
         "#tau = %.2f #pm %.2f ns\n"
         "slew-rate = %.2f #pm %.2f mV/ns\n"
         "#sigma_{jitter} = %.2f #pm %.2f ps\n",
         tau_val, tau_err,
         slope_avg, slope_avg_err,
         jitter_val*1000, jitter_err*1000
     );*/

     
     
     
     /*infoText.Form(
         "#tau = %.2f #pm %.2f ns\n"
         "S_{avgWF} = %.2f #pm %.2f mV/ns\n"
         "#sigma_{jitter} = %.2f #pm %.2f ps",
         tau_val, tau_err,
         slope_avg, slope_avg_err,        // slope dalla wf media
         jitter_val*1000, jitter_err*1000
     );*/



     // Creazione TPaveText
     /*TPaveText* paveInfo = Graphics::CreateInfoBox(0.6, 0.6, 0.9, 0.9, infoText, kWhite, kBlack, 0.03, 1);
     paveInfo->Draw();

     
     // 6️⃣ Scrivi TGraph
     c1->Write();
     
     dir_tau_histo->cd();
     // Istogrammi esistenti
     h2_amp_tau->Write();
     h2_tmax_tau->Write();
     h2_amp_t_max->Write();
     hProj_tau->Write();
     
     dir_slope_histo->cd();
     h2_slope_amp->Write();
     h2_tmax_slope->Write();
     hProj_slope->Write();
     
     dir_integral_histo->cd();
     h2_integral_amp->Write();
     h2_tmax_integral->Write();
     hProj_integral->Write();
     
     dir_baseline_histo->cd();
     h_amp_dist->Write();
     hProj_amp->Write();
     
     double t_min = 35.0;
     double t_max_int = 160.0;
     double integral = 0.0;
     double integral_err_avgwf = 0.0;

     for (int i = 0; i < g_avg->GetN() - 1; ++i) {
         double x1 = g_avg->GetX()[i];
         double x2 = g_avg->GetX()[i+1];
         double y1 = g_avg->GetY()[i];
         double y2 = g_avg->GetY()[i+1];
         double e1 = g_avg->GetEY()[i];
         double e2 = g_avg->GetEY()[i+1];

         // considera solo intervallo richiesto
         if (x2 < t_min) continue;
         if (x1 > t_max_int) break;

         // taglia i trapezi ai bordi del range
         double dx = x2 - x1;
         if (x1 < t_min) dx = x2 - t_min;
         if (x2 > t_max_int) dx = t_max_int - x1;

         // integrale con trapezio
         double dy = 0.5 * (y1 + y2);
         integral += dy * dx;

         // errore sul trapezio
         double trapezoid_err = std::sqrt( std::pow(dx/2.0 * e1, 2) + std::pow(dx/2.0 * e2, 2) );
         integral_err_avgwf += trapezoid_err * trapezoid_err; // somma quadratica
     }

     // radice quadrata finale per errore totale
     integral_err = std::sqrt(integral_err);

     std::cout << "Integrale da " << t_min << " ns a " << t_max_int
               << " ns = " << integral << " ± " << integral_err << std::endl;
     
     // Assumendo che sei dentro il ciclo per ogni canale
     TDirectory* dirWF = gDirectory->mkdir("wf-all");
     dirWF->cd();

     int wf_index = 0;  // contatore per dare nomi univoci
     for (const auto& wf : selected) {
         const auto& times = wf.first;
         const auto& amps  = wf.second;
         if (amps.empty()) continue;

         TString gname = TString::Format("wf_%03d_%s", wf_index, ch_name.Data());
         TGraph* g_wf = new TGraph(times.size(), times.data(), amps.data());
         g_wf->SetTitle(TString::Format("Waveform %03d - %s;Time (ns);Amplitude (mV)", wf_index, ch_name.Data()));
         g_wf->SetLineColor(kBlue);
         g_wf->SetLineWidth(1);

         g_wf->Write(gname);  // salva il TGraph nella directory wf-all
         wf_index++;
     }

     // Torna alla directory principale se vuoi salvare altri oggetti
     gDirectory->cd();*/
