#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TList.h>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>

#include "root-library/utils.h"
#include "root-library/waveform.h"
#include "root-library/gauss_stuff.h"
#include "root-library/timestuff.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TApplication.h>

// =====================================================
// Funzione principale
// =====================================================
void timing(const TString& filename, const TString& out_root_name)
{
    // =====================================================
    //                PARAMETRI MODIFICABILI
    // =====================================================
    const double t_min = 40.0;
    const double t_max = 42.0;

    const int    n_bins_max_amp      = 1000.0;   // numero di bin
    const double lower_range_max_amp = -200.0;    // valore minimo asse x
    const double upper_range_max_amp = 200.0;  // valore massimo asse x
    const int    n_bins_delta_t      = 10000.0;   // numero di bin
    const double lower_range_delta_t = -500;    // valore minimo asse x
    const double upper_range_delta_t = 500;  // valore massimo asse x
    const int n_peaks_requested = 4; // ad esempio
    constexpr bool use_linear_fit = true;  // oppure false
    const double laser_threshold = 70.0;
    // =====================================================

    // Apri file di input
    auto inFile = std::unique_ptr<TFile>(TFile::Open(filename));
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: cannot open input file " << filename << std::endl;
        return;
    }

    // Prepara laser tree & branches
    TTree* tree_laser = dynamic_cast<TTree*>(inFile->Get("laser"));
    std::vector<std::pair<std::vector<double>, std::vector<double>>> laser_waveforms;
    Double_t *laser_amp = nullptr, *laser_time = nullptr;
    if (tree_laser) {
        laser_amp  = Utils::setupBranch_dgz(tree_laser, "amplitude");
        laser_time = Utils::setupBranch_dgz(tree_laser, "time");

        // Numero di eventi laser
        Long64_t nEntries = tree_laser->GetEntries();

        for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
            tree_laser->GetEntry(iEntry);

            std::vector<double> time_vec(1024);
            std::vector<double> amp_vec(1024);

            for (int i = 0; i < 1024; ++i) {
                time_vec[i] = laser_time[i];
                amp_vec[i] = laser_amp[i];
            }

            laser_waveforms.emplace_back(std::make_pair(std::move(time_vec), std::move(amp_vec)));
        }

        std::cout << "Loaded " << laser_waveforms.size() << " laser waveforms\n";
    }


    // Carica tutti i canali
    std::map<TString, TTree*>           tree_map;
    std::map<TString, Double_t*>        time_map;
    std::map<TString, Double_t*>        amp_map_raw;

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
                    } else {
                        std::cerr << "Failed to setup branches for " << name << std::endl;
                    }
                }
            }
        }
    }

    // Allinea e corregge baseline
    Waveform waveAligner(tree_laser, laser_time, laser_amp);
    std::map<TString, TF1*> linearFitsMap;

    std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_aligned_waveforms;
    std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_baseline_corrected_waveforms;

    for (const auto& [ch_name, ch_tree] : tree_map) {
        auto ch_time = time_map[ch_name];
        auto ch_amp  = amp_map_raw[ch_name];

        if (ch_tree && ch_time && ch_amp) {
            auto aligned = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, true);
            auto corrected_aligned = waveAligner.correctBaseline(aligned);
            all_aligned_waveforms[ch_name] = std::move(corrected_aligned);
            
            // Ora prendi direttamente il baseline corrected originale non allineato:
            auto raw = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, false);
            auto corrected_raw = waveAligner.correctBaseline(raw);
            all_baseline_corrected_waveforms[ch_name] = std::move(corrected_raw);

            std::cout << "Aligned and baseline-corrected channel " << ch_name
                      << " with " << all_aligned_waveforms[ch_name].size() << " waveforms\n";
        } else {
            std::cerr << "Skipping channel " << ch_name << " due to missing data pointers\n";
        }
    }

    // File di output
    auto outFile = std::unique_ptr<TFile>(TFile::Open(out_root_name, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: cannot create output file " << out_root_name << std::endl;
        return;
    }
    TDirectory* dir_amp = outFile->mkdir("amp_max-distribution");
    TDirectory* dir_cal = outFile->mkdir("calibration");

    /** Fare il lavoro del fit multi gauss a mean fissa e sigma fissa**/
    
    // Calcolo e riempimento istogrammi usando la funzione
    for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
        TString unique_hist_name = TString::Format("%s", ch_name.Data());
        
        
        auto [hist_noise, hist_signal, hist_diff] = waveAligner.calibration(unique_hist_name, waveforms, t_min, t_max, n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        if (!hist_diff) continue;
        
        
        
        dir_amp->cd();
        TString fit_name = TString::Format("%s-multi-gauss", ch_name.Data());
        TF1* fit = GaussStuff::multi_gauss(hist_signal.get(), n_peaks_requested, fit_name);
        hist_signal->Write();
    

        if (!fit) continue;

        int n_par = fit->GetNpar();
        int n_peaks = n_par / 3;

        // vettori per punti e errori
        std::vector<double> x_vals, y_vals, y_errs;

        for (int i = 0; i < n_peaks; ++i) {
            double mean = fit->GetParameter(i * 3 + 1);
            double mean_err = fit->GetParError(i * 3 + 1);
            int index = i; // photo-electron number

            x_vals.push_back(index);
            y_vals.push_back(mean);
            y_errs.push_back(mean_err);
        }

        // Crea TGraphErrors
        TGraphErrors* graph = new TGraphErrors(n_peaks, x_vals.data(), y_vals.data(), nullptr, y_errs.data());
        graph->SetName((ch_name + "_mean_vs_pe").Data());
        graph->SetTitle("Mean vs Photo-electron number;#p.e.;amplitude (mV)");
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlue);

        dir_cal->cd();
        graph->Write();

        if (use_linear_fit) {
            // Fit lineare pol1
            TF1* linearFit = new TF1((ch_name + "_linearFit").Data(), "pol1", 1, n_peaks+1);
            graph->Fit(linearFit, "RQ"); // fit silenzioso
            linearFitsMap[ch_name] = linearFit;
            linearFit->Write();
        } if (!use_linear_fit) {
            // Prendi direttamente la media della prima gaussiana dal fit multigaussiano
            double first_mean = fit->GetParameter(4);       // parametro 1 = media della prima gaussiana
            double first_mean_err = fit->GetParError(4);   // errore della media della prima gaussiana

            
            TF1* constFit = new TF1((ch_name + "_constantFit").Data(),
                                    [first_mean](double*, double*) { return first_mean; },
                                    0, n_peaks+1, 0);
            constFit->SetLineColor(kRed);
            constFit->Write();
            linearFitsMap[ch_name] = constFit; // riutilizzo della mappa per uniformità
        }


        dir_amp->cd();
    }


    
   /// Crea directory per gli istogrammi 2D
    TDirectory* dir_fixed_th = outFile->mkdir("fixed_threshold");
    TDirectory* dir_wf = outFile->mkdir("waveform");
    TDirectory* dir_wf_sig = outFile->mkdir("waveform_clean");
    TDirectory* dir_cfd= outFile->mkdir("cfd");
    
    
    // Mappa calib_params da slope e intercept ricavati dai fit lineari o costanti
    std::map<TString, std::pair<double,double>> calib_params;
    for (const auto& [ch_name, fit] : linearFitsMap) {
        double slope = 0.0;
        double intercept = 0.0;

        if (use_linear_fit) {
            slope = fit->GetParameter(1);       // coefficiente lineare
            intercept = fit->GetParameter(0);   // intercetta
        } else {
            slope = fit->Eval(0);   // valore della costante (prima media del fit multigaussiano)
            intercept = 0.0;        // intercetta = 0
        }

        calib_params[ch_name] = std::make_pair(slope, intercept);
    }

    // Loop sui canali
    for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
        if (calib_params.find(ch_name) == calib_params.end()) {
            std::cerr << "No calibration params for " << ch_name << ", skipping.\n";
            continue;
        }

        auto [slope, intercept] = calib_params[ch_name];

        // soglia per l'istogramma delta t
        double threshold_fixed = 0.5 * slope + intercept;

        // crea istogramma delta t
        TString hist_name_fixed  = ch_name + "_fixed_th_diff_laser";
        TString hist_title_fixed = hist_name_fixed + ";amplitude (mV);#Delta t (t_{ch} - t_{laser}) (ns)";
        auto h_fixed_th_diff = std::make_unique<TH2F>(hist_name_fixed, hist_title_fixed,
                                                      n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
                                                      n_bins_delta_t, lower_range_delta_t, upper_range_delta_t);

        // soglie per classificazione waveform
        std::vector<double> thresholds = {
            0.5 * slope + intercept,
            1.5 * slope + intercept,
            2.5 * slope + intercept,
            3.5 * slope + intercept
        };

        // creazione istogrammi 2D per ogni classe
        std::vector<std::unique_ptr<TH2F>> h_class;
        for (size_t cls = 0; cls < thresholds.size() - 1; ++cls) {
            TString hist_name  = TString::Format("%s_class_%zu", ch_name.Data(), cls+1);
            TString hist_title = hist_name + ";time (ns);amplitude (mV)";
            h_class.push_back(std::make_unique<TH2F>(
                hist_name, hist_title, 1250, -50, 200,3000.0, -50, 250));
        }
        
        std::vector<std::unique_ptr<TH2F>> h_class_tc38_42;
        for (size_t cls = 0; cls < thresholds.size() - 1; ++cls) {
            TString hist_name  = TString::Format("%s_class_%zu_Tcrossing_38_42", ch_name.Data(), cls+1);
            TString hist_title = hist_name + ";time (ns);amplitude (mV)";
            h_class_tc38_42.push_back(std::make_unique<TH2F>(
                hist_name, hist_title, 1250, -50, 200, 3000, -50, 250));
        }

        size_t n_events = std::min(waveforms.size(), laser_waveforms.size());

        for (size_t i = 0; i < n_events; ++i) {
            const auto& [times_vec, amps_vec] = waveforms[i];
            const auto& [laser_times, laser_amps] = laser_waveforms[i];

            int size       = static_cast<int>(times_vec.size());
            int size_laser = static_cast<int>(laser_times.size());

            if (size < 2 || size_laser < 2) continue;

            double crossing_time_laser =
                TimeStuff::computeCrossingTime(laser_times.data(), laser_amps.data(), size_laser, laser_threshold);

            if (crossing_time_laser < 0) continue;

            // istogramma delta t
            auto res = TimeStuff::analyzeWaveform(times_vec.data(), amps_vec.data(), size, threshold_fixed);
            if (res.crossing_time >= 0) {
                double delta_t = res.crossing_time;
                h_fixed_th_diff->Fill(amps_vec[res.max_index], delta_t);
            }

            // classificazione waveform
            if (res.crossing_time < 0) continue;

            double max_amp = amps_vec[res.max_index];
            int cls_index = -1;
            for (size_t cls = 0; cls < thresholds.size() - 1; ++cls) {
                if (max_amp >= thresholds[cls] && max_amp < thresholds[cls + 1]) {
                    cls_index = cls;
                    break;
                }
            }

            if (cls_index >= 0) {
                for (int j = 0; j < size; ++j) {
                    h_class[cls_index]->Fill(times_vec[j], amps_vec[j]);
                }
            }
            if (cls_index >= 0 && res.crossing_time >= 39.0 && res.crossing_time <= 41.0) {
                for (int j = 0; j < size; ++j) {
                    h_class_tc38_42[cls_index]->Fill(times_vec[j], amps_vec[j]);
                }
            }
        }

        // salva istogramma delta t
        dir_fixed_th->cd();
        h_fixed_th_diff->Write();

        // salva istogrammi waveform classificati
        dir_wf->cd();
        for (auto& hist : h_class) {
            hist->Write();
        }
        
        dir_wf_sig->cd();
        for (auto& hist : h_class_tc38_42) {
            hist->Write();
        }
    }

   
        
        
    outFile->cd();
    // Creo TTree "info" per memorizzare numero di waveforms analizzate
    TTree info("info", "Info about analysis");

    // Numero di waveforms totali (qui prendo il totale di waveforms analizzate, esempio primo canale)
    Int_t n_wf = 0;
    if (!all_aligned_waveforms.empty()) {
        n_wf = static_cast<Int_t>(all_aligned_waveforms.begin()->second.size());
    }

    // Creo ramo "n_waveforms"
    info.Branch("n_waveforms", &n_wf, "n_waveforms/I");

    // Riempio e salvo
    info.Fill();

    
    outFile->Write();
    outFile->Close();
}




/*
for (const auto& [ch_name, waveforms] : all_baseline_corrected_waveforms) {
    if (calib_params.find(ch_name) == calib_params.end()) {
        std::cerr << "No calibration params for " << ch_name << ", skipping.\n";
        continue;
    }

    auto [slope, intercept] = calib_params[ch_name];
    double threshold = 0.5 * slope + intercept;
    std::cout << "slope " << slope << " " << intercept << " " << threshold << " " << "\n";

    TString hist_name_fixed  = ch_name + "_fixed_th_diff_laser";
    TString hist_title_fixed = hist_name_fixed + ";amplitude (mV);#Delta t (t_{ch} - t_{laser}) (ns)";

    auto h_fixed_th_diff = std::make_unique<TH2F>(
        hist_name_fixed,
        hist_title_fixed,
        n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
        n_bins_delta_t, lower_range_delta_t, upper_range_delta_t
    );
    
    // ---- HISTO per CFD
        TString hist_name_cfd  = ch_name + "_cfd50_diff_laser";
        TString hist_title_cfd = hist_name_cfd + ";amplitude (mV);#Delta t (t_{ch} - t_{laser}) (ns)";

        auto h_cfd_diff = std::make_unique<TH2F>(
            hist_name_cfd,
            hist_title_cfd,
            n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
            n_bins_delta_t, lower_range_delta_t, upper_range_delta_t
        );

    size_t n_events = std::min(waveforms.size(), laser_waveforms.size());

    for (size_t i = 0; i < n_events; ++i) {
        const auto& [times_vec, amps_vec]   = waveforms[i];
        const auto& [laser_times, laser_amps] = laser_waveforms[i];

        int size       = static_cast<int>(times_vec.size());
        int size_laser = static_cast<int>(laser_times.size());

        if (size < 2 || size_laser < 2) continue;

        double crossing_time_laser =
            TimeStuff::computeCrossingTime(
                laser_times.data(),
                laser_amps.data(),
                size_laser,
                laser_threshold
            );

        if (crossing_time_laser < 0) continue;

        auto res = TimeStuff::analyzeWaveform(
            times_vec.data(),
            amps_vec.data(),
            size,
            threshold
        );

        
        if (res.max_index != -1 && res.crossing_time >= 0) {
            double max_amp = amps_vec[res.max_index];
            double delta_t = res.crossing_time - crossing_time_laser;
            h_fixed_th_diff->Fill(max_amp, delta_t);

            // --- CFD 50% andando a ritroso e interpolando ---
            double threshold_cfd = 0.5 * max_amp;
            double crossing_time_cfd = -1.0;

            for (int i = res.max_index; i > 0; --i) {
                double y1 = amps_vec[i-1];
                double y2 = amps_vec[i];

                if (y1 < threshold_cfd && y2 >= threshold_cfd) {
                    double x1 = times_vec[i-1];
                    double x2 = times_vec[i];

                    // interpolazione lineare
                    crossing_time_cfd = x1 + (threshold_cfd - y1) * (x2 - x1) / (y2 - y1);
                    break;
                }
            }

            if (crossing_time_cfd >= 0) {
                double delta_t_cfd = crossing_time_cfd - crossing_time_laser;
                h_cfd_diff->Fill(max_amp, delta_t_cfd);
            }
        }
    }

    // ---------- SCRITTURA ISTOGRAMMI ----------
    dir_fixed_th->cd();
    h_fixed_th_diff->Write();
    dir_cfd->cd();
    h_cfd_diff->Write();

}



// Directory "ct"

// Parametri della scansione
float thr_min  = 0;    // soglia minima
float thr_max  = 50;   // soglia massima
float thr_step = 0.25; // step
double t_min_scan = 37.0; // tempo minimo in ns
double t_max_scan = 47.0; // tempo massimo in ns

TDirectory* dir_ct = outFile->mkdir("thr-scan-38-42ns");

for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {
    if (waveforms.empty()) continue;

    std::vector<double> thr_values;
    std::vector<double> counts_all;
    std::vector<double> err_all;

    for (float thr = thr_min; thr <= thr_max; thr += thr_step) {
        int n_crossings_total = 0;

        for (const auto& wf : waveforms) {
            const auto& times = wf.first;
            const auto& amps  = wf.second;

            // Costruisci un TGraph temporaneo per il waveform
            TGraph gr(times.size(), times.data(), amps.data());

            // Ottieni i crossing usando la funzione aggiornata
            auto crossings = Waveform::get_transitions(gr, thr, +1.0f);

            // Conta solo i crossing entro l'intervallo 38-42 ns
            for (const auto& [t, y] : crossings) {
                if (t < t_min_scan) continue;
                if (t > t_max_scan) break;
                n_crossings_total++;
            }
        }

        thr_values.push_back(thr);
        counts_all.push_back(n_crossings_total);
        err_all.push_back(std::sqrt(n_crossings_total));
    }

    // Grafico totale con errori
    TString gname_all  = "g_threshold_scan_38_42ns_" + ch_name;
    TString gtitle_all = "Threshold scan 38-42 ns " + ch_name;
    TGraphErrors* g_scan_all = new TGraphErrors(thr_values.size(),
                                               thr_values.data(), counts_all.data(),
                                               nullptr, err_all.data());
    g_scan_all->SetName(gname_all);
    g_scan_all->SetTitle(gtitle_all);

    // Scrive il grafico
    dir_ct->cd();
    g_scan_all->Write();
}




size_t max_to_save = 50;
size_t saved_count = 0;

TDirectory* dir_wf = outFile->mkdir("waveforms"); // crea una directory per le waveform
dir_wf->cd();

for (const auto& [ch_name, waveforms] : all_baseline_corrected_waveforms) {
    if (calib_params.find(ch_name) == calib_params.end()) continue;

    auto [slope, intercept] = calib_params[ch_name];
    double threshold = 0.5 * slope + intercept;

    for (size_t i = 0; i < waveforms.size() && saved_count < max_to_save; ++i) {
        const auto& [times_vec, amps_vec] = waveforms[i];
        int size = static_cast<int>(times_vec.size());
        if (size < 2) continue;

        // Analizzo la waveform
        auto res = TimeStuff::analyzeWaveform(times_vec.data(), amps_vec.data(), size, threshold);
        if (res.max_index == -1 || res.crossing_time < 0) continue;

        double max_amp = amps_vec[res.max_index];
        double crossing_time = res.crossing_time;

        // Creo TGraph
        TGraph* gr = new TGraph(size, times_vec.data(), amps_vec.data());
        gr->SetTitle(TString::Format("%s waveform %zu;Time (ns);Amplitude (mV)", ch_name.Data(), i));

        // Creo TPaveText
        TPaveText* info = new TPaveText(0.6, 0.7, 0.95, 0.85, "NDC");
        info->AddText(TString::Format("Max amplitude = %.2f mV", max_amp));
        info->AddText(TString::Format("Crossing threshold = %.2f ns", crossing_time));

        // Creo la canvas
        TCanvas* c = new TCanvas(TString::Format("%s_canvas_%zu", ch_name.Data(), i),
                                 TString::Format("%s waveform %zu", ch_name.Data(), i),
                                 800, 600);
        gr->Draw("AL");       // Linea con assi
        info->Draw();         // TPaveText sopra il grafico
        c->Write();           // Scrivo la canvas nel file ROOT

        // Pulizia memoria
        delete gr;
        delete info;
        delete c;

        ++saved_count;
    }

    if (saved_count >= max_to_save) break;
}*/
