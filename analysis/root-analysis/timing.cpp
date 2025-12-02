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
    const double t_min = 36.0;
    const double t_max = 41.0;
    
    const int    n_bins_max_amp      = 1000.0;   // numero di bin
    const double lower_range_max_amp = -200.0;    // valore minimo asse x
    const double upper_range_max_amp = 200.0;  // valore massimo asse x
    const int    n_bins_slope_amp      = 1000;   // numero di bin
    const double lower_range_slope_amp = 0;    // valore minimo asse x
    const double upper_range_slope_amp = 50;  // valore massimo asse x
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
    auto outFile = std::unique_ptr<TFile>(TFile::Open(out_root_name, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: cannot create output file " << out_root_name << std::endl;
        return;
    }
    std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_aligned_waveforms;
    std::map<TString, std::vector<std::pair<std::vector<double>, std::vector<double>>>> all_baseline_corrected_waveforms;
    for (const auto& [ch_name, ch_tree] : tree_map) {
        auto ch_time = time_map[ch_name];
        auto ch_amp  = amp_map_raw[ch_name];
        
        // Dentro al tuo loop for (const auto& [ch_name, ch_tree] : tree_map) {
        if (ch_tree && ch_time && ch_amp) {
            auto aligned = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, true);
            auto corrected_aligned = waveAligner.correctWaveforms(aligned, 30.0);
            
            all_aligned_waveforms[ch_name] = std::move(corrected_aligned);
            std::cout << "Aligned and baseline-corrected channel " << ch_name
            << " with " << all_aligned_waveforms[ch_name].size() << " waveforms\n";
            
            // === Crea due folder per salvare i TGraph ===
            TString dirBeforeName = Form("waveforms_before_%s", ch_name.Data());
            TString dirAfterName  = Form("waveforms_after_%s",  ch_name.Data());
            
            TDirectory* dir_before = outFile->mkdir(dirBeforeName);
            TDirectory* dir_after  = outFile->mkdir(dirAfterName);
            
            // Limita a 500 wf o al numero disponibile
            size_t nToSave = std::min<size_t>(10000, aligned.size());
            
            for (size_t i = 0; i < nToSave; i++) {
                const auto& wf_before = aligned[i];                    // prima della correzione
                const auto& wf_after  = all_aligned_waveforms[ch_name][i]; // dopo la correzione
                
                // --- Salva in "before" ---
                dir_before->cd();
                TGraph* gr_before = new TGraph(wf_before.first.size(),
                                               wf_before.first.data(),
                                               wf_before.second.data());
                gr_before->SetName(Form("wf_before_%s_%zu", ch_name.Data(), i));
                gr_before->SetTitle(Form("Waveform BEFORE correction ch %s #%zu", ch_name.Data(), i));
                gr_before->Write();
                
                // --- Salva in "after" ---
                dir_after->cd();
                TGraph* gr_after = new TGraph(wf_after.first.size(),
                                              wf_after.first.data(),
                                              wf_after.second.data());
                gr_after->SetName(Form("wf_after_%s_%zu", ch_name.Data(), i));
                gr_after->SetTitle(Form("Waveform AFTER correction ch %s #%zu", ch_name.Data(), i));
                gr_after->Write();
                // --- Salva laser corrispondente con scaling ---
                if (i < laser_waveforms.size()) {
                    const auto& laser_wf = laser_waveforms[i];
                    std::vector<double> laser_scaled(laser_wf.second.size());
                    
                    for (size_t j = 0; j < laser_wf.second.size(); ++j) {
                        laser_scaled[j] = laser_wf.second[j] / 10.0;  // riscalo
                    }
                    
                    TGraph* gr_laser = new TGraph(laser_wf.first.size(),
                                                  laser_wf.first.data(),
                                                  laser_scaled.data()); // uso il vettore scalato
                    gr_laser->SetName(Form("laser_wf_%s_%zu", ch_name.Data(), i));
                    gr_laser->SetTitle(Form("Laser waveform corresponding to ch %s #%zu", ch_name.Data(), i));
                    gr_laser->Write();
                }
                
            }
        } else {
            std::cerr << "Skipping channel " << ch_name << " due to missing data pointers\n";
        }
    }
    
    
    // File di output
    
    TDirectory* dir_amp = outFile->mkdir("amp_max-distribution");
    TDirectory* dir_cal = outFile->mkdir("calibration");
    
    //Fare il lavoro del fit multi gauss a mean fissa e sigma fissa
    
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
            // --- Stampa calibrazione ---
            double slope     = linearFit->GetParameter(1);
            double intercept = linearFit->GetParameter(0);
            
            std::cout << "[Calibrazione] " << ch_name
            << " slope = " << slope
            << " , intercept = " << intercept << std::endl;
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
    TDirectory* dir_wf_all = outFile->mkdir("all_waveforms");
    TDirectory* dir_wf = outFile->mkdir("classification_first_crossing");
    TDirectory* dir_wf_profile = outFile->mkdir("profile_classification_first_crossing");
    TDirectory* dir_wf_baseline = outFile->mkdir("baseline");
    TDirectory* cfd = outFile->mkdir("constant_fraction_discrimination");
    
    
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
        double threshold_fixed = 0.5 * slope + intercept;
        
        // istogrammi principali
        auto h_fixed_th_diff = std::make_unique<TH2F>(ch_name + "_fixed_th_diff_laser",
                                                      ch_name + "_fixed_th_diff_laser;amplitude (mV);#Delta t (ns)",
                                                      n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
                                                      n_bins_delta_t, lower_range_delta_t, upper_range_delta_t);
        
        // --- dichiarazioni istogrammi (fuori dal loop, insieme agli altri) ---
        auto h_cfd_dt_vs_amp = std::make_unique<TH2F>(ch_name + "_cfd_dt_vs_amp", ch_name +
                                                      "_cfd;amplitude (mV);#Deltat_{CFD} (ns)",n_bins_max_amp, lower_range_max_amp,upper_range_max_amp,n_bins_delta_t, lower_range_delta_t, upper_range_delta_t);
        
        auto h_linear_dt_vs_amp = std::make_unique<TH2F>(ch_name + "_linear_dt_vs_amp", ch_name +
                                                         "_linear;amplitude (mV);#Deltat_{lin.} (ns)",
                                                         n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
                                                         n_bins_delta_t, lower_range_delta_t, upper_range_delta_t);
        
        auto h_linear_dt_vs_slope = std::make_unique<TH2F>(ch_name + "_linear_dt_vs_slope", ch_name +
                                                           "_slope;linear slope;#Deltat_{lin.} (ns)",
                                                           n_bins_max_amp, lower_range_max_amp, upper_range_max_amp,
                                                           n_bins_delta_t, lower_range_delta_t, upper_range_delta_t);
      
        
        
        auto h_wf_all = std::make_unique<TH2F>(ch_name + "_all_waveforms",
                                               ch_name + "_all_waveforms;time (ns);amplitude (mV)",
                                               1250, -50, 200, 3000.0, -50, 250);
        
        auto h_wf_selected = std::make_unique<TH2F>(ch_name + "_waveforms_selected",
                                                    ch_name + "_waveforms_selected;time (ns);amplitude (mV)",
                                                    1250, -50, 200, 3000.0, -50, 250);
        
        
        // soglie classificazione
        std::vector<double> thresholds = {0.5*slope+intercept, 1.5*slope+intercept, 2.5*slope+intercept, 3.5*slope+intercept, 4.5*slope+intercept, 5.5*slope+intercept};
        size_t n_classes = thresholds.size() - 1;
        
        // istogrammi 2D per classificazioni
        std::vector<std::unique_ptr<TH2F>> h_class_2D;
        for (size_t cls=0; cls<n_classes; ++cls) {
            h_class_2D.push_back(std::make_unique<TH2F>(
                                                        TString::Format("%s_class_%zu", ch_name.Data(), cls+1),
                                                        TString::Format("%s_class_%zu;time (ns);amplitude (mV)", ch_name.Data(), cls+1),
                                                        1250, -50, 200, 3000.0, -50, 250
                                                        ));
        }
        
        // Istogrammi 1D delle ampiezze fino a 30 ns per baseline
        std::vector<std::unique_ptr<TH1F>> h_class_amp_baseline;
        for (size_t cls = 0; cls < n_classes; ++cls) {
            h_class_amp_baseline.push_back(std::make_unique<TH1F>(
                                                                  TString::Format("%s_class_%zu_amp_baseline", ch_name.Data(), cls+1),
                                                                  TString::Format("%s_class_%zu_amp_baseline;amplitude (mV);entries", ch_name.Data(), cls+1),
                                                                  500, -50, 50  // range ampiezza, regolabile se necessario
                                                                  ));
        }
        
        size_t n_events = std::min(waveforms.size(), laser_waveforms.size());
        // Vettori temporanei per spline cubica
        std::vector<double> slope_linear_all;
        std::vector<double> delta_t_linear_all;
        std::vector<double> max_amp_all;
        
        // --- Loop sui waveform ---
        for (size_t i = 0; i < waveforms.size(); ++i) {
            const auto& [times_vec, amps_vec] = waveforms[i];
            const auto& [laser_times, laser_amps] = laser_waveforms[i];
            if (times_vec.size() < 2 || laser_times.size() < 2) continue;
            size_t size = times_vec.size();
            
            for (size_t j = 0; j < size; ++j) h_wf_all->Fill(times_vec[j], amps_vec[j]);
            
            double crossing_time_laser = TimeStuff::computeCrossingTime(laser_times.data(), laser_amps.data(), laser_times.size(), laser_threshold);
            if (crossing_time_laser < 0) continue;
            
            auto res_no_constraint = TimeStuff::analyzeWaveform(times_vec.data(), amps_vec.data(), static_cast<int>(size), threshold_fixed);
            if (res_no_constraint.max_index < 0) continue;
            
            double delta_t = res_no_constraint.crossing_time;
            double max_amp = amps_vec[res_no_constraint.max_index];
            h_fixed_th_diff->Fill(max_amp, delta_t);
            
            auto res_cfd = TimeStuff::computeCFDCrossing(times_vec.data(),
                                                         amps_vec.data(),
                                                         static_cast<int>(size),
                                                         res_no_constraint.max_index,
                                                         0.5);
            double delta_t_cfd = res_cfd.crossing_time;
            h_cfd_dt_vs_amp->Fill(amps_vec[res_cfd.max_index], delta_t_cfd);
            
            auto [delta_t_linear, slope_linear] = TimeStuff::computeLinearFitWithSlope(
                                                                                       times_vec.data(),
                                                                                       amps_vec.data(),
                                                                                       static_cast<int>(size),
                                                                                       res_no_constraint.max_index
                                                                                       );
            // --- Classificazione ---
            int cls_selected = -1;
               
            for (size_t cls = 0; cls < n_classes; ++cls) {
                if (max_amp >= thresholds[cls] && max_amp < thresholds[cls + 1]) {
                    // Riempi sempre l'istogramma base della classe
                    for (size_t j = 0; j < size; ++j)
                        h_class_2D[cls]->Fill(times_vec[j], amps_vec[j]);
                    
                    // Riempie l'istogramma 1D per baseline solo fino a 30 ns
                    for (size_t j = 0; j < size; ++j) {
                        if (times_vec[j] <= 30.0)
                            h_class_amp_baseline[cls]->Fill(amps_vec[j]);
                    }
                    break; // trovato e gestito, si può uscire
                }
            }

                    
            h_linear_dt_vs_amp->Fill(max_amp, delta_t_linear);
            h_linear_dt_vs_slope->Fill(slope_linear, delta_t_linear);
            
            slope_linear_all.push_back(slope_linear);
            delta_t_linear_all.push_back(delta_t_linear);
            max_amp_all.push_back(max_amp);
            
            for (size_t j = 0; j < size; ++j) h_wf_selected->Fill(times_vec[j], amps_vec[j]);
        }
        
        // --- Scrivi tutti gli istogrammi su file ROOT ---
        dir_fixed_th->cd(); h_fixed_th_diff->Write();
        dir_wf_all->cd(); h_wf_all->Write(); h_wf_selected->Write();
        cfd->cd();
        h_cfd_dt_vs_amp->Write();
        h_linear_dt_vs_amp->Write();                 // pre-correzione
        h_linear_dt_vs_slope->Write();               // pre-correzione
        
        dir_wf_profile->cd();
        for (size_t cls = 0; cls < n_classes; ++cls) {
            // Istogramma base
            auto prof_class = h_class_2D[cls]->ProfileX(TString::Format("%s_class_%zu_profile", ch_name.Data(), cls + 1),1, -1, "e");
            prof_class->SetLineColor(kBlack);
            prof_class->SetLineWidth(2);
        }
     
        
        dir_wf->cd();
        for (size_t cls = 0; cls < n_classes; ++cls) {
            h_class_2D[cls]->Write();
        }
        dir_wf_baseline->cd();
        for (size_t cls = 0; cls < n_classes; ++cls) {
            // Istogramma 2D (tempo vs ampiezza)
            h_class_amp_baseline[cls]->Write();
        }

    }


        
    outFile->cd();
    // Creo TTree "info" per memorizzare numer o di waveforms analizzate
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


/*std::cout << "Aligned and baseline-corrected channel " << ch_name
 << " with " << all_aligned_waveforms[ch_name].size() << " waveforms\n";
 
 // === Crea due folder per salvare i TGraph ===
 TString dirBeforeName = Form("waveforms_before_%s", ch_name.Data());
 TString dirAfterName  = Form("waveforms_after_%s",  ch_name.Data());
 
 TDirectory* dir_before = outFile->mkdir(dirBeforeName);
 TDirectory* dir_after  = outFile->mkdir(dirAfterName);
 
 // Limita a 500 wf o al numero disponibile
 size_t nToSave = std::min<size_t>(500, aligned.size());
 
 for (size_t i = 0; i < nToSave; i++) {
     const auto& wf_before = aligned[i];                    // prima della correzione
     const auto& wf_after  = all_aligned_waveforms[ch_name][i]; // dopo la correzione
     
     // --- Salva in "before" ---
     dir_before->cd();
     TGraph* gr_before = new TGraph(wf_before.first.size(),
                                    wf_before.first.data(),
                                    wf_before.second.data());
     gr_before->SetName(Form("wf_before_%s_%zu", ch_name.Data(), i));
     gr_before->SetTitle(Form("Waveform BEFORE correction ch %s #%zu", ch_name.Data(), i));
     gr_before->Write();
     
     // --- Salva in "after" ---
     dir_after->cd();
     TGraph* gr_after = new TGraph(wf_after.first.size(),
                                   wf_after.first.data(),
                                   wf_after.second.data());
     gr_after->SetName(Form("wf_after_%s_%zu", ch_name.Data(), i));
     gr_after->SetTitle(Form("Waveform AFTER correction ch %s #%zu", ch_name.Data(), i));
     gr_after->Write();
 }*/



/*
 size_t n_canvas_saved = 0;
 if (res_no_constraint.crossing_time < 30.0 && cls_selected==0 && n_canvas_saved<10) {
     TString graph_name  = TString::Format("%s_wf_%zu", ch_name.Data(), i);
     TString graph_title = TString::Format("%s waveform %zu;time (ns);amplitude (mV)", ch_name.Data(), i);

     auto gr = std::make_unique<TGraph>(static_cast<int>(size), times_vec.data(), amps_vec.data());
     gr->SetName(graph_name);
     gr->SetTitle(graph_title);

     // crossing point
     double crossing_time = res_no_constraint.crossing_time;
     double crossing_amp = 0.0;
     for (size_t j = 1; j < size; ++j) {
         if (times_vec[j-1] <= crossing_time && crossing_time <= times_vec[j]) {
             double t0 = times_vec[j-1], t1 = times_vec[j];
             double a0 = amps_vec[j-1], a1 = amps_vec[j];
             crossing_amp = a0 + (crossing_time - t0)*(a1-a0)/(t1-t0);
             break;
         }
     }

     auto marker_crossing = std::make_unique<TMarker>(crossing_time, crossing_amp, 20);
     marker_crossing->SetMarkerColor(kRed);
     marker_crossing->SetMarkerSize(1.2);

     // massimo
     double max_time = times_vec[res_no_constraint.max_index];
     double max_val = amps_vec[res_no_constraint.max_index];
     auto marker_max = std::make_unique<TMarker>(max_time, max_val, 21);
     marker_max->SetMarkerColor(kBlue);
     marker_max->SetMarkerSize(1.2);

     // canvas
     auto c = std::make_unique<TCanvas>(TString::Format("c_%s_%zu", ch_name.Data(), i),
                                        TString::Format("Waveform %zu %s", i, ch_name.Data()),
                                        800, 600);
     c->cd();
     gr->Draw("AL");
     marker_crossing->Draw("P");
     marker_max->Draw("P");

     // linee soglia e classificazione
     TLine line_thresh( times_vec.front(), threshold_fixed, times_vec.back(), threshold_fixed );
     line_thresh.SetLineColor(kGreen); line_thresh.SetLineStyle(2); line_thresh.Draw("same");

     TLine line_cls_low( times_vec.front(), thresholds[0], times_vec.back(), thresholds[0] );
     TLine line_cls_high(times_vec.front(), thresholds[1], times_vec.back(), thresholds[1] );
     line_cls_low.SetLineColor(kMagenta); line_cls_low.SetLineStyle(2); line_cls_low.Draw("same");
     line_cls_high.SetLineColor(kMagenta); line_cls_high.SetLineStyle(2); line_cls_high.Draw("same");

     dir_wf->cd();
     c->Write();
     gr->Write();
     marker_crossing->Write(TString::Format("%s_crossing_marker", graph_name.Data()));
     marker_max->Write(TString::Format("%s_max_marker", graph_name.Data()));

     ++n_canvas_saved;
 }
}
 */
