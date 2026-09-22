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



#include <TFile.h>
#include <TH1D.h>
#include <vector>
#include <algorithm>


// =====================================================
// Funzione principale
// =====================================================
void timing(const TString& filename, const TString& out_root_name)
{
    // =====================================================
    //                PARAMETRI MODIFICABILI
    // =====================================================
    const double t_min = 38;
    const double t_max = 43;
    
    const int    n_bins_max_amp      = 2500;   // numero di bin
    const double lower_range_max_amp = -5.0;    // valore minimo asse x
    const double upper_range_max_amp = 495;  // valore massimo asse x
    const int    n_bins_slope_amp      = 1000;   // numero di bin
    const double lower_range_slope_amp = 0;    // valore minimo asse x
    const double upper_range_slope_amp = 50;  // valore massimo asse x
    const int    n_bins_delta_t      = 20000;   // numero di bin
    const double lower_range_delta_t = -500;    // valore minimo asse x
    const double upper_range_delta_t = 500;  // valore massimo asse x
    const int n_peaks_requested = 2; // ad esempio
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
            auto aligned = waveAligner.getWaveforms(ch_tree, ch_time, ch_amp, true,true);
            auto corrected_aligned = waveAligner.correctWaveforms(aligned, 30.0);
            
            all_aligned_waveforms[ch_name] = std::move(corrected_aligned);
            
            TimeStuff::save_waveforms(all_aligned_waveforms);

            
        }
        else {
            std::cerr << "Skipping channel " << ch_name << " due to missing data pointers\n";
        }
    }
    
    
    ///baseline
    // Istogramma della baseline
    TH1D* hBaseline = new TH1D(
        "hBaseline",
        "Baseline distribution;Amplitude;Counts",
        100,   // bin
        -10,  // min (adatta al tuo range)
        10    // max
    );
    for (const auto& [ch_name, waveforms] : all_aligned_waveforms) {

        for (const auto& wf : waveforms) {

            const auto& time = wf.first;
            const auto& amp  = wf.second;

            for (size_t i = 0; i < time.size(); ++i) {
                if (time[i] < 30.0) {
                    hBaseline->Fill(amp[i]);
                }
            }
        }
    }

    hBaseline->Write();
    
    
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
            TF1* linearFit = new TF1((ch_name + "_linearFit").Data(), "pol1", 0, n_peaks+1);
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
            //double first_mean = fit->GetParameter(4);       // parametro 1 = media della prima gaussiana
            //double first_mean_err = fit->GetParError(4);   // errore della media della prima gaussiana
            double first_mean = 1.5;       // parametro 1 = media della prima gaussiana
            double first_mean_err = 0;   // errore della media della prima gaussiana
            
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
    TDirectory* dir_wf_baseline = outFile->mkdir("baseline");
    TDirectory* cfd = outFile->mkdir("constant_fraction_discrimination");
    TDirectory* dir_crossing = outFile->mkdir("clean_wf");
    TDirectory* dir_mean_wf = outFile->mkdir("avarage_clean_wf");
    
    
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
        std::cout << "thr " << threshold_fixed << "\n";
        // istogrammi principali
        auto h_fixed_th_diff = std::make_unique<TH2F>(ch_name + "_fixed_th_diff_laser",
                                                      ch_name + "_fixed_th_diff_laser;amplitude (mV);#Delta t (ns)",
                                                      n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,
                                                      n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        
        // --- dichiarazioni istogrammi (fuori dal loop, insieme agli altri) ---
        auto h_cfd_dt_vs_amp = std::make_unique<TH2F>(ch_name + "_cfd_dt_vs_amp", ch_name +
                                                      "_cfd;amplitude (mV);#Deltat_{CFD} (ns)",n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp,upper_range_max_amp);
        
        auto h_linear_dt_vs_amp = std::make_unique<TH2F>(ch_name + "_linear_dt_vs_amp", ch_name +
                                                         "_linear;amplitude (mV);#Deltat_{lin.} (ns)",
                                                         n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        
        auto h_linear_dt_vs_slope = std::make_unique<TH2F>(ch_name + "_linear_dt_vs_slope", ch_name +
                                                           "_slope;linear slope;#Deltat_{lin.} (ns)",
                                                           n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        
        
        
        auto h_wf_all = std::make_unique<TH2F>(ch_name + "_all_waveforms",
                                               ch_name + "_all_waveforms;time (ns);amplitude (mV)",
                                               3000, -50, 250, n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        
        auto h_wf_selected = std::make_unique<TH2F>(ch_name + "_waveforms_selected",
                                                    ch_name + "_waveforms_selected;time (ns);amplitude (mV)",
                                                    n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp, upper_range_max_amp);
        
        
        // soglie classificazione
        std::vector<double> thresholds = {0.5*slope+intercept, 1.5*slope+intercept, 2.5*slope+intercept, 3.5*slope+intercept, 4.5*slope+intercept, 5.5*slope+intercept};
        size_t n_classes = thresholds.size() - 1;
        
        // istogrammi 2D per classificazioni
        std::vector<std::unique_ptr<TH2F>> h_class_2D;
        for (size_t cls=0; cls<n_classes; ++cls) {
            h_class_2D.push_back(std::make_unique<TH2F>(
                                                        TString::Format("%s_class_%zu", ch_name.Data(), cls+1),
                                                        TString::Format("%s_class_%zu;time (ns);amplitude (mV)", ch_name.Data(), cls+1),
                                                        n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp, upper_range_max_amp
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
        
        // Istogrammi 2D per classificazioni con crossing time selezionato
        std::vector<std::unique_ptr<TH2F>> h_crossing_2D;
        for (size_t cls = 0; cls < n_classes; ++cls) {
            h_crossing_2D.push_back(std::make_unique<TH2F>(
                                                           TString::Format("%s_crossing_class_%zu", ch_name.Data(), cls+1),
                                                           TString::Format("%s_crossing_class_%zu;time (ns);amplitude (mV)", ch_name.Data(), cls+1),
                                                           n_bins_delta_t, lower_range_delta_t, upper_range_delta_t,n_bins_max_amp, lower_range_max_amp, upper_range_max_amp
                                                           ));
        }
        
        static TH1D* h_max_amp_sel = nullptr;
        static TFile* f_max_amp = nullptr;

        if (!h_max_amp_sel) {
            f_max_amp = new TFile("max_amp_selected.root", "RECREATE");
            h_max_amp_sel = new TH1D(
                "h_max_amp_sel",
                "Selected max amplitudes;max_amp;Entries",
                500, -10, 90  // <-- modifica range se serve
            );
        }
        
        std::vector<TDirectory*> class_dirs(n_classes);
        std::vector<size_t> graph_counter(n_classes, 0);
        for (size_t cls = 0; cls < n_classes; ++cls) {

            // Se esiste già la directory la recupera, altrimenti la crea
            class_dirs[cls] = (TDirectory*)outFile->Get(Form("class_%zu", cls));
            
            if (!class_dirs[cls])
                class_dirs[cls] = outFile->mkdir(Form("%s_class_%zu", ch_name.Data(), cls));
        }
        
        size_t n_events = std::min(waveforms.size(), laser_waveforms.size());
        // Vettori temporanei per spline cubica
        std::vector<double> slope_linear_all;
        std::vector<double> delta_t_linear_all;
        std::vector<double> max_amp_all;
        
        // --- Fuori dal loop sui waveform, inizializza vettore per classi ---
        // Imposta binning
        int n_bins_wf = 1024;
        double t_min_wf = -50.0;
        double t_max_wf = 200.0;
        double bin_width = (t_max_wf - t_min_wf) / n_bins_wf;
        std::vector<std::vector<std::vector<double>>> wf_selected_by_class(n_classes);
        for (size_t cls = 0; cls < n_classes; ++cls)
            wf_selected_by_class[cls].resize(n_bins_wf);
        
        
        
        // --- Loop unificato sui waveform ---
        for (size_t i = 0; i < waveforms.size(); ++i) {
            const auto& [times, amps_vec] = waveforms[i];
            
            size_t size = times.size();
            
            // Riempio gli histogrammi WF → tempi già shiftati
            for (size_t j = 0; j < size; ++j)
                h_wf_all->Fill(times[j], amps_vec[j]);
            
            
            
            auto res_no_constraint = TimeStuff::analyzeWaveform(
                                                                times.data(),
                                                                amps_vec.data(),
                                                                static_cast<int>(size),
                                                                threshold_fixed
                                                                );
            
            
            if (res_no_constraint.max_index < 0 || !std::isfinite(res_no_constraint.crossing_time)) {
                continue; // oppure continue
            }

            double max_amp = amps_vec[res_no_constraint.max_index];
            double delta_t = res_no_constraint.crossing_time;

            h_fixed_th_diff->Fill(max_amp, delta_t);
            
            
            auto res_cfd = TimeStuff::computeCFDCrossing(
                                                         times.data(),
                                                         amps_vec.data(),
                                                         static_cast<int>(size),
                                                         res_no_constraint.max_index,
                                                         0.5
                                                         );
            double delta_t_cfd = res_cfd.crossing_time;
            h_cfd_dt_vs_amp->Fill(amps_vec[res_cfd.max_index], delta_t_cfd);
            
            
            auto [t_lin_abs, slope_linear] = TimeStuff::computeLinearFitWithSlope(
                                                                                  times.data(),
                                                                                  amps_vec.data(),
                                                                                  static_cast<int>(size),
                                                                                  res_no_constraint.max_index
                                                                                  );
            double delta_t_linear = t_lin_abs;
            h_linear_dt_vs_amp->Fill(max_amp, delta_t_linear);
            h_linear_dt_vs_slope->Fill(slope_linear, delta_t_linear);
            
            slope_linear_all.push_back(slope_linear);
            delta_t_linear_all.push_back(delta_t_linear);
            max_amp_all.push_back(max_amp);
            
            
            
            
            
            for (size_t j = 0; j < size; ++j)
                h_wf_selected->Fill(times[j], amps_vec[j]);
            
            for (size_t cls = 0; cls < n_classes; ++cls) {
                wf_selected_by_class[cls].resize(n_bins_wf);
            }
            
            // --- FILTRO crossing 36-39 ns ---
            if (delta_t < 37 || delta_t > 40.0) continue;
            double true_max_amp = *std::max_element(amps_vec.begin(), amps_vec.end());
            if (max_amp < true_max_amp || max_amp > true_max_amp) continue;
            h_max_amp_sel->Fill(max_amp);
            
            // --- Selezione classe ---
            size_t cls_selected = n_classes; // default = non classificato
            for (size_t cls = 0; cls < n_classes; ++cls) {
                if (max_amp >= thresholds[cls] && max_amp < thresholds[cls + 1]) {
                    cls_selected = cls;
                    break;
                }
            }
            if (cls_selected == n_classes) continue;
            
            // --- Popola wf_selected_by_class ---
            for (size_t j = 0; j < times.size(); ++j) {
                int bin_idx = std::clamp(int((times[j] - t_min_wf) / bin_width), 0, n_bins_wf - 1);
                wf_selected_by_class[cls_selected][bin_idx].push_back(amps_vec[j]);
            }
            
            // --- Popola anche l'istogramma 2D dei crossing ---
            for (size_t j = 0; j < times.size(); ++j) {
                h_crossing_2D[cls_selected]->Fill(times[j], amps_vec[j]);
            }
            
            for (size_t cls = 0; cls < n_classes; ++cls) {

                if (max_amp >= thresholds[cls] && max_amp < thresholds[cls + 1]) {

                    // ====== Istogrammi (come già fai) ======
                    for (size_t j = 0; j < size; ++j)
                        h_class_2D[cls]->Fill(times[j], amps_vec[j]);

                    for (size_t j = 0; j < size; ++j)
                        if (times[j] <= 30.0)
                            h_class_amp_baseline[cls]->Fill(amps_vec[j]);

                    // ====== Salvataggio waveform come TGraph ======
                    class_dirs[cls]->cd();

                    TGraph gr(size, times.data(), amps_vec.data());
                    gr.SetName(Form("wf_%zu", graph_counter[cls]++));
                    gr.SetTitle(Form("Waveform class %zu", cls));

                    gr.Write();

                    outFile->cd();
                    break;
                }
            }
            
        }
        
        
        
        
        // --- Scrivi istogrammi 2D ---
        dir_crossing->cd();
        for (size_t cls = 0; cls < n_classes; ++cls) {
            h_crossing_2D[cls]->Write();
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
        
        /// --- Scrivi tutti gli istogrammi su file ROOT ---
        dir_fixed_th->cd();
        h_fixed_th_diff->Write();
        
        dir_wf_all->cd();
        h_wf_all->Write();
        h_wf_selected->Write();
        
        cfd->cd();
        h_cfd_dt_vs_amp->Write();
        h_linear_dt_vs_amp->Write();        // pre-correzione
        h_linear_dt_vs_slope->Write();      // pre-correzione
        
        if (f_max_amp && h_max_amp_sel) {
            f_max_amp->cd();
            h_max_amp_sel->Write();
            f_max_amp->Close();
        }
        
        // --- Costruzione dei TGraph (baseline-corrected) ---
        // --- Costruzione TGraphErrors (baseline-corrected) con 1024 bin ---
        dir_mean_wf->cd();
        
        double t_baseline_min = -40.0;
        double t_baseline_max = 0.0;
        
        for (size_t cls = 0; cls < n_classes; ++cls) {
            std::vector<double> times(n_bins_wf);
            std::vector<double> amps(n_bins_wf, 0.0);
            std::vector<double> errs(n_bins_wf, 0.0);
            
            // tempi dei bin (centri)
            for (int b = 0; b < n_bins_wf; ++b)
                times[b] = t_min_wf + (b + 0.5) * bin_width;
            
            // --- baseline dai bin pre-trigger ---
            std::vector<double> baseline_samples;
            int bmin = std::max(0, int(std::floor((t_baseline_min - t_min_wf) / bin_width)));
            int bmax = std::min(n_bins_wf - 1, int(std::floor((t_baseline_max - t_min_wf) / bin_width)));
            for (int b = bmin; b <= bmax; ++b)
                baseline_samples.insert(baseline_samples.end(),
                                        wf_selected_by_class[cls][b].begin(),
                                        wf_selected_by_class[cls][b].end());
            
            double baseline_median = 0.0;
            if (!baseline_samples.empty()) {
                std::sort(baseline_samples.begin(), baseline_samples.end());
                size_t m = baseline_samples.size() / 2;
                baseline_median = (baseline_samples.size() % 2 == 0) ?
                0.5 * (baseline_samples[m - 1] + baseline_samples[m]) :
                baseline_samples[m];
            }
            
            // --- calcola mediana e errore per ciascun bin ---
            for (int b = 0; b < n_bins_wf; ++b) {
                auto &v = wf_selected_by_class[cls][b];
                if (v.empty()) continue;
                
                // mediana del bin
                std::sort(v.begin(), v.end());
                size_t mid = v.size() / 2;
                double median_val = (v.size() % 2 == 0) ? 0.5 * (v[mid - 1] + v[mid]) : v[mid];
                
                // sottrai baseline
                amps[b] = median_val - baseline_median;
                
                // errore: deviazione standard / sqrt(N)
                double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
                double var = 0.0;
                for (double vv : v) var += (vv - mean) * (vv - mean);
                if (v.size() > 1) var /= (v.size() - 1);
                errs[b] = (v.size() > 0) ? std::sqrt(var) / std::sqrt(v.size()) : 0.0;
            }
            
            // --- crea TGraphErrors ---
            TGraphErrors *tg = new TGraphErrors(n_bins_wf,
                                                times.data(),
                                                amps.data(),
                                                nullptr,   // ex = 0
                                                errs.data());
            tg->SetName(TString::Format("%s_mean_wf_class_%zu", ch_name.Data(), cls + 1));
            tg->SetTitle(TString::Format("%s_mean_wf_class_%zu;time (ns);amplitude (mV)", ch_name.Data(), cls + 1));
            tg->SetLineColor(kRed);
            tg->SetLineWidth(2);
            tg->Write();
            delete tg;
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

    
    
    outFile->Close();
}


