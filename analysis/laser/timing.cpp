#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>
// graphics.h  (o dove usi TH2 / TH2F)
#include <TH2.h>     // definisce TH2, TH2F, TH2D …

#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/laser/laser_stuff/gauss_stuff.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/utils.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis//root/graphics.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/laser/laser_stuff/AmplitudeProcessor.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/laser/laser_stuff/DataExtractor.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/laser/laser_stuff/timing-crossing.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/laser/laser_stuff/timing-analysis.h"



// -----------------------------------------------------------------------------
//  Salva istogrammi e fit in un file ROOT esterno (nome passato dall’utente)
// -----------------------------------------------------------------------------
void timing(const TString& filename,
                  const std::map<TString, std::pair<float, float>>& baseline_time_windows,
                  const TString& out_root_name)          // <‑‑ nuovo argomento
{
    // ---------- apertura file dati -------------------------------------------------
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Errore: impossibile aprire file " << filename << std::endl;
        return;
    }

    // ---------- branch laser (opzionale) -------------------------------------------
    TTree*    tree_laser  = dynamic_cast<TTree*>(file->Get("laser"));
    Double_t* laser_amp   = nullptr;
    Double_t* laser_time  = nullptr;
    if (tree_laser) {
        laser_amp  = Utils::setupBranch_dgz(tree_laser, "amplitude");
        laser_time = Utils::setupBranch_dgz(tree_laser, "time");
    }

    // ---------- estrazione dati canale ---------------------------------------------
    std::map<TString, TTree*>     tree_map;
    std::map<TString, Double_t*>  time_map, amp_map;
    DataExtractor extractor(file);
    if (!extractor.extract(tree_map, time_map, amp_map)) {
        file->Close();
        return;
    }

    // ---------- correzione baseline ------------------------------------------------
    AmplitudeResults results = compute_corrected_waveforms(tree_map,time_map,amp_map,baseline_time_windows);

    auto& wf_map   = results.corrected_waveforms_map;
    auto& med_map  = results.baseline_medians_map;

    // ---------- calcolo massimi veri ----------------------------------------------
    PeakResults pk = compute_peak_info(wf_map, tree_map, time_map, med_map, baseline_time_windows);
    auto& max_map     = pk.peak_amplitudes_map;
    auto& pos_map     = pk.peak_positions_map;
    auto& rms_histos  = pk.rms_hist_map;

    // ---------- fit e istogrammi ---------------------------------------------------
    std::map<TString, std::pair<TH1F*, TF1*>> amp_dist;

    for (const auto& [ch, arr] : max_map) {
        if (!arr || tree_map.find(ch) == tree_map.end()) continue;

        auto [hist, fit] = amplitude_histo(arr, tree_map.at(ch),ch);
        if (hist && fit) amp_dist[ch] = {hist, fit};
    }
    
    // ---------- fit per half p.e. ---------------------------------------------------
    auto [graph, linear_fit] = linear_fit_correction(amp_dist);
    if (!graph || !linear_fit) {
        std::cerr << "Errore nel creare il grafico o il fit." << std::endl;
        return;
    }
    double thr_halfPE = linear_fit->GetParameter(0)*0.5 + linear_fit->GetParameter(1);

    /* ===== 2. crossing per i canali (soglia half‑p.e.) ========================= */
    std::map<TString, std::pair<Double_t*, Double_t*>> t_cross_map =
        TimingProcessor::compute_threshold_crossings(wf_map, time_map, tree_map, thr_halfPE);
    auto all_crossings = TimingProcessor::compute_all_crossings(wf_map, time_map, tree_map, thr_halfPE);

    /* ===== 3. crossing per il laser (soglia half‑max evento per evento) ======= */
    Double_t* t_laser = TimingProcessor::compute_laser_crossings(tree_laser, laser_amp, laser_time);
    if (!t_laser) {
        std::cerr << "Laser crossing non disponibile." << std::endl;
    }

    /* ===== 4. istogrammi Δt (canale – laser) ================================== */
    TimingAnalysis::Result dt_results = TimingAnalysis::Fixed_thr(
            t_cross_map, all_crossings, t_laser, tree_map, tree_laser);

    // ---------- salvataggio su file ROOT -------------------------------------------
    TFile outFile(out_root_name, "RECREATE");
    if (outFile.IsZombie()) {
        std::cerr << "Errore: impossibile creare " << out_root_name << std::endl;
    } else {
        // Salva gli istogrammi e i fit esistenti
        for (const auto& [ch, pairHF] : amp_dist) {
            pairHF.first->Write();   // TH1F
        }
        // Salva anche gli istogrammi RMS
        for (const auto& [ch, rms_hist] : rms_histos) {
            if (rms_hist) rms_hist->Write();
        }
        // Salva il TGraphErrors e il fit lineare
        graph->Write("amplitude_mean_vs_point");
        // Salva gli istogrammi Δt calcolati dalla classe TimingAnalysis
        for (const auto& [ch, h_first] : dt_results.dt_first_map) {
            if (h_first) h_first->Write();
        }
        for (const auto& [ch, h_all] : dt_results.dt_all_map) {
            if (h_all) h_all->Write();
        }
        if (dt_results.h_laser_time) dt_results.h_laser_time->Write();
        // SALVA i grafici 2D crossing time vs amplitude
        for (const auto& [ch, h2] : dt_results.h2_amp_vs_dt_map) {
            if (h2) h2->Write();
        }
        
        outFile.Close();
    }

    // ---------- cleanup memoria dinamica ------------------------------------------
    for (auto& [ch, arr] : max_map) delete[] arr;
    for (auto& [ch, arr] : wf_map)  delete[] arr;
    for (auto& [ch, arr] : med_map) delete[] arr;

    file->Close();
}
