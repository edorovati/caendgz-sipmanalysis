#ifndef TIMING_ANALYSIS_H
#define TIMING_ANALYSIS_H

#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <TH2.h>  
class TimingAnalysis {
public:
    
    struct Result {
        std::map<TString, TH1F*> dt_first_map;
        std::map<TString, TH1F*> dt_all_map;
        TH1F* h_laser_time = nullptr;

        std::map<TString, TH2F*> h2_amp_vs_dt_map;  // <-- istogrammi 2D

        ~Result() {
            for (auto& [_, h] : dt_first_map) if (h) delete h;
            for (auto& [_, h] : dt_all_map) if (h) delete h;
            if (h_laser_time) delete h_laser_time;
            for (auto& [_, h2] : h2_amp_vs_dt_map) if (h2) delete h2;
        }
    };

    
    /**
     * Calcola gli istogrammi Δt tra crossing di canale e laser.
     * @param t_cross_map primo crossing per canale (t_cross_map[ch] è Double_t* per evento)
     * @param all_crossings tutti i crossing per canale (all_crossings[ch][evento] vettore tempi)
     * @param t_laser puntatore array tempi crossing laser per evento
     * @param tree_map mappa canale -> TTree segnale
     * @param tree_laser TTree laser (per numero eventi)
     * @return Result con gli istogrammi dt_first_map, dt_all_map e istogramma distribuzione laser
     */
    static Result Fixed_thr(const std::map<TString, std::pair<Double_t*, Double_t*>>& t_cross_map,
                            const std::map<TString, std::vector<std::vector<double>>>& all_crossings,
                            const Double_t* t_laser,
                            const std::map<TString, TTree*>& tree_map,
                            TTree* tree_laser)
    {
        Result res;

        if (!t_laser || !tree_laser) {
            std::cerr << "TimingAnalysis: t_laser o tree_laser non validi." << std::endl;
            return res;
        }

        Long64_t Nlaser = tree_laser->GetEntries();

        // Istogramma distribuzione tempi laser
        res.h_laser_time = new TH1F("laser_time", "Distribuzione tempi crossing laser;tempo (ns);entries", 4000, 0, 200);
        res.h_laser_time->SetDirectory(nullptr);

        for (Long64_t ie = 0; ie < Nlaser; ++ie) {
            double tl = t_laser[ie];
            if (!std::isnan(tl)) res.h_laser_time->Fill(tl);
        }

        for (int bin = 1; bin <= res.h_laser_time->GetNbinsX(); ++bin) {
            double content = res.h_laser_time->GetBinContent(bin);
            res.h_laser_time->SetBinError(bin, std::sqrt(content));
        }

        for (const auto& [ch, t_amp_pair] : t_cross_map) {
            if (tree_map.find(ch) == tree_map.end()) {
                std::cerr << "TimingAnalysis: canale " << ch << " non trovato in tree_map." << std::endl;
                continue;
            }
            Long64_t Nsig = tree_map.at(ch)->GetEntries();
            if (Nsig != Nlaser) {
                std::cerr << "Warning: #eventi laser != #eventi segnale per canale " << ch << std::endl;
                continue;
            }

            // Istogrammi Δt
            TH1F* h_first = new TH1F(Form("dt_first_crossing_%s", ch.Data()),
                                     Form("#Delta t (1st ‑ laser) %s;#Delta t (ns);entries", ch.Data()),
                                     4000, 0, 200);
            h_first->SetDirectory(nullptr);

            TH1F* h_all = new TH1F(Form("dt_all_crossing_%s", ch.Data()),
                                   Form("#Delta t (all ‑ laser) %s;#Delta t (ns);entries", ch.Data()),
                                   4000, 0, 200);
            h_all->SetDirectory(nullptr);
            const double AMP_MIN = 0.0;
            const double AMP_MAX = 80.0;
            const double DT_MIN  = 30.0;
            const double DT_MAX  = 50.0;
            // Istogramma 2D Ampiezza vs Δt (esempio range ampiezza 0-20, Δt 37-41)
            TH2F* h2_amp_vs_dt = new TH2F(
                Form("h2_amp_vs_dt_fixed_thr_%s", ch.Data()),
                Form("amp_vs_dt - %s;Amplitude (mV);#Delta t (ns)", ch.Data()),
                400, AMP_MIN, AMP_MAX,   // 100 bins amplitude da 0 a 20 mV
                400, DT_MIN, DT_MAX    // 40 bins Δt da 37 a 41 ns
            );
            h2_amp_vs_dt->SetDirectory(nullptr);
            
            // limiti del TH2F (coerenti con la sua definizione)
            
            for (Long64_t ie = 0; ie < Nsig; ++ie) {
                double tl = t_laser[ie];
                if (std::isnan(tl)) continue;

                double ts_first = t_amp_pair.first[ie];
                double amp_first = t_amp_pair.second[ie];

                if (!std::isnan(ts_first) && !std::isnan(amp_first) && amp_first >= AMP_MIN && amp_first <= AMP_MAX && ts_first - tl >= DT_MIN  && ts_first - tl <= DT_MAX)
                {
                    h2_amp_vs_dt->Fill(amp_first, ts_first - tl);
                }

                if (!std::isnan(ts_first)) h_first->Fill(ts_first - tl);

                for (double ts : all_crossings.at(ch)[ie]) {
                    if (!std::isnan(ts)) h_all->Fill(ts - tl);
                }
            }

            // Imposta errori bin istogrammi
            for (int bin = 1; bin <= h_first->GetNbinsX(); ++bin) {
                double content = h_first->GetBinContent(bin);
                h_first->SetBinError(bin, std::sqrt(content));
            }
            for (int bin = 1; bin <= h_all->GetNbinsX(); ++bin) {
                double content = h_all->GetBinContent(bin);
                h_all->SetBinError(bin, std::sqrt(content));
            }

            res.dt_first_map[ch] = h_first;
            res.dt_all_map[ch] = h_all;
            res.h2_amp_vs_dt_map[ch] = h2_amp_vs_dt;
        }

        return res;
    }




};

#endif // TIMING_ANALYSIS_H
