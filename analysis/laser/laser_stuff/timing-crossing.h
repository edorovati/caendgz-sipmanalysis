#ifndef TIMING_PROCESSOR_H          // ─── inizio guardia ───
#define TIMING_PROCESSOR_H

#include <TString.h>
#include <TTree.h>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

/*========================================================================*\
|                       Classe utility per il timing                       |
\*========================================================================*/
class TimingProcessor {
public:
    // ---------------------------------------------------------------------
    //  1) Primo crossing di soglia per evento
    // ---------------------------------------------------------------------
    static std::map<TString, std::pair<Double_t*, Double_t*>> compute_threshold_crossings(
        const std::map<TString, Double_t*>& wf_map,      // waveform corrette
        const std::map<TString, Double_t*>& time_map,    // tempo campioni
        const std::map<TString, TTree*>&    tree_map,    // alberi canali
        double thr_mV)                                   // soglia in mV
    {
        const int nsamp = 1024;
        std::map<TString, std::pair<Double_t*, Double_t*>> out;

        for (const auto& [ch, wf] : wf_map) {
            if (!tree_map.count(ch)) continue;
            Long64_t N = tree_map.at(ch)->GetEntries();

            Double_t* t_arr = new Double_t[N];  // crossing times
            Double_t* a_arr = new Double_t[N];  // amplitudes associate

            for (Long64_t ie = 0; ie < N; ++ie) {
                int j_cross = -1;

                // Cerca crossing soglia (interpolato)
                for (int j = 1; j < nsamp; ++j) {
                    double v0 = wf[ie*nsamp + j - 1];
                    double v1 = wf[ie*nsamp + j];
                    if (v0 < thr_mV && v1 >= thr_mV) {
                        j_cross = j;
                        double t0 = time_map.at(ch)[j - 1];
                        double t1 = time_map.at(ch)[j];
                        t_arr[ie] = t0 + (thr_mV - v0) * (t1 - t0) / (v1 - v0);
                        break;
                    }
                }

                if (j_cross < 0) {
                    // Nessun crossing trovato
                    t_arr[ie] = std::numeric_limits<Double_t>::quiet_NaN();
                    a_arr[ie] = std::numeric_limits<Double_t>::quiet_NaN();
                    continue;
                }

                // Cerca massimo locale a partire da crossing
                int idx_start = j_cross;
                int idx_max = idx_start;
                double val_max = wf[ie*nsamp + idx_start];

                // Scorri in avanti per trovare massimo locale
                // Condizione: il valore deve crescere fino a max, poi decrescere senza risalire
                bool descending = false;

                for (int k = idx_start + 1; k < nsamp; ++k) {
                    double val = wf[ie*nsamp + k];
                    if (!descending) {
                        if (val > val_max) {
                            val_max = val;
                            idx_max = k;
                        } else if (val < val_max) {
                            descending = true; // ora deve solo scendere
                        }
                    } else {
                        // Se dopo discesa risale --> non è massimo valido
                        if (val > val_max) {
                            // Massimo non valido, abort
                            idx_max = -1;
                            break;
                        }
                    }
                }

                if (idx_max < 0) {
                    // Non trovato massimo valido
                    a_arr[ie] = std::numeric_limits<Double_t>::quiet_NaN();
                } else {
                    a_arr[ie] = wf[ie*nsamp + idx_max];
                }
            }

            out[ch] = {t_arr, a_arr};
        }

        return out;
    }

    // ---------------------------------------------------------------------
    //  2) Tutti i crossing sopra soglia per evento
    // ---------------------------------------------------------------------
    static std::map<TString, std::vector<std::vector<double>>> compute_all_crossings(
        const std::map<TString, Double_t*>& wf_map,
        const std::map<TString, Double_t*>& time_map,
        const std::map<TString, TTree*>&   tree_map,
        double thr_mV)
    {
        const int nsamp = 1024;
        std::map<TString, std::vector<std::vector<double>>> out;

        for (const auto& [ch, wf] : wf_map) {
            if (!tree_map.count(ch)) continue;
            Long64_t N = tree_map.at(ch)->GetEntries();
            std::vector<std::vector<double>> v_all(N);

            for (Long64_t ie = 0; ie < N; ++ie) {
                for (int j = 1; j < nsamp; ++j) {
                    double v0 = wf[ie*nsamp + j-1];
                    double v1 = wf[ie*nsamp + j];
                    if (v0 < thr_mV && v1 >= thr_mV) {
                        double t0 = time_map.at(ch)[j-1];
                        double t1 = time_map.at(ch)[j];
                        double tc = t0 + (thr_mV - v0)*(t1 - t0)/(v1 - v0);
                        v_all[ie].push_back(tc);
                    }
                }
            }
            out[ch] = std::move(v_all);
        }
        return out;
    }

    // ---------------------------------------------------------------------
    //  3) Half‑maximum crossing per la traccia laser
    // ---------------------------------------------------------------------
    static Double_t* compute_laser_crossings(TTree* tree_laser,
                                             Double_t* laser_amp,
                                             Double_t* laser_time)
    {
        if (!tree_laser || !laser_amp || !laser_time) return nullptr;

        const int nsamp = 1024;
        Long64_t  N     = tree_laser->GetEntries();
        auto* t_laser = new Double_t[N];

        for (Long64_t ie = 0; ie < N; ++ie) {
            tree_laser->GetEntry(ie);

            double vmin =  1e9, vmax = -1e9;
            for (int j = 0; j < nsamp; ++j) {
                double v = laser_amp[j];
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
            double thr = vmin + 0.5*(vmax - vmin);

            int j_cross = -1;
            for (int j = 1; j < nsamp; ++j) {
                if (laser_amp[j-1] < thr && laser_amp[j] >= thr) { j_cross = j; break; }
            }
            if (j_cross < 0) {
                t_laser[ie] = std::numeric_limits<Double_t>::quiet_NaN();
                continue;
            }

            double t0 = laser_time[j_cross-1];
            double t1 = laser_time[j_cross];
            double v0 = laser_amp [j_cross-1];
            double v1 = laser_amp [j_cross];
            t_laser[ie] = t0 + (thr - v0)*(t1 - t0)/(v1 - v0);
        }
        return t_laser;
    }
};

#endif /* TIMING_PROCESSOR_H */    // ─── fine guardia ───
