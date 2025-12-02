#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include "TTree.h"
#include "timestuff.h"
#include "utils.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>

class Waveform {
public:
    Waveform(TTree* laserTree, Double_t* laserTime, Double_t* laserAmp, double threshold = 70.0)
        : laser_tree(laserTree), laser_time(laserTime), laser_amp(laserAmp), threshold(threshold)
    {
        waveform_size = 1024;
        computeLaserCrossings();
    }

    // Funzione 1: Allineamento
    std::vector<std::pair<std::vector<double>, std::vector<double>>> getWaveforms(
        TTree* tree, Double_t* time, Double_t* amp, bool align = true) const
    {
        int n_events = std::min((int)laser_tree->GetEntries(), (int)tree->GetEntries());
        std::vector<std::pair<std::vector<double>, std::vector<double>>> waveforms;
        for (int i = 0; i < n_events; ++i) {
            tree->GetEntry(i);
            std::vector<double> t_vec(waveform_size);
            std::vector<double> a_vec(waveform_size);
            for (int j = 0; j < waveform_size; ++j) {
                t_vec[j] = time[j];
                if (align)
                    t_vec[j] -= laser_crossings[i];
                a_vec[j] = amp[j];
            }
            waveforms.emplace_back(std::move(t_vec), std::move(a_vec));
        }
        return waveforms;
    }

   

    std::vector<std::pair<std::vector<double>, std::vector<double>>> correctWaveforms(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
    double pre_signal_end = 30.0, size_t bunch_size = 1024
                                                                                      ) {
        std::vector<std::pair<std::vector<double>, std::vector<double>>> corrected_waveforms;
        size_t n_wf = waveforms.size();
        size_t n_bunches = (n_wf + bunch_size - 1) / bunch_size;
        
        for (size_t b = 0; b < n_bunches; ++b) {
            size_t start_idx = b * bunch_size;
            size_t end_idx = std::min(start_idx + bunch_size, n_wf);
            
            if (end_idx - start_idx <= 1) continue; // niente da fare se solo la wf da scartare
            
            // --- media del bunch (escludendo la prima waveform) ---
            std::vector<double> mean_amp = waveforms[start_idx + 1].second;
            std::vector<double> time_ref = waveforms[start_idx + 1].first;
            
            for (size_t wf_idx = start_idx + 2; wf_idx < end_idx; ++wf_idx) {
                const auto& amp = waveforms[wf_idx].second;
                for (size_t i = 0; i < mean_amp.size(); ++i)
                    mean_amp[i] += amp[i];
            }
            
            for (auto& a : mean_amp) a /= double(end_idx - start_idx - 1);
            
            // --- fit pre-signal con pol0 ---
            std::vector<double> pre_time, pre_amp;
            for (size_t i = 0; i < time_ref.size(); ++i)
                if (time_ref[i] <= pre_signal_end) {
                    pre_time.push_back(time_ref[i]);
                    pre_amp.push_back(mean_amp[i]);
                }
            
            TGraph g(pre_time.size(), pre_time.data(), pre_amp.data());
            TF1 fit("fit", "[0]", 0, pre_signal_end);
            g.Fit(&fit, "Q"); // Q per silenzioso
            double offset = fit.GetParameter(0);
            
            // --- correggi tutte le waveform del bunch (tranne la prima) ---
            for (size_t wf_idx = start_idx + 1; wf_idx < end_idx; ++wf_idx) {
                const auto& [time, amp] = waveforms[wf_idx];
                std::vector<double> corrected_amp = amp;
                for (auto& a : corrected_amp) a -= offset;
                corrected_waveforms.emplace_back(time, corrected_amp);
            }
        }
        
        return corrected_waveforms;
    }

    // Funzione 3: calibrazione
    static std::tuple<std::unique_ptr<TH1F>, std::unique_ptr<TH1F>, std::unique_ptr<TH1F>>
    calibration(
        const TString& ch_name,
        const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
        double t_min,
        double t_max,
        int n_bins,
        double lower_range,
        double upper_range,
        bool apply_threshold = false,
        double amp_threshold = 0.0
    ) {
        std::vector<double> max_signal_vals;
        std::vector<double> max_noise_vals;

        max_signal_vals.reserve(waveforms.size());
        max_noise_vals.reserve(waveforms.size());

        // durata delle finestre
        double dt_signal = t_max - t_min;
        double dt_total  = waveforms.front().first.back() - waveforms.front().first.front();
        double dt_noise  = dt_total - dt_signal;

        for (const auto& wf : waveforms) {
            const auto& times = wf.first;
            const auto& amps  = wf.second;

            double max_signal = std::numeric_limits<double>::lowest();
            double max_noise  = std::numeric_limits<double>::lowest();

            bool found_signal = false;
            bool found_noise  = false;

            for (size_t i = 0; i < times.size(); ++i) {
                if (times[i] >= t_min && times[i] <= t_max) {
                    max_signal = std::max(max_signal, amps[i]);
                    found_signal = true;
                } else {
                    max_noise = std::max(max_noise, amps[i]);
                    found_noise = true;
                }
            }

            if (found_signal && (!apply_threshold || max_signal > amp_threshold)) {
                max_signal_vals.push_back(max_signal);
            }
            if (found_noise) {
                max_noise_vals.push_back(max_noise);
            }
        }

        if (max_signal_vals.empty() || max_noise_vals.empty()) {
            return {nullptr, nullptr, nullptr};
        }

        // Istogrammi
        auto hist_signal = std::make_unique<TH1F>((ch_name + "_signal").Data(),
                                                   (ch_name + "_signal").Data(),
                                                   n_bins, lower_range, upper_range);
        auto hist_noise  = std::make_unique<TH1F>((ch_name + "_noise").Data(),
                                                  (ch_name + "_noise").Data(),
                                                  n_bins, lower_range, upper_range);

        for (double v : max_signal_vals) hist_signal->Fill(v);
        for (double v : max_noise_vals)  hist_noise->Fill(v);

        // Normalizzazione per la lunghezza delle finestre
        
        // Costruzione istogramma differenza: clone e sottrazione
        auto hist_diff = std::make_unique<TH1F>(*hist_signal); // clone del segnale
        hist_diff->SetName((ch_name + "_diff").Data());
        hist_diff->SetTitle((ch_name + "_diff").Data());
        hist_diff->Add(hist_noise.get(), -1.0); // sottraggo il rumore

        return {std::move(hist_noise), std::move(hist_signal), std::move(hist_diff)};
    }



    struct AvgWaveform {
        std::vector<double> times;
        std::vector<double> amplitudes;
        std::vector<double> errors;
    };

    // Calcola la waveform media (baseline corretta) da un insieme di waveform selezionate
    enum class StackMethod { Median = 0, SigmaClip = 1, WeightedNoise = 2, TrimmedMean = 3 };

    static double median_of(std::vector<double> v) {
        if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
        std::sort(v.begin(), v.end());
        size_t n = v.size();
        return (n % 2) ? v[n/2] : 0.5*(v[n/2-1] + v[n/2]);
    }

    static double mad_of(const std::vector<double>& v, double med) {
        std::vector<double> d; d.reserve(v.size());
        for (double x : v) d.push_back(std::fabs(x - med));
        return median_of(d);
    }

    // linear interpolation of single point t onto (tsrc, ysrc), returns NaN if out-of-range
    static double interp_point(const std::vector<double>& tsrc, const std::vector<double>& ysrc, double t) {
        if (tsrc.size() < 2) return std::numeric_limits<double>::quiet_NaN();
        if (t <= tsrc.front()) return ysrc.front();
        if (t >= tsrc.back()) return ysrc.back();
        auto it = std::upper_bound(tsrc.begin(), tsrc.end(), t);
        size_t j = std::distance(tsrc.begin(), it);
        if (j == 0) return ysrc.front();
        if (j >= tsrc.size()) return ysrc.back();
        size_t j0 = j - 1, j1 = j;
        double t0 = tsrc[j0], t1 = tsrc[j1];
        double y0 = ysrc[j0], y1 = ysrc[j1];
        double alpha = (t - t0) / (t1 - t0);
        return y0 + alpha * (y1 - y0);
    }
    struct WaveformInfo {
        std::vector<double> times;
        std::vector<double> amplitudes;
        double crossing_time;
        double max_amplitude;
    };
    /*
     single unified function:
     - selected: vector of pairs (times, values) per event
     - method: stacking method
     - tmin/tmax/dt: grid (if dt <= 0, derived from first event)
     - clip_sigma/clip_iters: for SigmaClip
     - trim_frac: for TrimmedMean (fraction trimmed from each tail, 0..0.49)
    */
    static AvgWaveform computeAverageWaveform(
        const std::vector<WaveformInfo>& selected,
        StackMethod method = StackMethod::Median,
        double tmin = -1000.0, double tmax = 1000.0, double dt = -1.0,
        double clip_sigma = 3.0, int clip_iters = 3,
        double trim_frac = 0.1)
    {
        AvgWaveform result;
        if (selected.empty()) return result;

        // infer dt if not provided (use median of diffs of first event if available)
        if (dt <= 0.0) {
            const auto& ts0 = selected.front().times;
            if (ts0.size() >= 2) {
                std::vector<double> diffs;
                diffs.reserve(ts0.size()-1);
                for (size_t i = 1; i < ts0.size(); ++i) diffs.push_back(ts0[i] - ts0[i-1]);
                dt = median_of(diffs);
                if (!(dt > 0.0)) dt = 1.0; // fallback
            } else {
                dt = 1.0;
            }
        }

        // build common grid of centers
        std::vector<double> tgrid;
        for (double t = tmin; t <= tmax + 1e-12; t += dt) tgrid.push_back(t + dt * 0.5);

        struct Event {
            std::vector<double> interp; // one value per grid point (NaN if originally out-of-range)
            double noise_var; // estimated from baseline region (sigma^2)
        };
        std::vector<Event> events;
        events.reserve(selected.size());

        for (const auto& wf : selected) {
            const auto& ts = wf.times;
            const auto& ys = wf.amplitudes;
            if (ts.size() != ys.size() || ts.size() < 2) continue;

            // baseline: median of ys where ts <= 35 ns
            std::vector<double> base_pts;
            for (size_t i = 0; i < ts.size(); ++i) if (ts[i] <= 35.0) base_pts.push_back(ys[i]);
            double base_med = base_pts.empty() ? 0.0 : median_of(base_pts);

            // estimate baseline noise variance via MAD -> sigma
            double noise_var = 0.0;
            if (!base_pts.empty()) {
                double med = median_of(base_pts);
                double mad = mad_of(base_pts, med);
                double sigma = (mad > 0.0) ? 1.4826 * mad : 0.0;
                noise_var = sigma * sigma;
            }

            // interpolate onto grid (centers)
            std::vector<double> interp;
            interp.reserve(tgrid.size());
            for (double tc : tgrid) {
                double v = interp_point(ts, ys, tc);
                if (std::isnan(v)) interp.push_back(std::numeric_limits<double>::quiet_NaN());
                else interp.push_back(v - base_med); // baseline-correct here
            }

            events.push_back(Event{std::move(interp), noise_var});
        }

        if (events.empty()) return result;

        // for each time bin, collect values and compute stack
        for (size_t i = 0; i < tgrid.size(); ++i) {
            std::vector<double> vals;
            vals.reserve(events.size());
            std::vector<double> wts; wts.reserve(events.size());
            for (const auto& ev : events) {
                double v = ev.interp[i];
                if (!std::isnan(v)) {
                    vals.push_back(v);
                    double w = (ev.noise_var > 0.0) ? (1.0 / ev.noise_var) : 1.0;
                    wts.push_back(w);
                }
            }

            if (vals.empty()) continue;
            double amp = 0.0;
            double err = 0.0;
            size_t N = vals.size();

            if (method == StackMethod::Median) {
                double med = median_of(vals);
                double mad = mad_of(vals, med);
                double sigma = 1.4826 * mad;
                amp = med;
                err = sigma / std::sqrt((double)N);
            }
            else if (method == StackMethod::SigmaClip) {
                std::vector<double> cur = vals;
                for (int it = 0; it < clip_iters && cur.size() >= 2; ++it) {
                    double mean = std::accumulate(cur.begin(), cur.end(), 0.0) / cur.size();
                    double var = 0.0;
                    for (double x : cur) var += (x - mean)*(x - mean);
                    var /= (cur.size() > 1 ? (cur.size()-1) : 1);
                    double sd = std::sqrt(var);
                    std::vector<double> next; next.reserve(cur.size());
                    for (double x : cur) if (std::fabs(x - mean) <= clip_sigma * sd) next.push_back(x);
                    if (next.size() == cur.size()) break;
                    cur.swap(next);
                }
                double mean = std::accumulate(cur.begin(), cur.end(), 0.0) / cur.size();
                double var = 0.0;
                for (double x : cur) var += (x - mean)*(x - mean);
                var /= (cur.size() > 1 ? (cur.size()-1) : 1);
                amp = mean;
                err = std::sqrt(var / (double)cur.size());
            }
            else if (method == StackMethod::WeightedNoise) {
                double sw = 0.0, swx = 0.0;
                for (size_t k = 0; k < vals.size(); ++k) {
                    double w = wts[k];
                    sw += w;
                    swx += w * vals[k];
                }
                amp = (sw > 0.0) ? (swx / sw) : 0.0;
                err = (sw > 0.0) ? std::sqrt(1.0 / sw) : 0.0;
            }
            else if (method == StackMethod::TrimmedMean) {
                std::sort(vals.begin(), vals.end());
                double f = std::clamp(trim_frac, 0.0, 0.49);
                size_t trim = static_cast<size_t>(std::floor(f * vals.size()));
                size_t lo = trim, hi = vals.size() - trim;
                if (lo >= hi) {
                    amp = median_of(vals);
                    err = 0.0;
                } else {
                    double s = 0.0;
                    for (size_t k = lo; k < hi; ++k) s += vals[k];
                    double mean = s / (double)(hi - lo);
                    double var = 0.0;
                    for (size_t k = lo; k < hi; ++k) var += (vals[k] - mean)*(vals[k] - mean);
                    var /= ((hi - lo) > 1 ? (hi - lo - 1) : 1);
                    amp = mean;
                    err = std::sqrt(var / (double)(hi - lo));
                }
            }

            result.times.push_back(tgrid[i]);
            result.amplitudes.push_back(amp);
            result.errors.push_back(err);
        }

        return result;
    }

    
   

    static std::vector<WaveformInfo>
    selection(const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
              double lower_thr, double upper_thr, double thr,
              double& crossing_time_avg, double& crossing_sigma)
    {
        std::vector<std::pair<std::vector<double>, std::vector<double>>> preselected;
        std::vector<double> pre_crossing_times;

        // ---- PRE-SELEZIONE ----
        for (const auto& wf_raw : waveforms) {
            const auto& times = wf_raw.first;
            const auto& amps_raw = wf_raw.second;
            if (amps_raw.empty()) continue;

            // Baseline
            std::vector<double> baseline_window;
            for (size_t i = 0; i < times.size(); ++i) {
                if (times[i] <= 30.0) baseline_window.push_back(amps_raw[i]);
                else break;
            }
            if (baseline_window.empty()) continue;

            std::nth_element(baseline_window.begin(),
                             baseline_window.begin() + baseline_window.size()/2,
                             baseline_window.end());
            double baseline_median = baseline_window[baseline_window.size()/2];

            std::vector<double> amps = amps_raw;
            for (auto& a : amps) a -= baseline_median;

            // Max in range 38-50 ns
            double maxAmpRange = -1e9;
            for (size_t i = 0; i < times.size(); ++i) {
                if (times[i] >= 38.0 && times[i] <= 50.0) {
                    if (amps[i] > maxAmpRange) maxAmpRange = amps[i];
                }
            }

            if (maxAmpRange >= lower_thr && maxAmpRange <= upper_thr) {
                preselected.emplace_back(times, amps);
                pre_crossing_times.push_back(maxAmpRange);
            }
        }

        if (preselected.empty()) return {};

        // ---- Crossing e massimo reale ----
        std::vector<WaveformInfo> selected;
        std::vector<double> crossing_times;
        for (const auto& wf : preselected) {
            const auto& times = wf.first;
            const auto& amps  = wf.second;

            TimeStuff::CrossingResult result = TimeStuff::analyzeWaveform(times.data(), amps.data(),
                                                                          static_cast<int>(times.size()), thr);
            // Se crossing_time < 0 → waveform scartata dalla funzione stessa
            if (result.crossing_time < 0) continue;

            // Qui sappiamo che non ci sono altri massimi più grandi dopo max_index
            double maxAfterCrossing = amps[result.max_index];

            // Controlla che il massimo sia nel range desiderato
            if (maxAfterCrossing < lower_thr || maxAfterCrossing > upper_thr) continue;

            // Se passa il controllo → seleziona la waveform
            selected.push_back({times, amps, result.crossing_time, maxAfterCrossing});
            crossing_times.push_back(result.crossing_time);
        }

        if (selected.empty()) return {};

        // ---- Istogramma, bin più probabile, media e sigma ----
        double bin_width = 0.1;
        std::map<int,int> hist;
        for (double t : crossing_times) hist[static_cast<int>(t/bin_width)]++;

        int max_bin = 0, max_count = 0;
        for (auto [b,c] : hist) if (c>max_count) { max_count=c; max_bin=b; }

        double bin_center = max_bin * bin_width;
        std::vector<double> local_crossings;
        for (double t : crossing_times)
            if (t >= bin_center - 3*bin_width && t <= bin_center + 3*bin_width)
                local_crossings.push_back(t);

        if (local_crossings.empty()) return {};

        crossing_time_avg = std::accumulate(local_crossings.begin(), local_crossings.end(), 0.0) / local_crossings.size();
        double sqsum=0;
        for (double t : local_crossings) sqsum += (t-crossing_time_avg)*(t-crossing_time_avg);
        crossing_sigma = std::sqrt(sqsum/local_crossings.size());

        // ---- Filtro finale su crossing window ----
        std::vector<WaveformInfo> final_selected;
        for (size_t i=0;i<selected.size();++i)
            if (selected[i].crossing_time >= crossing_time_avg - 4 &&
                selected[i].crossing_time <= crossing_time_avg + 4)
                final_selected.push_back(selected[i]);

        return final_selected;
    }


    

    static std::pair<double, double> integral_waveform_with_error(
        const std::vector<std::pair<std::vector<double>, std::vector<double>>>& selected_wfs,
        double /*crossing_time*/)   // crossing_time non serve più
    {
        std::vector<double> gains;

        for (const auto& wf : selected_wfs) {
            const auto& times = wf.first;
            const auto& amps  = wf.second;

            if (times.empty() || amps.empty()) continue;

            // ====== Integrazione con trapezi su tutto l’intervallo ======
            double integral_wf = 0.0;
            for (size_t i = 0; i + 1 < times.size(); ++i) {
                double t0 = times[i],   t1 = times[i + 1];
                double a0 = amps[i];
                double a1 = amps[i + 1];
                double dt = t1 - t0;

                integral_wf += 0.5 * (a0 + a1) * dt;
            }

            gains.push_back(integral_wf);
        }

        // ====== Calcolo media ed errore standard ======
        double mean_gain = 0.0, sum_gain = 0.0, sum_gain_sq = 0.0;
        int n = gains.size();

        for (auto g : gains) {
            sum_gain += g;
            sum_gain_sq += g * g;
        }

        if (n > 0) mean_gain = sum_gain / n;

        double stddev_gain = 0.0;
        if (n > 1) {
            double variance = (sum_gain_sq / n) - (mean_gain * mean_gain);
            if (variance > 0) stddev_gain = std::sqrt(variance);
        }

        double err_mean_gain = (n > 0) ? (stddev_gain / std::sqrt(n)) : 0.0;

        return {mean_gain, err_mean_gain};
    }

    static std::pair<double, double> estimate_tau_with_error(
        const std::vector<std::pair<std::vector<double>, std::vector<double>>>& selected_wfs)
    {
        std::vector<double> taus;

        for (const auto& wf : selected_wfs) {
            const auto& times = wf.first;
            const auto& amps  = wf.second;

            if (times.empty() || amps.empty()) continue;

            // ====== Trova massimo tra 38 e 42 ns ======
            double max_amp = -1e9;
            size_t idx_max = 0;
            for (size_t i = 0; i < times.size(); ++i) {
                if (times[i] >= 38.0 && times[i] <= 42.0) {
                    if (amps[i] > max_amp) {
                        max_amp = amps[i];
                        idx_max = i;
                    }
                }
            }

            if (max_amp < 0) continue; // nessun massimo trovato

            // ====== Trova punto dove l'ampiezza scende a 37% del massimo ======
            double target = 0.37 * max_amp;
            size_t idx_tau = idx_max;

            for (size_t i = idx_max; i < times.size(); ++i) {
                if (amps[i] <= target) {
                    idx_tau = i;
                    break;
                }
            }

            // ====== Calcola tau come differenza di tempo ======
            double tau = times[idx_tau] - times[idx_max];
            if (tau > 0) taus.push_back(tau);
        }

        // ====== Calcola media ed errore standard ======
        double mean_tau = 0.0, sum_tau = 0.0, sum_tau_sq = 0.0;
        int n = taus.size();

        for (auto t : taus) {
            sum_tau += t;
            sum_tau_sq += t * t;
        }

        if (n > 0) mean_tau = sum_tau / n;

        double stddev_tau = 0.0;
        if (n > 1) {
            double variance = (sum_tau_sq / n) - (mean_tau * mean_tau);
            if (variance > 0) stddev_tau = std::sqrt(variance);
        }

        double err_mean_tau = (n > 0) ? (stddev_tau / std::sqrt(n)) : 0.0;

        return {mean_tau, err_mean_tau};
    }


    
    // ====================================
    // Metodo 2: Integrale da t0 a 1% maxAmp
    // ====================================
    static double computeIntegralFromT0Minus5ToEnd(const std::vector<double>& times,
                                                   const std::vector<double>& amps,
                                                   double t0)
    {
        double integral = 0.0;
        bool start_integrate = false;

        for (size_t i = 1; i < times.size(); ++i) {  // partiamo da i=1 perché usiamo i-1
            if (!start_integrate && times[i] >= t0 - 5.0) {
                start_integrate = true;
            }
            if (start_integrate) {
                double dt = times[i] - times[i - 1];
                integral += 0.5 * (amps[i] + amps[i - 1]) * dt;  // trapezoidale
            }
        }

        return integral;
    }



    // ====================================
    // Metodo 3: Distribuzione ampiezze <35ns
    // ====================================
    static void fillAmplitudeDistributionBefore35ns(const std::vector<double>& times,
                                             const std::vector<double>& amps,
                                             TH1D* hist)
    {
        if (!hist) return;
        for (size_t i = 0; i < times.size(); ++i) {
            if (times[i] <= 35.0) {
                hist->Fill(amps[i]);
            } else {
                break; // tempi ordinati
            }
        }
    }
    
    static std::vector<std::pair<float, float>> get_transitions(const TGraph& graph, float threshold, float sign)
    {
        std::vector<std::pair<float, float>> values;
        bool armed = false;

        for (int i = 0; i < graph.GetN(); ++i) {
            double x = graph.GetX()[i];
            double y = graph.GetY()[i];

            if (!armed && sign * y > sign * threshold * 0.5)
                continue;

            armed = true;

            if (sign * y < sign * threshold)
                continue;

            values.emplace_back(x, y);
            armed = false;
        }

        return values;
    }


    
private:
    TTree* laser_tree;
    Double_t* laser_time;
    Double_t* laser_amp;
    double threshold;
    int waveform_size;
    std::vector<double> laser_crossings;

    void computeLaserCrossings() {
        int n_events = laser_tree ? laser_tree->GetEntries() : 0;
        laser_crossings.clear();
        laser_crossings.reserve(n_events);

        for (int i = 0; i < n_events; ++i) {
            laser_tree->GetEntry(i);
            double t_laser = TimeStuff::computeCrossingTime(laser_time, laser_amp, waveform_size, threshold);
            laser_crossings.push_back(t_laser);
        }
    }
};

#endif // WAVEFORM_H
