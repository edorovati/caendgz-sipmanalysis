#ifndef FITTING_H
#define FITTING_H

#include <vector>
#include <tuple>
#include <cmath>
#include <TGraph.h>
#include <TF1.h>

class Fitting {
public:
    /**
     * Stima la costante di decadimento τ di una waveform media tramite fit esponenziale.
     *
     * @param times       Vettore dei tempi (ns)
     * @param amps        Vettore delle ampiezze (mV)
     * @param t_min_max   Tempo minimo per la ricerca del massimo (ns)
     * @param t_max_max   Tempo massimo per la ricerca del massimo (ns)
     * @return tuple<tau, tau_err, idx_max>
     */
    static std::tuple<double, double, size_t>
    estimateTau(const std::vector<double>& times,
                const std::vector<double>& amps,
                double t_min_max = 35.0,
                double t_max_max = 45.0)
    {
        double max_amp = -1e9;
        double t_max = 0.0;
        size_t idx_max = 0;

        // 🔍 1. Trova massimo nel range definito
        for (size_t i = 0; i < times.size(); ++i) {
            if (times[i] >= t_min_max && times[i] <= t_max_max && amps[i] > max_amp) {
                max_amp = amps[i];
                t_max   = times[i];
                idx_max = i;
            }
        }

        // 📉 2. Prepara dati decadimento fino al 5% del massimo
        std::vector<double> times_decay, amps_decay;
        double threshold = 0.05 * max_amp;

        for (size_t i = idx_max; i < times.size(); ++i) {
            if (amps[i] < threshold) break;  // fermati al 5% del massimo
            times_decay.push_back(times[i] - t_max); // shift massimo a zero
            amps_decay.push_back(amps[i]);
        }

        double tau = 0.0, tau_err = 0.0;

        // 📈 3. Fit esponenziale se abbastanza punti
        if (times_decay.size() > 3) {
            double amp_est = max_amp;
            double tau_est = 5.0;

            // Stima iniziale di tau
            size_t i1 = 0;
            while (i1 < times_decay.size() && amps_decay[i1] <= 0) i1++;
            size_t i2 = i1 + 1;
            while (i2 < times_decay.size() && amps_decay[i2] <= 0) i2++;
            if (i2 < times_decay.size()) {
                amp_est = amps_decay[i1];
                tau_est = (times_decay[i2] - times_decay[i1]) /
                          std::log(amps_decay[i1] / amps_decay[i2]);
            }

            TGraph decay_graph(times_decay.size(), &times_decay[0], &amps_decay[0]);
            TF1 decay_fit("decay_fit", "[0]*exp(-x/[1])", 0, times_decay.back());
            decay_fit.SetParameters(amp_est, tau_est);
            decay_fit.SetParNames("amp","tau");
            decay_graph.Fit(&decay_fit, "RQ");

            tau     = decay_fit.GetParameter(1);
            tau_err = decay_fit.GetParError(1);
        }

        return {tau, tau_err, idx_max};
    }

    
     // ==============================
    // Metodo 1: Calcolo slope t80-t20
    // ==============================
    // Modifica computeSlopeT80T20 per ritornare anche il fit
    static std::tuple<double, double, double>
    computeSlopeT80T20(const std::vector<double>& times,
                       const std::vector<double>& amps,
                       size_t idx_max)
    {
        const double ampMax = amps[idx_max];
        const double t_max  = times[idx_max];

        const double level80 = 0.8 * ampMax;
        const double level20 = 0.2 * ampMax;

        double t80 = NAN, t20 = NAN;
        int i80 = -1, i20 = -1;

        // Trova t80 andando indietro dal massimo
        for (int i = static_cast<int>(idx_max); i > 0; --i) {
            if (amps[i] >= level80 && amps[i-1] < level80) {
                double frac = (level80 - amps[i-1]) / (amps[i] - amps[i-1]);
                t80 = times[i-1] + frac * (times[i] - times[i-1]);
                i80 = i;
                break;
            }
        }
        if (i80 < 0) return {}; // non trovato

        // Trova t20 andando indietro da i80
        for (int i = i80 - 1; i > 0; --i) {
            if (amps[i] >= level20 && amps[i-1] < level20) {
                double frac = (level20 - amps[i-1]) / (amps[i] - amps[i-1]);
                t20 = times[i-1] + frac * (times[i] - times[i-1]);
                i20 = i;
                break;
            }
        }
        if (i20 < 0) return {}; // non trovato

        // Interpolazione fine tra t20 e t80 per fit lineare
        std::vector<double> tVals;
        std::vector<double> yVals;
        const int nInterp = 5;

        for (int i = i20; i <= i80; ++i) {
            double t1 = times[i-1], t2 = times[i];
            double y1 = amps[i-1],  y2 = amps[i];

            for (int j = 0; j <= nInterp; ++j) {
                double f = double(j) / nInterp;
                double t = t1 + f * (t2 - t1);
                double y = y1 + f * (y2 - y1);
                if (y >= level20 && y <= level80) {
                    tVals.push_back(t);
                    yVals.push_back(y);
                }
            }
        }
        if (tVals.size() < 2) return {};

        // Fit lineare con ROOT
        TGraph gr(tVals.size(), tVals.data(), yVals.data());
        TF1 fitFunc("fitFunc", "[0]*x + [1]", t20, t80);
        fitFunc.SetParameters((level80 - level20) / (t80 - t20), 0);

        if (gr.Fit(&fitFunc, "QNR") != 0) return {};

        double slope     = fitFunc.GetParameter(0);
        double slope_err = fitFunc.GetParError(0);
        double t0        = (slope != 0) ? -fitFunc.GetParameter(1) / slope : t20;

        return {slope, slope_err, t0};
    }


    static void gauss_fit(
    TGraphErrors* g,
    int peakIndex,
    double& mean, double& meanErr,
    double& sigma, double& sigmaErr)
    {
        if (!g || peakIndex < 0 || peakIndex >= g->GetN()) {
            std::cerr << "Indice picco non valido o grafico nullo." << std::endl;
            return;
        }
        
        // Ottieni punto di picco
        double peakX, peakY;
        g->GetPoint(peakIndex, peakX, peakY);
        double thresholdY = peakY * 0.5; // 5% del massimo
        
        // Trova estremi sinistro e destro
        int leftIndex = peakIndex;
        double xTmp, yTmp;
        while (leftIndex > 0) {
            g->GetPoint(leftIndex, xTmp, yTmp);
            if (yTmp <= thresholdY) break;
            --leftIndex;
        }
        
        int rightIndex = peakIndex;
        while (rightIndex < g->GetN() - 1) {
            g->GetPoint(rightIndex, xTmp, yTmp);
            if (yTmp <= thresholdY) break;
            ++rightIndex;
        }
        
        double fitMin, fitMax;
        g->GetPoint(leftIndex, fitMin, yTmp);
        g->GetPoint(rightIndex, fitMax, yTmp);
        
        std::cout << "Fit range: " << fitMin << " - " << fitMax << std::endl;
        
        // Costruisci sottografo con errori
        TGraphErrors gFitSegment;
        for (int i = leftIndex; i <= rightIndex; ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            double ex = g->GetErrorX(i);
            double ey = g->GetErrorY(i);
            gFitSegment.SetPoint(gFitSegment.GetN(), x, y);
            gFitSegment.SetPointError(gFitSegment.GetN()-1, ex, ey);
        }
        
        // Fit gaussiano
        TF1 gaus("gaus", "gaus", fitMin, fitMax);
        gaus.SetParameters(peakY, peakX, (fitMax - fitMin) / 6.0);
        gFitSegment.Fit(&gaus, "IMREQ"); // I: improve, M: better params, R: range, Q: quiet
        
        mean     = gaus.GetParameter(1);
        meanErr  = gaus.GetParError(1);
        sigma    = gaus.GetParameter(2);
        sigmaErr = gaus.GetParError(2);
    }
};

#endif // FITTING_H
