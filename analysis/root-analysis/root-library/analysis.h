#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TFile.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include "../root-library/graphics.h"

class Analysis {
public:
    struct SignalBackground {
        double S;
        double B;
        double diff;
    };

    static SignalBackground CalculateSignalAndBackground(TH1D* hist, int b) {
        int nBins = hist->GetNbinsX();
        auto* xaxis = hist->GetXaxis();

        // Valore centrale in tempo del bin b
        double center_time = xaxis->GetBinCenter(b);

        // Range in tempo (ns)
        double S_min = center_time - 2.0;
        double S_max = center_time + 2.0;

        double B1_min = center_time - 4.0;
        double B1_max = center_time - 2.0;

        double B2_min = center_time + 2.0;
        double B2_max = center_time + 4.0;

        // Converti i range temporali in bin usando FindBin
        int S_start = std::max(1, xaxis->FindBin(S_min));
        int S_end   = std::min(nBins, xaxis->FindBin(S_max));

        int B1_start = std::max(1, xaxis->FindBin(B1_min));
        int B1_end   = std::min(nBins, xaxis->FindBin(B1_max));

        int B2_start = std::max(1, xaxis->FindBin(B2_min));
        int B2_end   = std::min(nBins, xaxis->FindBin(B2_max));

        // Stampa info
        std::cout << "Bin centrale: " << b << " (posizione x: " << center_time << " ns)\n";
        std::cout << "Signal range time: [" << S_min << ", " << S_max << "] ns -> bins [" << S_start << ", " << S_end << "]\n";
        std::cout << "Background1 range time: [" << B1_min << ", " << B1_max << "] ns -> bins [" << B1_start << ", " << B1_end << "]\n";
        std::cout << "Background2 range time: [" << B2_min << ", " << B2_max << "] ns -> bins [" << B2_start << ", " << B2_end << "]\n";

        // Calcolo integrali
        double S = hist->Integral(S_start, S_end);
        double B = hist->Integral(B1_start, B1_end) + hist->Integral(B2_start, B2_end);

        std::cout << "Signal count: " << S << ", Background count: " << B << "\n";

        SignalBackground result = {S, B, S - B};
        return result;
    }
    
    // Calcola DCR (Dark Count Rate) e relativo errore
    static std::pair<double, double> CalculateDCR(double counts, double time_window, int n_wf) {
        if (time_window <= 0 || n_wf == 0) {
            std::cerr << "Errore: time_window <= 0 o n_wf == 0" <<std::endl;
            return {0., 0.};
        }
        double rate = counts / (time_window * n_wf);          // Hz
        double error = std::sqrt(counts) / (time_window * n_wf);
        rate *= 1e-3;  // kHz
        error *= 1e-3; // kHz
        return {rate, error};
    }
    
    // Calcola crosstalk e errore binomiale (in %)
    static std::pair<double, double> CalculateCrosstalk(double counts1,double counts2) {
        if (counts1 <= 0) {
            std::cerr << "Errore: counts1 <= 0 per calcolo crosstalk" <<std::endl;
            return {0., 0.};
        }
        double ct = counts1 / counts2;
        double error = std::sqrt(ct * (1 - ct) / counts2);
        return {ct * 100.0, error * 100.0};
    }
    
    // Calcola probabilità e errore binomiale
    static std::pair<double, double> CalculateProbability(double counts,int n_wf) {
        if (n_wf == 0) {
            std::cerr << "Errore: n_wf == 0 per calcolo probabilità" <<std::endl;
            return {0., 0.};
        }
        double p = counts / n_wf;
        double error = std::sqrt(p * (1 - p) / n_wf);
        return {p, error};
    }
};

#endif // ANALYSIS_H
