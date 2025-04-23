#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TF1.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#include "utils.h"

class Analysis {
public:
    static const double wf;
    static const double timeWindow;
    static std::tuple<double, double, double, double, double> fitGaussian(TGraph* graph,
                                                                            double fitMin, double fitMax,
                                                                            double initialAmplitude, double initialMean, double initialSigma) {
        TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
        gausFit->SetParameters(initialAmplitude, initialMean, initialSigma);
        graph->Fit(gausFit, "R");
        std::cout << "Fit completato!" << std::endl;
        double mean = gausFit->GetParameter(1);
        double std = gausFit->GetParameter(2);
        double meanErr = gausFit->GetParError(1);
        double stdErr = gausFit->GetParError(2);
        double chi2 = gausFit->GetChisquare();
        double ndf = gausFit->GetNDF();
        double reducedChi2 = chi2 / ndf;
        return std::make_tuple(mean, std, meanErr, stdErr, reducedChi2);
    }

    static std::tuple<double, double> calculateCounts(const std::vector<double>& x, const std::vector<double>& y, double thr) {
        double counts = 0.0;
        double counts_errors = 0.0;

        for (int i = 0; i < x.size() - 1; ++i) {
            if (x[i] <= thr && x[i + 1] > thr) {
                std::tie(counts, counts_errors) = Utils::interpolate(
                    thr,
                    x[i], y[i],
                    x[i + 1], y[i + 1],
                    sqrt(y[i]), sqrt(y[i + 1])
                );
                break;
            }
        }

        
        return std::make_tuple(counts, counts_errors);
    }



    static std::tuple<double, double> calculateCTDCRRatio(double DCR_value, double DCR_error, double CT_value, double CT_error) {
        double CT_DCR = 0;
        double CT_DCR_error = 0;

        if (DCR_value != 0) {
            CT_DCR = (CT_value / DCR_value) * 100;
            CT_DCR_error = CT_DCR * std::sqrt(
                std::pow(CT_error / CT_value, 2) + std::pow(DCR_error / DCR_value, 2)
            );
        }

        return std::make_tuple(CT_DCR, CT_DCR_error);
    }
    
    static std::pair<double, double> calculateDarkRateAndError(double DCR, double DCR_Err, double wf, double timeWindow) {
        double darkRate = (DCR / (wf * timeWindow)) / 1000;
        double darkRateErr = (DCR_Err / DCR) * darkRate;
        return std::make_pair(darkRate, darkRateErr);
    }

};

#endif 
