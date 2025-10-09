/******************************************************************************
 * @file Analysis.h
 * @brief Declaration of the Analysis class: static methods for waveform and
 *        statistical analyses in ROOT, including Gaussian fits, threshold
 *        interpolation, CT/DCR ratio, and dark-rate estimation.
 *
 * This header provides:
 *   - Analysis::fitGaussian            : Gaussian fit on TGraph
 *   - Analysis::calculateCounts        : Linear interpolation at threshold
 *   - Analysis::calculateCTDCRRatio    : Cross-talk to dark-count ratio
 *   - Analysis::calculateDarkRateAndError : Dark rate (kHz) with error
 ******************************************************************************/
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "utils.h"
#include <TF1.h>
#include <TGraph.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

/******************************************************************************
 *                                Analysis Class
 * ----------------------------------------------------------------------------
 * @class Analysis
 * @brief Provides a collection of static methods for waveform and statistical
 *        analysis in ROOT, such as Gaussian fitting, threshold counting, and
 *        DCR/CT ratio evaluation.
 *
 * This class includes the following public static functions:
 * ----------------------------------------------------------------------------
 * | **Public Functions**                   | Description                         |
 * ----------------------------------------------------------------------------
 * | **fitGaussian**                        | Performs a Gaussian fit on a graph  |
 *                                          | and returns fit parameters and chi² |
 * ----------------------------------------------------------------------------
 * | **calculateCounts**                    | Finds threshold crossing point in   |
 *                                          | histogram and interpolates count    |
 * ----------------------------------------------------------------------------
 * | **calculateCTDCRRatio**                | Computes ratio between CT and DCR   |
 *                                          | with propagated uncertainty         |
 * ----------------------------------------------------------------------------
 * | **calculateDarkRateAndError**          | Calculates the dark count rate      |
 *                                          | normalized per ns and its error     |
 * ----------------------------------------------------------------------------
 ******************************************************************************/
class Analysis {
public:
    /// Default number of waveforms (used in dark-rate calculation)
    static constexpr double defaultNumWaveforms = 1.0;
    /// Default time window per waveform [s]
    static constexpr double defaultTimeWindow = 1.0;

    /**
     * @brief Fit a Gaussian function to a TGraph within [fitMin, fitMax].
     *
     * This method fits the function "gaus" (ROOT builtin) to the provided graph,
     * using the specified initial parameters, and returns the fit results.
     *
     * @param graph             Pointer to the TGraph to be fitted.
     * @param fitMin            Lower x-bound for the fit window.
     * @param fitMax            Upper x-bound for the fit window.
     * @param initialAmplitude  Initial guess for the Gaussian amplitude.
     * @param initialMean       Initial guess for the Gaussian mean.
     * @param initialSigma      Initial guess for the Gaussian sigma (standard deviation).
     * @return A tuple containing:
     *   - mean      : Fitted Gaussian mean
     *   - sigma     : Fitted Gaussian sigma
     *   - meanErr   : Uncertainty on the mean
     *   - sigmaErr  : Uncertainty on the sigma
     *   - redChi2   : Reduced chi-square of the fit
     */
    static std::tuple<double, double, double, double, double>
    fitGaussian(TGraph* graph,
                 double fitMin,
                 double fitMax,
                 double initialAmplitude,
                 double initialMean,
                 double initialSigma) {
        TF1* gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
        gausFit->SetParameters(initialAmplitude, initialMean, initialSigma);
        graph->Fit(gausFit, "R+Q");

        double mean     = gausFit->GetParameter(1);
        double sigma    = gausFit->GetParameter(2);
        double meanErr  = gausFit->GetParError(1);
        double sigmaErr = gausFit->GetParError(2);
        double chi2     = gausFit->GetChisquare();
        double ndf      = gausFit->GetNDF();
        double redChi2  = (ndf > 0) ? (chi2 / ndf) : -1.0;

        return std::make_tuple(mean, sigma, meanErr, sigmaErr, redChi2);
    }

    /**
     * @brief Linearly interpolate counts at a given threshold.
     *
     * Scans the x-vector for the interval containing thr and performs
     * linear interpolation of y, propagating statistical errors.
     *
     * @param x    Vector of x-values (monotonic).
     * @param y    Corresponding y-values (counts).
     * @param thr  Threshold x-position for interpolation.
     * @return A tuple:
     *   - counts    : Interpolated y-value at thr
     *   - countsErr : Propagated statistical error
     */
    static std::tuple<double, double>
    calculateCounts(const std::vector<double>& x,
                    const std::vector<double>& y,
                    double thr) {
        double counts = 0.0;
        double countsErr = 0.0;
        for (size_t i = 0; i + 1 < x.size(); ++i) {
            if (x[i] <= thr && thr <= x[i+1]) {
                double err1 = std::sqrt(y[i]);
                double err2 = std::sqrt(y[i+1]);
                std::tie(counts, countsErr) =
                    Utils::interpolate(thr,
                                       x[i], y[i], err1,
                                       x[i+1], y[i+1], err2);
                break;
            }
        }
        return std::make_tuple(counts, countsErr);
    }

    /**
     * @brief Compute CT to DCR ratio (in percent) with error propagation.
     *
     * Efficiency = CT / DCR, expressed as percentage. Propagated uncertainty
     * uses binomial approximation: sqrt(eff*(1-eff)/DCR).
     *
     * @param DCR_value  Measured dark count rate.
     * @param DCR_error  Uncertainty on the dark count rate.
     * @param CT_value   Measured cross-talk count.
     * @param CT_error   Uncertainty on the cross-talk count.
     * @return A tuple:
     *   - ratio    : CT/DCR ratio [%]
     *   - ratioErr : Propagated error on ratio
     */
    static std::tuple<double, double>
    calculateCTDCRRatio(double DCR_value,
                        double DCR_error,
                        double CT_value,
                        double CT_error) {
        double ratio = 0.0;
        double ratioErr = 0.0;

        if (DCR_value > 0 && CT_value <= DCR_value) {
            double eff = CT_value / DCR_value;
            ratio = eff * 100.0;
            ratioErr = 100.0 * std::sqrt(eff * (1.0 - eff) / DCR_value);
        } else if (CT_value > DCR_value) {
            std::cerr << "[Analysis] Warning: CT_value > DCR_value. "
                      << "Unphysical efficiency>1." << std::endl;
        }
        return std::make_tuple(ratio, ratioErr);
    }

    /**
     * @brief Calculate dark count rate (kHz) and its uncertainty.
     *
     * Converts raw DCR count and error into rate per kHz, normalized by
     * number of waveforms and time window.
     *
     * @param DCR         Raw dark count.
     * @param DCR_Err     Uncertainty on the raw count.
     * @param numWaveforms Number of waveforms (triggers).
     * @param timeWindow  Duration of each waveform [s].
     * @return A pair:
     *   - darkRate    : Rate in kHz
     *   - darkRateErr : Propagated uncertainty
     */
    static std::pair<double, double>
    calculateDarkRateAndError(double DCR,
                               double DCR_Err,
                               double numWaveforms = defaultNumWaveforms,
                               double timeWindow = defaultTimeWindow) {
        double normFactor = 1.0 / (numWaveforms * timeWindow * 1e3);
        double darkRate    = DCR * normFactor;
        double darkRateErr = DCR_Err * normFactor;
        return std::make_pair(darkRate, darkRateErr);
    }
};

#endif // ANALYSIS_H
