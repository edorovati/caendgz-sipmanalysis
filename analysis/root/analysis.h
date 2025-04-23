/******************************************************************************
 * @file Analysis.h
 * @brief Declaration of the Analysis class, offering a collection of static
 *        methods for common waveform and statistical analyses within ROOT.
 *
 * This header defines:
 *   - Analysis: static routines for Gaussian fitting, threshold-based counting,
 *               cross-talk/DCR ratio calculation, and dark-rate estimation.
 *
 ******************************************************************************/

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "utils.h"
#include <TF1.h>
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
  static const double wf;
  static const double timeWindow;
  static std::tuple<double, double, double, double, double>
    
// ------------------------------------------------------------------------------------ //
/// @brief  Fit a Gaussian to the provided graph over [fitMin, fitMax].
/// @param  graph             Pointer to the TGraph to fit.
/// @param  fitMin            Lower bound of the fit interval.
/// @param  fitMax            Upper bound of the fit interval.
/// @param  initialAmplitude  Initial guess for Gaussian amplitude.
/// @param  initialMean       Initial guess for Gaussian mean.
/// @param  initialSigma      Initial guess for Gaussian sigma.
/// @return Tuple of (mean, sigma, meanErr, sigmaErr, redChi2):
      ///         • mean      : Fitted mean
      ///         • sigma     : Fitted standard deviation
      ///         • meanErr   : Uncertainty on the mean
      ///         • sigmaErr  : Uncertainty on the sigma
      ///         • redChi2   : Reduced χ² of the fit
// ------------------------------------------------------------------------------------ //
  fitGaussian(TGraph *graph, double fitMin, double fitMax,
              double initialAmplitude, double initialMean,
              double initialSigma) {
    TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
    gausFit->SetParameters(initialAmplitude, initialMean, initialSigma);
    graph->Fit(gausFit, "R");
    double mean = gausFit->GetParameter(1);
    double std = gausFit->GetParameter(2);
    double meanErr = gausFit->GetParError(1);
    double stdErr = gausFit->GetParError(2);
    double chi2 = gausFit->GetChisquare();
    double ndf = gausFit->GetNDF();
    double reducedChi2 = chi2 / ndf;
    return std::make_tuple(mean, std, meanErr, stdErr, reducedChi2);
  }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Interpolate the number of counts at a specified threshold.
/// @param  x    Vector of x-values (e.g., bin boundaries or thresholds).
/// @param  y    Vector of y-values (e.g., counts per bin).
/// @param  thr  Threshold at which to interpolate count.
/// @return Tuple of (counts, countsErr):
     ///         • counts    : Interpolated count at thr
     ///         • countsErr : Statistical error on the interpolated count
// ------------------------------------------------------------------------------------ //
  static std::tuple<double, double>
  calculateCounts(const std::vector<double> &x, const std::vector<double> &y,
                  double thr) {
    double counts = 0.0;
    double counts_errors = 0.0;

    for (int i = 0; i < x.size() - 1; ++i) {
      if (x[i] <= thr && x[i + 1] > thr) {
        std::tie(counts, counts_errors) = Utils::interpolate(
            thr, x[i], y[i], x[i + 1], y[i + 1], sqrt(y[i]), sqrt(y[i + 1]));
        break;
      }
    }

    return std::make_tuple(counts, counts_errors);
  }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Compute the cross-talk (CT) to dark count rate (DCR) ratio in percent.
/// @param  DCR_value  Measured dark count rate.
/// @param  DCR_error  Uncertainty on the dark count rate.
/// @param  CT_value   Measured cross-talk count.
/// @param  CT_error   Uncertainty on the cross-talk count.
/// @return Tuple of (ratio, ratioErr):
      ///         • ratio    : CT/DCR ratio [%]
      ///         • ratioErr : Propagated error on the ratio
// ------------------------------------------------------------------------------------ //
  static std::tuple<double, double> calculateCTDCRRatio(double DCR_value,
                                                        double DCR_error,
                                                        double CT_value,
                                                        double CT_error) {
    double CT_DCR = 0;
    double CT_DCR_error = 0;

    if (DCR_value != 0) {
      CT_DCR = (CT_value / DCR_value) * 100;
      CT_DCR_error = CT_DCR * std::sqrt(std::pow(CT_error / CT_value, 2) +
                                        std::pow(DCR_error / DCR_value, 2));
    }

    return std::make_tuple(CT_DCR, CT_DCR_error);
  }
// ------------------------------------------------------------------------------------ //
    
// ------------------------------------------------------------------------------------ //
/// @brief  Calculate the dark rate (in kHz) and its propagated error.
/// @param  DCR         Raw dark count value.
/// @param  DCR_Err     Uncertainty on the dark count.
/// @param  wf          Number of waveforms (or triggers).
/// @param  timeWindow  Time window per waveform [s].
/// @return Pair of (darkRate, darkRateErr):
      ///         • darkRate    : Dark rate in kHz
      ///         • darkRateErr : Propagated uncertainty on dark rate
// ------------------------------------------------------------------------------------ //
  static std::pair<double, double>
  calculateDarkRateAndError(double DCR, double DCR_Err, double wf,
                            double timeWindow) {
    double darkRate = (DCR / (wf * timeWindow)) / 1000;
    double darkRateErr = (DCR_Err / DCR) * darkRate;
    return std::make_pair(darkRate, darkRateErr);
  }
};
// ------------------------------------------------------------------------------------ //

#endif
