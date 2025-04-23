/******************************************************************************
 * @file Utils.h
 * @brief Declaration of the Utils class, offering static helper functions for
 *        file parsing, interpolation, graph creation, and basic numerical operations.
 *
 * This header defines:
 *   - Utils: routines for extracting overvoltage, interpolation, reading data
 *            files, creating TGraphErrors, finding maxima, computing derivatives,
 *            and writing graphs to ROOT TTrees.
 *
 ******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <TString.h>
#include <TGraphErrors.h>
#include <vector>
#include <utility>
#include <fstream>
#include <cmath>
#include <iostream>

/******************************************************************************
 *                                Utils Class
 * ----------------------------------------------------------------------------
 * @class Utils
 * @brief Collection of static utility functions for data I/O, interpolation,
 *        graph creation, and basic numerical operations.
 *
 * Public Functions:
 * ----------------------------------------------------------------------------
 * | Name                          | Description                                           |
 * ----------------------------------------------------------------------------
 * | extractOvervoltage            | Extract overvoltage value from a filename string      |
 * ----------------------------------------------------------------------------
 * | interpolate                   | Perform linear interpolation with error propagation   |
 * ----------------------------------------------------------------------------
 * | readDataFile                  | Read threshold/count pairs from a text file           |
 * ----------------------------------------------------------------------------
 * | createGraph                   | Build a TGraphErrors from vectors of x/y and errors   |
 * ----------------------------------------------------------------------------
 * | findMaximum                   | Find index of maximum y above a threshold             |
 * ----------------------------------------------------------------------------
 * | computeDerivativeWithErrors   | Compute derivative and error from discrete data       |
 * ----------------------------------------------------------------------------
 * | writeGraphTTree               | Write a TGraphErrors into a ROOT TTree                |
 * ----------------------------------------------------------------------------
 ******************************************************************************/

class Utils {
public:

// ------------------------------------------------------------------------------------ //
/// @brief  Extract the overvoltage (in volts) encoded in a filename.
/// @param  filename  TString containing an overvoltage label (e.g. "5,0OV.txt").
/// @return Overvoltage value as double (commas converted to periods).
// ------------------------------------------------------------------------------------ //
    static double extractOvervoltage(const TString &filename) {
        TString name = filename;
        name.ReplaceAll(".txt", "");
        name.ReplaceAll("OV", "");
        name.ReplaceAll(",", ".");
        return name.Atof();
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Perform linear interpolation between two points with propagated error.
/// @param  x      The x-value at which to interpolate.
/// @param  x1     First known x-value.
/// @param  y1     First known y-value.
/// @param  x2     Second known x-value.
/// @param  y2     Second known y-value.
/// @param  yErr1  Statistical error on y1.
/// @param  yErr2  Statistical error on y2.
/// @return Pair (y, yErr): interpolated y and its propagated error.
// ------------------------------------------------------------------------------------ //
    static std::pair<double, double>
    interpolate(double x,
                double x1, double y1,
                double x2, double y2,
                double yErr1, double yErr2)
    {
        if (x2 == x1) {
            return std::make_pair(y1, yErr1);
        }
        double w1 = (x2 - x) / (x2 - x1);
        double w2 = (x - x1) / (x2 - x1);
        double y  = w1 * y1 + w2 * y2;
        double yErr  = std::sqrt(std::pow(w1 * yErr1, 2) +
                                 std::pow(w2 * yErr2, 2));
        return std::make_pair(y, yErr);
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Read threshold and count pairs from a whitespace-delimited text file.
/// @param  fullPath   Path to the input text file.
/// @param  soglie     Output vector for threshold values.
/// @param  conteggi   Output vector for count values.
/// @return True if file was opened and read successfully; false otherwise.
// ------------------------------------------------------------------------------------ //
    static bool readDataFile(const TString &fullPath,
                             std::vector<double> &soglie,
                             std::vector<double> &conteggi)
    {
        std::ifstream infile(fullPath.Data());
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << fullPath.Data() << std::endl;
            return false;
        }
        double soglia;
        int conteggio;
        while (infile >> soglia >> conteggio) {
            soglie.push_back(soglia);
            conteggi.push_back(conteggio);
        }
        infile.close();
        return true;
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Construct a TGraphErrors from vectors of x, y, and their errors.
/// @param  xValues   Vector of x-coordinate values.
/// @param  yValues   Vector of y-coordinate values.
/// @param  xError    Vector of x-axis errors.
/// @param  yErrors   Vector of y-axis errors.
/// @return Pointer to a newly created TGraphErrors object.
// ------------------------------------------------------------------------------------ //
    static TGraphErrors*
    createGraph(const std::vector<double> &xValues,
                const std::vector<double> &yValues,
                const std::vector<double> &xError,
                const std::vector<double> &yErrors)
    {
        size_t N = xValues.size();
        TGraphErrors* graph = new TGraphErrors(N);
        for (size_t i = 0; i < N; ++i) {
            graph->SetPoint(i, xValues[i], yValues[i]);
            graph->SetPointError(i, xError[i], yErrors[i]);
        }
        return graph;
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Find the index of the maximum y-value above a given x-threshold.
/// @param  x          Vector of x-values.
/// @param  y          Vector of y-values.
/// @param  threshold  Minimum x-value to consider.
/// @return Index of the maximum y above threshold, or -1 if none found.
// ------------------------------------------------------------------------------------ //
    static int findMaximum(const std::vector<double> &x,
                           const std::vector<double> &y,
                           double threshold)
    {
        int peakIndex = -1;
        double maxVal = -1;
        for (size_t i = 0; i < x.size(); ++i) {
            if (x[i] >= threshold && y[i] > maxVal) {
                maxVal = y[i];
                peakIndex = static_cast<int>(i);
            }
        }
        return peakIndex;
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Compute the pointwise derivative of y w.r.t. x, with error propagation.
/// @param  x           Vector of x-values.
/// @param  y           Vector of y-values.
/// @param  y_err       Vector of y-value errors.
/// @param  correlation If true, errors are summed; otherwise added in quadrature.
/// @return Tuple of (x_midpoints, derivatives, derivativeErrors).
// ------------------------------------------------------------------------------------ //
    static std::tuple<std::vector<double>,
                      std::vector<double>,
                      std::vector<double>>
    computeDerivativeWithErrors(const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<double> &y_err,
                                bool correlation)
    {
        std::vector<double> x_mean, deriv, derivErr;
        size_t N = std::min({x.size(), y.size(), y_err.size()});
        for (size_t i = 0; i + 1 < N; ++i) {
            double dx = x[i+1] - x[i];
            double dy = y[i+1] - y[i];
            double derivVal = std::fabs(dy / dx);
            double err1 = y_err[i], err2 = y_err[i+1];
            double dyErr = correlation ? (err1 + err2)
                                       : std::sqrt(err1*err1 + err2*err2);
            double derivErrVal = std::fabs(dyErr / dx);
            x_mean.push_back((x[i] + x[i+1]) / 2.0);
            deriv.push_back(derivVal);
            derivErr.push_back(derivErrVal);
        }
        return std::make_tuple(x_mean, deriv, derivErr);
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Write a TGraphErrors into the current ROOT file as a TTree branch.
/// @param  graph  Pointer to the TGraphErrors to write.
/// @param  name   Name to assign to the graph object in the TFile.
// ------------------------------------------------------------------------------------
    static void writeGraphTTree(TGraphErrors* graph,
                                const char* name)
    {
        graph->SetName(name);
        graph->Write();
    }
// ------------------------------------------------------------------------------------ //
};

#endif
