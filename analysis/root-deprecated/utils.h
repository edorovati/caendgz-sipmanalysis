/******************************************************************************
 * @file Utils.h
 * @brief Declaration of the Utils class: static helper functions for file parsing,
 *        interpolation, graph creation, numerical analysis, and ROOT I/O.
 *
 * This header provides:
 *   - Utils::extractOvervoltage         : Extract overvoltage value from filename
 *   - Utils::interpolate                : Linear interpolation with error propagation
 *   - Utils::readDataFile               : Read two-column data from a text file
 *   - Utils::createGraph                : Build TGraphErrors from data and errors
 *   - Utils::findMaximum                : Find index of maximum y above a threshold
 *   - Utils::computeDerivative          : Compute discrete derivative dx/dy via midpoints
 *   - Utils::writeGraphTTree            : Write a TGraphErrors into the current ROOT file
 *   - Utils::findLocalMaxima            : Identify all local maxima in a dataset
 ******************************************************************************/
#ifndef UTILS_H
#define UTILS_H

#include <TString.h>
#include <TGraphErrors.h>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

/******************************************************************************
 *                                Utils Class
 * ----------------------------------------------------------------------------
 * @class Utils
 * @brief Collection of static utility functions for data I/O, interpolation,
 *        graph creation, and basic numerical operations in ROOT.
 *
 * This class includes the following public static functions:
 * ----------------------------------------------------------------------------
 * | **Public Functions**             | **Description**                                 |
 * ----------------------------------------------------------------------------
 * | **extractOvervoltage**           | Extracts overvoltage value from filename        |
 * ----------------------------------------------------------------------------
 * | **interpolate**                  | Performs linear interpolation with error
 *                                      propagation                                   |
 * ----------------------------------------------------------------------------
 * | **readDataFile**                 | Reads two-column data from a text file          |
 * ----------------------------------------------------------------------------
 * | **createGraph**                  | Builds a TGraphErrors from data and error vectors|
 * ----------------------------------------------------------------------------
 * | **findMaximum**                  | Finds index of maximum y above a threshold      |
 * ----------------------------------------------------------------------------
 * | **computeDerivative**            | Computes discrete derivative dy/dx
 *                                      with midpoints                                 |
 * ----------------------------------------------------------------------------
 * | **writeGraphTTree**              | Writes a TGraphErrors into the current ROOT file|
 * ----------------------------------------------------------------------------
 * | **findLocalMaxima**              | Identifies all local maxima in a y-vector       |
 * ----------------------------------------------------------------------------
 ******************************************************************************/
class Utils {
public:
    /**
     * @brief Extract overvoltage (V) from filename containing "vbias_".
     * @param filename  TString with pattern "vbias_<value>" (commas allowed).
     * @return Parsed voltage as double, or -1 on failure.
     */
    static double extractOvervoltage(const TString &filename) {
        TString name = filename;
        Ssiz_t start = name.Index("vbias_");
        if (start == kNPOS) return -1.0;
        start += 6;
        Ssiz_t end = start;
        while (end < name.Length() &&
               (isdigit(name[end]) || name[end] == '.' || name[end] == ',')) {
            ++end;
        }
        TString voltStr = name(start, end - start);
        voltStr.ReplaceAll(",", ".");
        return voltStr.Atof();
    }

    /**
     * @brief Perform linear interpolation between (x1,y1) and (x2,y2).
     * @param x       X-value at which to interpolate.
     * @param x1      First known x.
     * @param y1      First known y.
     * @param x2      Second known x.
     * @param y2      Second known y.
     * @param yErr1   Error on first y.
     * @param yErr2   Error on second y.
     * @return Pair (yInterp, yErrInterp) with rounded y and propagated error.
     */
    static std::pair<double, double>
    interpolate(double x,
                double x1, double y1,
                double x2, double y2,
                double yErr1, double yErr2) {
        if (x2 == x1) {
            return {std::round(y1), yErr1};
        }
        double w1 = (x2 - x) / (x2 - x1);
        double w2 = (x - x1) / (x2 - x1);
        double yInterp = w1 * y1 + w2 * y2;
        double errInterp = std::fabs(w1 * yErr1 + w2 * yErr2);
        return {std::round(yInterp), errInterp};
    }

    /**
     * @brief Read two-column data from a whitespace-delimited text file.
     * @param fullPath  Path to input file.
     * @param col1      Output vector for first column values.
     * @param col2      Output vector for second column values.
     * @return True if read successfully, false on error.
     */
    static bool readDataFile(const TString &fullPath,
                             std::vector<double> &col1,
                             std::vector<double> &col2) {
        std::ifstream infile(fullPath.Data());
        if (!infile) {
            std::cerr << "Error opening file: " << fullPath.Data() << std::endl;
            return false;
        }
        std::string line;
        std::getline(infile, line);  // skip header
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            double v1, v2;
            if (iss >> v1 >> v2) {
                col1.push_back(v1);
                col2.push_back(v2);
            } else {
                std::cerr << "Warning: malformed line skipped: " << line << std::endl;
            }
        }
        return true;
    }

    /**
     * @brief Build a TGraphErrors from vectors of values and errors.
     * @param xVals    X-coordinate values.
     * @param yVals    Y-coordinate values.
     * @param xErr     X errors.
     * @param yErr     Y errors.
     * @return Pointer to new TGraphErrors with points set.
     */
    static TGraphErrors* createGraph(const std::vector<double> &xVals,
                                     const std::vector<double> &yVals,
                                     const std::vector<double> &xErr,
                                     const std::vector<double> &yErr) {
        size_t n = xVals.size();
        TGraphErrors* g = new TGraphErrors(n);
        for (size_t i = 0; i < n; ++i) {
            g->SetPoint(i, xVals[i], yVals[i]);
            g->SetPointError(i, xErr[i], yErr[i]);
        }
        return g;
    }

    /**
     * @brief Find the index of the maximum y for x >= threshold.
     * @param x           X-values.
     * @param y           Y-values.
     * @param threshold   Minimum x to consider.
     * @return Index of max y above threshold, or -1 if none.
     */
    static int findMaximum(const std::vector<double> &x,
                           const std::vector<double> &y,
                           double threshold) {
        int idxMax = -1;
        double yMax = -INFINITY;
        for (size_t i = 0; i < x.size(); ++i) {
            if (x[i] >= threshold && y[i] > yMax) {
                yMax = y[i];
                idxMax = static_cast<int>(i);
            }
        }
        return idxMax;
    }

    /**
     * @brief Compute discrete derivative dy/dx with midpoints.
     * @param x        X-values (size N).
     * @param y        Y-values (size N).
     * @return Tuple of:
     *   - vector of midpoint x-values
     *   - vector of derivative values (dy/dx)
     */
    static std::tuple<std::vector<double>, std::vector<double>>
    computeDerivative(const std::vector<double> &x,
                      const std::vector<double> &y) {
        std::vector<double> xmid, dydx;
        size_t N = std::min(x.size(), y.size());
        for (size_t i = 0; i + 1 < N; ++i) {
            double dx = x[i+1] - x[i];
            double dy = y[i+1] - y[i];
            xmid.push_back(0.5*(x[i] + x[i+1]));
            dydx.push_back(-dy / dx);
        }
        return {xmid, dydx};
    }

    /**
     * @brief Write a TGraphErrors into the current ROOT file.
     * @param graph  Graph to write.
     * @param name   Name key in the TFile.
     */
    static void writeGraphTTree(TGraphErrors* graph,
                                const char* name) {
        graph->SetName(name);
        graph->Write();
    }

    /**
     * @brief Identify all local maxima in a y-vector.
     * @param y  Vector of values.
     * @return Vector of indices i where y[i-1] < y[i] > y[i+1].
     */
    static std::vector<int> findLocalMaxima(const std::vector<double> &y) {
        std::vector<int> peaks;
        for (size_t i = 1; i + 1 < y.size(); ++i) {
            if (y[i] > y[i-1] && y[i] > y[i+1]) {
                peaks.push_back(static_cast<int>(i));
            }
        }
        return peaks;
    }

  static Double_t* setupBranch_dgz(TTree* tree, const char* branchName) {
        Double_t* arr = new Double_t[1024];
        tree->SetBranchAddress(branchName, arr);
        return arr;
    }
    // -----------------------------------------------------------------------------
    //  Trova il picco "vero" della waveform su tutta la finestra temporale.
    //  Ritorna NaN se non trova cluster significativi.
    // -----------------------------------------------------------------------------
    static double find_signal_peak_anywhere(const double* wf,              // waveform corretta (1024 campioni)
                                    const double* time,            // vettore tempi (1024)
                                    int nsamples,
                                    double baseline_median,        // già calcolata
                                    double noise_rms,              // stimata sulla baseline
                                    double kSigma   = 5.0,         // soglia = median + k*sigma
                                    int    min_width = 3,          // # campioni consecutivi > soglia
                                    int    min_width50 = 2)        // larghezza a metà picco
    {
        const double threshold = baseline_median + kSigma * noise_rms;

        int j = 0;
        double best_peak = std::numeric_limits<double>::quiet_NaN();

        while (j < nsamples) {
            // Non controllo il tempo, cerco in tutta la waveform
            if (wf[j] <= threshold) { ++j; continue; }
            int cluster_start = j;
            double cluster_max = wf[j];

            // Avanza finché resta sopra soglia
            while (j < nsamples && wf[j] > threshold) {
                cluster_max = std::max(cluster_max, wf[j]);
                ++j;
            }
            int cluster_end = j;                // primo sample sotto soglia
            int width = cluster_end - cluster_start;

            if (width >= min_width) {
                // Misura larghezza a metà picco (grezzo)
                int left  = cluster_start;
                int right = cluster_end - 1;
                double half = cluster_max * 0.5;
                while (left  < cluster_end && wf[left]  < half) ++left;
                while (right > cluster_start && wf[right] < half) --right;
                int width50 = right - left + 1;

                if (width50 >= min_width50) {
                    // Cluster accettato: aggiorna picco migliore
                    if (std::isnan(best_peak) || cluster_max > best_peak)
                        best_peak = cluster_max;
                }
            }
        }
        return best_peak;
    }
};

#endif // UTILS_H
