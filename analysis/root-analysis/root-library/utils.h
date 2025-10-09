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
#include <TH2D.h>
#include <TH1D.h>

class Utils {
public:
    static Double_t* setupBranch_dgz(TTree* tree, const char* branchName) {
        Double_t* arr = new Double_t[1024];
        tree->SetBranchAddress(branchName, arr);
        return arr;
    }
    
    static int FindMaxBin(TH1D* hist) {
        int maxBin = hist->GetMaximumBin();
        return maxBin;
    }
    
    static double CountEvents(TH1D* hist, int bin_start, int bin_end) {
        int nBins = hist->GetNbinsX();
        bin_start = std::max(1, std::min(bin_start, nBins));
        bin_end = std::max(1, std::min(bin_end, nBins));
        return hist->Integral(bin_start, bin_end);
    }
    
    static double median(std::vector<double> vec) {
        if (vec.empty()) return 0.0;
        std::sort(vec.begin(), vec.end());
        size_t n = vec.size();
        if (n % 2 == 0)
            return 0.5 * (vec[n/2 - 1] + vec[n/2]);
        else
            return vec[n/2];
    }
    static double median_error(const std::vector<double>& vec) {
        if (vec.empty()) return 0.0;

        // Calcola la mediana
        std::vector<double> sorted = vec;
        std::sort(sorted.begin(), sorted.end());
        size_t n = sorted.size();
        double med = (n % 2 == 0) ? 0.5 * (sorted[n/2 - 1] + sorted[n/2]) : sorted[n/2];

        // Calcola deviazione standard del campione
        double sum_sq = 0.0;
        for (double v : vec) sum_sq += (v - med) * (v - med);
        double sigma = std::sqrt(sum_sq / (n - 1));

        // Errore standard della mediana
        double sigma_med = 1.253 * sigma / std::sqrt(n);
        return sigma_med;
    }
    
    static TH1D* projectionY(TH2D* h2, double lower_x, double upper_x, const char* name = "_proj") {
        if (!h2) return nullptr;

        // Troviamo i bin corrispondenti ai limiti in X
        int binx1 = h2->GetXaxis()->FindBin(lower_x);
        int binx2 = h2->GetXaxis()->FindBin(upper_x);

        // Proiettiamo solo i bin selezionati su Y
        TH1D* hProj = h2->ProjectionY(name, binx1, binx2);

        return hProj;
    }
    
    static TH1D* projectionX(TH2D* h2, double lower_y, double upper_y, const char* name = "_projX") {
        if (!h2) return nullptr;

        // Troviamo i bin corrispondenti ai limiti in Y
        int biny1 = h2->GetYaxis()->FindBin(lower_y);
        int biny2 = h2->GetYaxis()->FindBin(upper_y);

        // Proiettiamo solo i bin selezionati su X
        TH1D* hProj = h2->ProjectionX(name, biny1, biny2);

        return hProj;
    }



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
    
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    computeDerivative(const std::vector<double> &x,
                      const std::vector<double> &y) {
        std::vector<double> xmid, dydx, err;
        size_t N = std::min(x.size(), y.size());
        for (size_t i = 0; i + 1 < N; ++i) {
            double dx = x[i+1] - x[i];
            double dy = y[i+1] - y[i];
            xmid.push_back((x[i] + x[i+1]) / 2.0 + dx); // solo la media
            dydx.push_back(-dy / dx);

            // Calcolo errore come nel tuo codice, evita sqrt di numero negativo
            double e = std::sqrt(std::max(0.0, -dy / dx));
            err.push_back(e);
        }
        return {xmid, dydx, err}; // <- qui serve il punto e virgola
    }

    
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

    
};

#endif // UTILS_H
