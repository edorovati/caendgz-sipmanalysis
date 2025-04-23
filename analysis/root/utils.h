#ifndef UTILS_H
#define UTILS_H

#include <TString.h>
#include <vector>
#include <utility>
#include <fstream>
#include <cmath>
#include <iostream>

class Utils {
public:
    static double extractOvervoltage(const TString &filename) {
        TString name = filename;
        name.ReplaceAll(".txt", "");
        name.ReplaceAll("OV", "");
        name.ReplaceAll(",", ".");
        return name.Atof();
    }

    static std::pair<double, double> interpolate(double x, double x1, double y1, double x2, double y2, double yErr1, double yErr2) {
        if (x2 == x1) {
            return std::make_pair(y1, yErr1);
        }

        double w1 = (x2 - x) / (x2 - x1);
        double w2 = (x - x1) / (x2 - x1);

        double y = w1 * y1 + w2 * y2;

        double yErr1Contribution = w1 * yErr1;
        double yErr2Contribution = w2 * yErr2;
        double yErr = std::sqrt(std::pow(yErr1Contribution, 2) + std::pow(yErr2Contribution, 2));

       
        return std::make_pair(y, yErr);
    }


    static bool readDataFile(const TString &fullPath, std::vector<double>& soglie, std::vector<double>& conteggi) {
        std::ifstream infile(fullPath.Data());
        if (!infile.is_open()) {
            std::cerr << "Errore nell'apertura del file " << fullPath.Data() << std::endl;
            return false;
        }
        std::string line;
        std::getline(infile, line);
        double soglia;
        int conteggio;
        while (infile >> soglia >> conteggio) {
            soglie.push_back(soglia);
            conteggi.push_back(conteggio);
        }
        infile.close();
        return true;
    }
    static TGraphErrors* createGraph(const std::vector<double>& xValues, const std::vector<double>& yValues, const std::vector<double>& xError, const std::vector<double>& yErrors) {
        size_t N = xValues.size();
        TGraphErrors* graph = new TGraphErrors(N);

        for (size_t i = 0; i < N; ++i) {
            graph->SetPoint(i, xValues[i], yValues[i]);
            graph->SetPointError(i, xError[i], yErrors[i]);
        }

        return graph;
    }
    
    static int findMaximum(const std::vector<double>& x, const std::vector<double>& y, double threshold) {
            int peakIndex = -1;
            double maxVal = -1;
            
            for (int i = 0; i < x.size(); ++i) {
                if (x[i] >= threshold && y[i] > maxVal) {
                    maxVal = y[i];
                    peakIndex = i;
                }
            }
            
            return peakIndex;
        }
    
    
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> computeDerivativeWithErrors(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& y_err,
        bool correlation)
    {
        std::vector<double> x_mean;
        std::vector<double> derivate;
        std::vector<double> derivateErr;

        int N = std::min({x.size(), y.size(), y_err.size()});

        for (int i = 0; i < N - 1; ++i) {
            double deltaS = x[i] - x[i + 1];
            double deltaC = y[i] - y[i + 1];

            double derivata = std::fabs(deltaC / deltaS);

            double err1 = y_err[i];
            double err2 = y_err[i + 1];

            double deltaC_err;
            if (correlation) {
                deltaC_err = err1 + err2;
            } else {
                deltaC_err = std::sqrt(err1 * err1 + err2 * err2);
            }

            double derivataErr = std::fabs(deltaC_err / deltaS);
            double xMedia = (x[i] + x[i + 1]) / 2.0;

            x_mean.push_back(xMedia);
            derivate.push_back(derivata);
            derivateErr.push_back(derivataErr);
        }

        return std::make_tuple(x_mean, derivate, derivateErr);
    }
    
    static void writeGraphTTree(TGraphErrors* graph, const char* name) {
            graph->SetName(name);
            graph->Write();
        }



    };
#endif 
