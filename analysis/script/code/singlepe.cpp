/******************************************************************************
 * @file compareSinglePE.cpp
 * @brief Visualization and linear fitting of single photoelectron (1 PE) amplitude vs. voltage.
 *
 * This script is used to compare the response of different sensors or measurement sets
 * by plotting the amplitude of the single photoelectron (1 PE) signal as a function of
 * applied voltage.
 *
 * Functionality:
 *   - Loads TTree objects from multiple ROOT files.
 *   - Extracts voltage, single PE amplitude, standard deviation, and error.
 *   - Builds graphs with error bars and shaded uncertainty bands.
 *   - Fits a custom linear function of the form:
 *         y = slope * (x - intercept)
 *     to the points after a given voltage threshold (`fitStartX`).
 *   - The intercept with the x-axis represents an *estimated breakdown voltage* (V_break).
 *   - Displays graphs and fits together with legends.
 *
 * Usage:
 *   - Open ROOT and run interactively, or from command line:
 *       root -l
 *       .x compareSinglePE.cpp
 *
 *   - OR from batch mode:
 *       root -l -b -q 'compareSinglePE.cpp'
 *
 *   - Inside the script (or your macro), call the function as:
 *
  *       std::vector<std::pair<TString, TString>> files = {
  *         {"sn31-A1.root", "Commercial 50 um"},
  *         {"sn41-A1.root", "Custom 50 um"}
  *       };
  *       compareSinglePE(files, 54);  // Start fit from 54 V
 *
 *   - Make sure:
 *       • The ROOT files exist in the current directory (or use full paths).
 *       • Each file contains a TTree named "Info" with branches:
 *           - voltage
 *           - pe (single PE amplitude)
 *           - stdDev (standard deviation for band)
 *           - peError (error bar)
 *       • All custom headers (analysis.h, graphics.h, etc.) are reachable by the compiler.
 *
 * Notes:
 *   - This analysis is useful for estimating roughly the breakdown voltage (V_break)
 *     by extrapolating the single PE amplitude to zero via linear fit.
 *   - It supports multiple sensors/data sets for visual comparison.
 ******************************************************************************/

#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/analysis.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/graphics.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/utils.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/colorpalettemanager.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/script/code/Helper_Function/RootUtils.h"


#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TLegend.h"
#include "TF1.h"
#include "TBox.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

// ----------------------------------------------------------------------------
// @brief Main function to compare single PE amplitude from multiple ROOT files
// ----------------------------------------------------------------------------
void compareSinglePE(const std::vector<std::pair<TString,TString>>& inputFiles, double fitStartX) {
  TCanvas* c = GraphicsUtils::createCanvas("ComparePE", "Single p.e. amplitude", 800,600, false);
  TLegend* legend = GraphicsUtils::createLegend(0.6,0.7,0.9,0.9);

  for (size_t idx = 0; idx < inputFiles.size(); ++idx) {
    const auto& [filePath, sensorName] = inputFiles[idx];

    TTree* tree = nullptr;
    if (!RootUtils::loadTreeFromFile(filePath, "Info", tree)) {
      continue;
    }

    double vol, pe, std, peErr;
    tree->SetBranchAddress("voltage", &vol);
    tree->SetBranchAddress("pe", &pe);
    tree->SetBranchAddress("stdDev", &std);
    tree->SetBranchAddress("peError", &peErr);

    std::vector<BandPoint> band;
    auto* grPE = new TGraphErrors();
    int pt = 0;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      grPE->SetPoint(pt, vol, pe);
      grPE->SetPointError(pt, 0, peErr);
      band.push_back({vol, pe + std, pe - std});
      ++pt;
    }

    Int_t color = ColorPaletteManager::getColorByIndex(idx);
    std::sort(band.begin(), band.end(), [](auto& a, auto& b){ return a.x < b.x; });
    auto* grBand = GraphicsUtils::createBandGraph(band, color, 0.3, color, 2);
    GraphicsUtils::setGraphStyle(grPE, 20, color);
    grBand->SetFillColorAlpha(color, 0.3);

    if (idx == 0) {
      grBand->Draw("AF");
      grPE->Draw("P SAME");
      GraphicsUtils::setAxesBoldTitles(grBand, "Voltage (V)", "Single p.e. (mV)", true, 50.0, 60.0, 0.0, 30.0);
    } else {
      grBand->Draw("F SAME");
      grPE->Draw("P SAME");
    }

    int nFitPoints = 0;
    double xMinFit = 1e6, xMaxFit = -1e6;
    for (int i = 0; i < grPE->GetN(); ++i) {
      double x, y;
      grPE->GetPoint(i, x, y);
      if (x >= fitStartX) {
        if (x < xMinFit) xMinFit = x;
        if (x > xMaxFit) xMaxFit = x;
        ++nFitPoints;
      }
    }

    if (nFitPoints >= 2) {
      TF1* fitFunc = new TF1(Form("fitFunc_%zu", idx), "[0]*(x - [1])", xMinFit, xMaxFit);
      fitFunc->SetLineColor(color);
      fitFunc->SetLineStyle(2);
      grPE->Fit(fitFunc, "QR");

      double slope     = fitFunc->GetParameter(0);
      double intercept = fitFunc->GetParameter(1);
      double slopeErr  = fitFunc->GetParError(0);
      double intErr    = fitFunc->GetParError(1);

      std::cout << "File: " << filePath << "\n";
      std::cout << "  Label: " << sensorName << "\n";
      std::cout << "  Fit: y = slope*(x - intercept)\n";
      std::cout << "    Slope     = " << slope << " ± " << slopeErr << "\n";
      std::cout << "    Intercept = " << intercept << " ± " << intErr << "\n\n";

      double xExtensionMin = xMaxFit - 10.0;
      double xExtensionMax = xMaxFit + 10.0;
      TF1* fitExtrapolated = new TF1(Form("fitExtrapolated_%zu", idx),
                                     Form("%f*(x - %f)", slope, intercept),
                                     xExtensionMin, xExtensionMax);
      fitExtrapolated->SetLineColor(color);
      fitExtrapolated->SetLineStyle(2);
      fitExtrapolated->Draw("SAME");

      fitFunc->Draw("SAME");
    } else {
      std::cerr << "WARNING: Not enough points >= " << fitStartX << " for fitting " << sensorName << "\n";
    }

    GraphicsUtils::addEntryToLegend(legend, grPE, sensorName.Data(), "lep");
  }

  legend->Draw();
  RootUtils::drawExplanationLegend();
}
