/******************************************************************************
 * @file compareCT.cpp
 * @brief Cross-talk (CT) comparison utilities for SiPM sensors.
 *
 * This script provides two main plotting utilities:
 *
 *  1. `compareCT`:
 *     - Plots Cross-Talk (CT) probability as a function of bias voltage.
 *
 *  2. `compareCTvsPE`:
 *     - Plots Cross-Talk (CT) probability as a function of single photoelectron
 *       (1 PE) amplitude.
 *     - Optionally includes horizontal error bands based on PE standard deviation.
 *
 * Usage:
 *   - Compile and run in ROOT environment:
 *       root -l
 *       .x compareCT.cpp
 *
 *   - OR in batch mode:
 *       root -l -b -q 'compareCT.cpp'
 *
 * Example:
 *   std::vector<std::pair<TString, TString>> files = {
 *     {"sn31-A1.root", "Commercial 50 um"},
 *     {"sn41-A1.root", "Custom 50 um"}
 *   };
 *
 *   compareCT(files);         // Plots CT vs Voltage
 *   compareCTvsPE(files);     // Plots CT vs PE amplitude
 *
 * Requirements:
 *   - ROOT files must contain a TTree named "Info" with relevant branches:
 *       For `compareCT`: voltage, CT, CT_Error
 *       For `compareCTvsPE`: pe, stdDev, CT, CT_Error
 *   - Custom headers (graphics, utils, etc.) must be accessible and compiled.
 *   - The data should be cleaned and calibrated prior to visualization.
 *
 * Notes:
 *   - Useful for evaluating CT trends across sensors and correlating CT
 *     with gain (via single p.e. amplitude).
 *   - Helps visualize systematic behavior of sensors as a function of operating voltage.
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
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

// ----------------------------------------------------------------------------
// @brief Plot CT (%) vs Bias Voltage for multiple sensors
// ----------------------------------------------------------------------------
void compareCT(const std::vector<std::pair<TString,TString>>& inputFiles) {
  TCanvas* c = GraphicsUtils::createCanvas("CompareCT", "CT vs Bias Voltage", 800,600, false);
  TLegend* legend = GraphicsUtils::createLegend(0.6, 0.7, 0.9, 0.9);

  for (size_t idx = 0; idx < inputFiles.size(); ++idx) {
    const auto& [filePath, sensorName] = inputFiles[idx];

    TTree* tree = nullptr;
    if (!RootUtils::loadTreeFromFile(filePath, "Info", tree)) {
      continue;
    }

    double vol, ct, ctErr;
    tree->SetBranchAddress("voltage", &vol);
    tree->SetBranchAddress("CT", &ct);
    tree->SetBranchAddress("CT_Error", &ctErr);

    auto* grCT = new TGraphErrors();
    int pt = 0;

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      grCT->SetPoint(pt, vol, ct);
      grCT->SetPointError(pt, 0, ctErr);
      ++pt;
    }

    Int_t color = ColorPaletteManager::getColorByIndex(idx);
    GraphicsUtils::setGraphStyle(grCT, 20, color);

    if (idx == 0) {
      grCT->Draw("AP");
      GraphicsUtils::setAxesBoldTitles(grCT, "Voltage (V)", "X-Talk probability (%)", true, 50.0, 60.0, 0.0, 100);
    } else {
      grCT->Draw("P SAME");
    }

    GraphicsUtils::addEntryToLegend(legend, grCT, sensorName.Data(), "lep");
  }

  legend->Draw();
}

// ----------------------------------------------------------------------------
// @brief Plot CT (%) vs Single p.e. amplitude for multiple sensors
// ----------------------------------------------------------------------------
void compareCTvsPE(const std::vector<std::pair<TString,TString>>& inputFiles) {
  TCanvas* c = GraphicsUtils::createCanvas("CompareCTvsPE", "CT vs Single p.e. amplitude", 800, 600, false);
  TLegend* legend = GraphicsUtils::createLegend(0.6, 0.7, 0.9, 0.9);

  for (size_t idx = 0; idx < inputFiles.size(); ++idx) {
    const auto& [filePath, sensorName] = inputFiles[idx];

    TTree* tree = nullptr;
    if (!RootUtils::loadTreeFromFile(filePath, "Info", tree)) {
      continue;
    }

    double pe, stdDev, ct, ctErr;
    tree->SetBranchAddress("pe", &pe);
    tree->SetBranchAddress("stdDev", &stdDev);
    tree->SetBranchAddress("CT", &ct);
    tree->SetBranchAddress("CT_Error", &ctErr);

    auto* grCTvsPE = new TGraphErrors();
    std::vector<BandPoint> horizontalBand;
    int pt = 0;

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      grCTvsPE->SetPoint(pt, pe, ct);
      grCTvsPE->SetPointError(pt, 0, ctErr);
      horizontalBand.push_back({ct, pe + stdDev, pe - stdDev});
      ++pt;
    }

    Int_t color = ColorPaletteManager::getColorByIndex(idx);
    GraphicsUtils::setGraphStyle(grCTvsPE, 20, color);

    std::sort(horizontalBand.begin(), horizontalBand.end(), [](const BandPoint& a, const BandPoint& b) {
      return a.x < b.x;
    });

    auto* bandGraph = GraphicsUtils::createHorizontalBandGraph(horizontalBand, color, 0.3, color, 2);

    if (idx == 0) {
      bandGraph->Draw("AF");
      grCTvsPE->Draw("P SAME");
      GraphicsUtils::setAxesBoldTitles(bandGraph, "Single p.e. amplitude (mV)", "X-Talk probability (%)",
                                       true, 0.0, 30.0, 0.0, 100);
    } else {
      bandGraph->Draw("F SAME");
      grCTvsPE->Draw("P SAME");
    }

    GraphicsUtils::addEntryToLegend(legend, grCTvsPE, sensorName.Data(), "lep");
  }

  legend->Draw();

  RootUtils::drawExplanationLegend();
}
