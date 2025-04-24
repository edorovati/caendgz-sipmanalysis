/******************************************************************************
 * @file dcr_analysis.cpp
 * @brief A script for DCR analysis, including waveform processing, Gaussian fitting,
 *        threshold counting, and DCR/CT ratio evaluation for different overvoltage levels.
 *
 * This script performs the following operations:
 *   - Reads data files containing waveform information.
 *   - Computes derivatives of the data to identify peaks.
 *   - Fits a Gaussian to the data and calculates the dark count rate (DCR) and
 *     cross-talk (CT) ratio.
 *   - Dynamically assigns colors to graphs based on overvoltage values.
 *   - Outputs the results to a ROOT file, including graphs and information about
 *     overvoltages, counts, errors, DCR, and CT ratios.
 *
 ******************************************************************************/

#include "../../root/analysis.h"
#include "../../root/graphics.h"
#include "../../root/utils.h"
#include "../../root/colorpalettemanager.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TRandom.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

// ----------------------------------------------------------------------------
// @brief Structure to hold legend entry data
// ----------------------------------------------------------------------------
struct LegendEntry {
  double overvoltage;
  TString label;
  int color;
  TGraphErrors *graph;
};
std::vector<LegendEntry> legendEntries;

// ----------------------------------------------------------------------------
// @brief Main function for DCR analysis
// ----------------------------------------------------------------------------
void dcr() {
  // Create output ROOT file
  TFile *outFile = new TFile("output.root", "RECREATE");

  // Set up directory for file scanning
  TSystemDirectory dir("Scan", "Scan");
  TMultiGraph *mgOrig = new TMultiGraph();  // MultiGraph for original data
  std::vector<double> allOvervoltages;     // List of all overvoltages

  // Create output directory and tree for storing results
  outFile->mkdir("Summary");
  outFile->cd("Summary");
  TTree *infoTree = new TTree("Info", "Information on overvoltages [current cycle]");
  
  // Define variables to store information for each overvoltage
  double ov, pe, std, peErr, stdErr, DCR, CT, reducedChi2;
  double DCR_Error, CT_Error;
  int colorIndex = 0;

  // Create branches for storing data in the tree
  infoTree->Branch("overvoltage", &ov, "overvoltage/D");
  infoTree->Branch("pe", &pe, "pe/D");
  infoTree->Branch("stdDev", &std, "stdDev/D");
  infoTree->Branch("peError", &peErr, "peError/D");
  infoTree->Branch("stdDevError", &stdErr, "stdDevError/D");
  infoTree->Branch("DCR", &DCR, "DCR/D");
  infoTree->Branch("CT", &CT, "CT/D");
  infoTree->Branch("DCR_Error", &DCR_Error, "DCR_Error/D");
  infoTree->Branch("CT_Error", &CT_Error, "CT_Error/D");

  // Scan files in the "Scan" directory
  TList *files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "No files found in the Scan folder!" << std::endl;
    return;
  }

  // Loop over all files in the directory
  TSystemFile *fileItem = nullptr;
  TIter next(files);
  while ((fileItem = (TSystemFile *)next())) {
    TString fname = fileItem->GetName();
    if (fileItem->IsDirectory() || fname.BeginsWith(".") || fname == "0OV.txt")
      continue;
    if (!fname.Contains("OV") || !fname.EndsWith(".txt"))
      continue;

    std::vector<double> thresholds;     // Threshold values
    std::vector<double> counts;   // Count values

    // Read data from the file
    if (!Utils::readDataFile(TString("Scan/") + fname, thresholds, counts))
      continue;
    
    int N = thresholds.size();
    if (N == 0)
      continue;

    // Extract overvoltage from the filename
    double overvoltage = Utils::extractOvervoltage(fname);
    if (overvoltage != 0.0) {
      allOvervoltages.push_back(overvoltage);
    }

    // Calculate error values
    std::vector<double> yErrors;
    for (size_t i = 0; i < counts.size(); ++i) {
      yErrors.push_back(sqrt(counts[i]));
    }
    
    // Define x-error as a fixed value
    std::vector<double> xErrors(thresholds.size(), 0.25);
    
    // Create graphs for raw data and derivatives
    TGraphErrors *gRawOrig = Utils::createGraph(thresholds, counts, xErrors, yErrors);
    auto [sMedie, derivate, derivateErr] = Utils::computeDerivativeWithErrors(thresholds, counts, yErrors, false);
    std::vector<double> xErrors_deriv(sMedie.size(), 0.0);
    TGraphErrors *gDerivOrig = Utils::createGraph(sMedie, derivate, xErrors_deriv, derivateErr);

    // Find peak index in the derivative
    int peakIndex = Utils::findMaximum(sMedie, derivate, 6.0);

    if (peakIndex != -1) {
      // Fit Gaussian to the derivative and calculate DCR and CT ratio
      std::tie(pe, std, peErr, stdErr, reducedChi2) = Analysis::fitGaussian(
          gDerivOrig, sMedie[peakIndex] - 4.0, sMedie[peakIndex] + 4.0,
          derivate[peakIndex], sMedie[peakIndex], 1.0);
      std::tie(DCR, DCR_Error) = Analysis::calculateCounts(thresholds, counts, 0.5 * pe);
      auto [CT_value, CT_error] = Analysis::calculateCounts(thresholds, counts, 1.5 * pe);
      std::tie(CT, CT_Error) = Analysis::calculateCTDCRRatio(DCR, DCR_Error, CT_value, CT_error);
      
      // Store data in the info tree
      ov = overvoltage;
      infoTree->Fill();
    }

    // Assign color dynamically based on overvoltage value
    Int_t dynamicColor = ColorPaletteManager::getColorFromOvervoltage(overvoltage);
    GraphicsUtils::setGraphStyle(gRawOrig, 20, dynamicColor, 1);
    ++colorIndex;

    // Set style for the derivative graph
    GraphicsUtils::setGraphStyle(gDerivOrig, 20, gRawOrig->GetLineColor(), 1);

    // Create subdirectory for storing graphs
    TString dirName = fname;
    dirName.ReplaceAll(".txt", "");
    TDirectory *subDir = outFile->mkdir(dirName.Data());
    subDir->cd();

    // Write graphs to ROOT file
    Utils::writeGraphTTree(gRawOrig, "gRaw");
    Utils::writeGraphTTree(gDerivOrig, "gDeriv");
    outFile->cd();

    // Create legend entry
    TString legendLabel;
    legendLabel.Form("@ %.1f OV", overvoltage);
    gRawOrig->SetTitle(legendLabel.Data());
    mgOrig->Add(gRawOrig, "LPE");

    // Add entry to legend
    LegendEntry entry;
    entry.overvoltage = overvoltage;
    entry.label = legendLabel;
    entry.color = dynamicColor;
    entry.graph = gRawOrig;
    legendEntries.push_back(entry);
  }

  // Create and style canvas
  auto cOrig = GraphicsUtils::createCanvas("cOrig", "all curves", 800, 600, true);
  GraphicsUtils::setAxisAndDraw(mgOrig, "Threshold (mV)", "Counts", true);

  // Sort legend entries by overvoltage
  std::sort(legendEntries.begin(), legendEntries.end(),
            [](const LegendEntry &a, const LegendEntry &b) {
              return a.overvoltage < b.overvoltage;
            });

  // Create legend and draw it
  auto legOrig = GraphicsUtils::createLegend(0.7, 0.7, 0.89, 0.89);
  for (const auto &entry : legendEntries) {
    legOrig->AddEntry(entry.graph, entry.label.Data(), "pe");
  }
  legOrig->Draw();

  // Write canvas to file
  cOrig->Write();
  outFile->Write();
  outFile->Close();
}
