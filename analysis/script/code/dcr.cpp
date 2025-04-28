/******************************************************************************
 * @file dcr_analysis.cpp
 * @brief A script for DCR analysis, including waveform processing, Gaussian fitting,
 *        threshold counting, and DCR/CT ratio evaluation for different overvoltage levels.
 *
 * This script performs the following operations:
 *   - Reads data files containing waveform information.
 *   - Subtracts 0OV baseline measurements.
 *   - Computes derivatives of the data to identify peaks.
 *   - Fits a Gaussian to the data and calculates the dark count rate (DCR) and
 *     cross-talk (CT) ratio.
 *   - Dynamically assigns colors to graphs based on overvoltage values.
 *   - Outputs the results to a ROOT file, including graphs and information about
 *     overvoltages, counts, errors, DCR, and CT ratios.
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
  TMultiGraph *mgOrig = new TMultiGraph();
  std::vector<double> allOvervoltages;

  // Create output directory and tree
  outFile->mkdir("Summary");
  outFile->cd("Summary");

  TTree *infoTree = new TTree("Info", "Information on overvoltages [current cycle]");

  // Variables for TTree
  double ov, pe, std, peErr, stdErr, DCR, CT, reducedChi2;
  double DCR_Error, CT_Error;

  // Create branches
  infoTree->Branch("overvoltage", &ov, "overvoltage/D");
  infoTree->Branch("pe", &pe, "pe/D");
  infoTree->Branch("stdDev", &std, "stdDev/D");
  infoTree->Branch("peError", &peErr, "peError/D");
  infoTree->Branch("stdDevError", &stdErr, "stdDevError/D");
  infoTree->Branch("DCR", &DCR, "DCR/D");
  infoTree->Branch("CT", &CT, "CT/D");
  infoTree->Branch("DCR_Error", &DCR_Error, "DCR_Error/D");
  infoTree->Branch("CT_Error", &CT_Error, "CT_Error/D");

  // --- First, read the 0OV baseline ---
  std::vector<double> thresholds0V, counts0V;
  bool baselineAvailable = Utils::readDataFile("Scan/0OV.txt", thresholds0V, counts0V);

  if (!baselineAvailable) {
    std::cerr << "WARNING: 0OV.txt baseline not found! Continuing without baseline subtraction." << std::endl;
  }

  std::vector<double> errors0V;
  if (baselineAvailable) {
    for (auto c : counts0V) {
      errors0V.push_back(sqrt(c));
    }
  }

  // --- Now, read all OV files ---
  TList *files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "No files found in the Scan folder!" << std::endl;
    return;
  }

  TSystemFile *fileItem = nullptr;
  TIter next(files);
  while ((fileItem = (TSystemFile *)next())) {
    TString fname = fileItem->GetName();

    if (fileItem->IsDirectory() || fname.BeginsWith("."))
      continue;
    if (!fname.Contains("OV") || !fname.EndsWith(".txt"))
      continue;
    if (fname == "0OV.txt")  // Skip baseline (already loaded)
      continue;

    std::vector<double> thresholds, counts;
    if (!Utils::readDataFile("Scan/" + fname, thresholds, counts))
      continue;

    if (thresholds.empty())
      continue;

    // Extract overvoltage value
    double overvoltage = Utils::extractOvervoltage(fname);
    allOvervoltages.push_back(overvoltage);

    // Prepare error arrays
    std::vector<double> yErrors(counts.size());
    for (size_t i = 0; i < counts.size(); ++i) {
      yErrors[i] = sqrt(counts[i]);
    }

    // --- Subtract 0V baseline if available ---
    if (baselineAvailable && thresholds.size() == thresholds0V.size()) {
      for (size_t i = 0; i < counts.size(); ++i) {
        counts[i] -= counts0V[i];
        yErrors[i] = sqrt(pow(yErrors[i], 2) + pow(errors0V[i], 2));
      }
    }

    // Fixed X errors
    std::vector<double> xErrors(thresholds.size(), 0.25);

    // Create graphs
    TGraphErrors *gRawOrig = Utils::createGraph(thresholds, counts, xErrors, yErrors);

    // Derivative and derivative graph
    auto [sMedie, derivate, derivateErr] = Utils::computeDerivativeWithErrors(thresholds, counts, yErrors, false);
    TGraphErrors *gDerivOrig = Utils::createGraph(sMedie, derivate, std::vector<double>(sMedie.size(), 0.0), derivateErr);
      
      // --- Find all local maximums ---
      std::vector<int> peakIndices = Utils::findLocalMaxima(derivate);


      // Checking that there are maximum
      if (!peakIndices.empty()) {
          // Find the maximum with the smallest X value => first spad
          int firstPeakIndex = peakIndices[0];
          for (int idx : peakIndices) {
              if (sMedie[idx] < sMedie[firstPeakIndex]) {
                  firstPeakIndex = idx;
              }
          }

          // --- Fit on ALL peaks ---
          for (int idx : peakIndices) {
              auto [peTemp, stdTemp, peErrTemp, stdErrTemp, reducedChi2Temp] = Analysis::fitGaussian(
                  gDerivOrig,
                  sMedie[idx] - 0.5, sMedie[idx] + 0.5,
                  derivate[idx],
                  sMedie[idx],
                  1.0
              );

              // Save results only of first spad
              if (idx == firstPeakIndex) {
                  pe = peTemp;
                  std = stdTemp;
                  peErr = peErrTemp;
                  stdErr = stdErrTemp;
                  reducedChi2 = reducedChi2Temp;

                  // --- DCR and CT calculation ONLY for the first spad ---
                  std::tie(DCR, DCR_Error) = Analysis::calculateCounts(thresholds, counts, 0.5 * pe);
                  auto [CT_value, CT_error] = Analysis::calculateCounts(thresholds, counts, 1.5 * pe);
                  std::tie(CT, CT_Error) = Analysis::calculateCTDCRRatio(DCR, DCR_Error, CT_value, CT_error);

                  ov = overvoltage;
                  infoTree->Fill();
              }
          }
      }



    // Set graph styles
    Int_t dynamicColor = ColorPaletteManager::getColorFromOvervoltage(overvoltage);
    GraphicsUtils::setGraphStyle(gRawOrig, 20, dynamicColor, 1);
    GraphicsUtils::setGraphStyle(gDerivOrig, 20, dynamicColor, 1);

    // Write graphs
    TString dirName = fname;
    dirName.ReplaceAll(".txt", "");
    TDirectory *subDir = outFile->mkdir(dirName.Data());
    subDir->cd();

    Utils::writeGraphTTree(gRawOrig, "gRaw");
    Utils::writeGraphTTree(gDerivOrig, "gDeriv");
    outFile->cd();

    // Store legend
    TString legendLabel;
    legendLabel.Form("@ %.1f OV", overvoltage);
    gRawOrig->SetTitle(legendLabel.Data());
    mgOrig->Add(gRawOrig, "LPE");

    LegendEntry entry = {overvoltage, legendLabel, dynamicColor, gRawOrig};
    legendEntries.push_back(entry);
  }

  // Create and style canvas
  auto cOrig = GraphicsUtils::createCanvas("cOrig", "All curves", 800, 600, true);
  GraphicsUtils::setAxisAndDraw(mgOrig, "Threshold (mV)", "Counts", true);

  // Sort legends
  std::sort(legendEntries.begin(), legendEntries.end(),
            [](const LegendEntry &a, const LegendEntry &b) {
              return a.overvoltage < b.overvoltage;
            });

  auto legOrig = GraphicsUtils::createLegend(0.7, 0.7, 0.89, 0.89);
  for (const auto &entry : legendEntries) {
    legOrig->AddEntry(entry.graph, entry.label.Data(), "pe");
  }
  legOrig->Draw();

  // Write everything
  cOrig->Write();
  outFile->Write();
  outFile->Close();
}
