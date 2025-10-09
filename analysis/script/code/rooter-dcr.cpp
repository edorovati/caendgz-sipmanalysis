/******************************************************************************
 * @file dcr_analysis.cpp
 * @brief Dark Count Rate (DCR) analysis macro.
 *
 * This script performs:
 *   - Reading of threshold‐count text files (data and optional reference).
 *   - Baseline subtraction using a reference file.
 *   - Numerical derivative to locate the single‐photon transition peak.
 *   - Gaussian fit of the derivative to extract 1 PE amplitude and σ.
 *   - Calculation of DCR and CT (cross‐talk) via threshold‐based counts.
 *   - Dynamic color assignment based on bias voltage (overvoltage).
 *   - Saving of graphs, fit parameters and DCR/CT results into a ROOT file.
 *
 * Usage:
 *   // Prepare a vector of (filepath, voltage-string, type) tuples:
 *   std::vector<std::tuple<std::string, std::string, std::string>> inputFiles = {
 *     {"data/sn41-A1/vbias_52.5.txt", "52.5", "data"},
 *     {"data/sn41-A1/vbias_53.0.txt", "53.0", "data"},
 *     // … add as many data files as needed …
 *     {"data/sn41-A1/vbias_52.0_ref.txt", "52.0", "reference"}  // optional baseline
 *   };
 *
 *   // In ROOT prompt or in a .C macro:
 *   .L dcr_analysis.cpp
 *   dcr("output.root", inputFiles, "MySensorName");
 *
 *   // Arguments:
 *   //   outputFileName : name of the ROOT file to create
 *   //   input          : vector of tuples (filePath, voltage, fileType)
 *   //                    fileType == "data" or "reference"
 *   //   sensorName     : optional label printed on all canvases
 *
 * Notes:
 *   - The reference file (first tuple with type "reference") is used
 *     to subtract background counts before fitting.
 *   - The 1 PE fit range is chosen around the derivative peak ± some window.
 *   - DCR is taken at 0.5·PE, CT at 1.5·PE, then CT/DCR ratio is computed.
 ******************************************************************************/


#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/analysis.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/graphics.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/utils.h"
#include "/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root/colorpalettemanager.h"

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
  double voltage;
  TString label;
  int color;
  TGraphErrors* graph;
};

std::vector<LegendEntry> legendEntries;

// ----------------------------------------------------------------------------
// @brief Main function for DCR analysis
// ----------------------------------------------------------------------------
void dcr(const TString& outputFileName,
         const std::vector<std::tuple<std::string, std::string, std::string>>& input,
         const TString& sensorName = "") {
  
  TMultiGraph* mgOrig = new TMultiGraph();
  std::vector<LegendEntry> legendEntries;
  TFile* outFile = new TFile(outputFileName, "RECREATE");
  outFile->mkdir("Summary");
  outFile->cd("Summary");

  TTree* infoTree = new TTree("Info", "Information on overvoltages [current cycle]");
  double ov, pe, std, peErr, stdErr, DCR, CT, reducedChi2;
  double DCR_Error, CT_Error;

  infoTree->Branch("voltage", &ov, "voltage/D");
  infoTree->Branch("pe", &pe, "pe/D");
  infoTree->Branch("stdDev", &std, "stdDev/D");
  infoTree->Branch("peError", &peErr, "peError/D");
  infoTree->Branch("stdDevError", &stdErr, "stdDevError/D");
  infoTree->Branch("DCR", &DCR, "DCR/D");
  infoTree->Branch("CT", &CT, "CT/D");
  infoTree->Branch("DCR_Error", &DCR_Error, "DCR_Error/D");
  infoTree->Branch("CT_Error", &CT_Error, "CT_Error/D");

  std::vector<double> thresholdsRef, countsRef, errorsRef;
  bool baselineAvailable = false;

  for (const auto& [filePathStr, voltageStr, fileType] : input) {
    if (fileType == "reference") {
        TString filePath = filePathStr;
        if (!Utils::readDataFile(filePath, thresholdsRef, countsRef)) {
          std::cerr << "WARNING: Reference file " << filePath << " could not be read.\n";
        } else {
          baselineAvailable = true;
          errorsRef.resize(countsRef.size());
          for (size_t i = 0; i < countsRef.size(); ++i)
            errorsRef[i] = std::round(sqrt(countsRef[i]));
          std::cout << "✅ Reference baseline loaded from " << filePath << "\n";
        }
        break;
    }
  }


  // Process data files
  for (const auto& [filePathStr, voltageStr, fileType] : input) {
    if (fileType != "data") continue;

    TString filePath = filePathStr;
    TString fname = gSystem->BaseName(filePath);
    double voltage = atof(voltageStr.c_str());

    std::cout << "➡️  Processing file: " << filePath << " @ " << voltage << "V\n";

    std::vector<double> thresholds, counts;
    if (!Utils::readDataFile(filePath, thresholds, counts) || thresholds.empty()) {
      std::cerr << "❌ Error reading file " << filePath << "\n";
      continue;
    }

    std::vector<double> xErrors(thresholds.size(), 0.0);
    std::vector<double> yErrors(counts.size());
    for (size_t i = 0; i < counts.size(); ++i)
        yErrors[i] = std::round(sqrt(std::max(0.0, counts[i])));


    TGraphErrors* gRawOrig = Utils::createGraph(thresholds, counts, xErrors, yErrors);

    auto [sMedie, derivate] = Utils::computeDerivative(thresholds, counts);
    std::vector<double> derivErrors;
    size_t N = std::min(counts.size(), thresholds.size());
    for (size_t i = 0; i < N; ++i) {
      double dx = thresholds[i+1] - thresholds[i];
      double err = std::round( std::sqrt( std::abs(derivate[i]) ) / std::abs(dx) );
      derivErrors.push_back(err);
    }

    TGraphErrors* gDerivOrig = Utils::createGraph(sMedie, derivate, std::vector<double>(sMedie.size(), 0.0), derivErrors);

    int peakIndex = Utils::findMaximum(sMedie, derivate, 4.5);
    if (peakIndex != -1) {
      std::tie(pe, std, peErr, stdErr, reducedChi2) = Analysis::fitGaussian(
        gDerivOrig, sMedie[peakIndex] - 1.5, sMedie[peakIndex] + 2.0,
        derivate[peakIndex], sMedie[peakIndex], 1.0);

      if (baselineAvailable && thresholds.size() == thresholdsRef.size()) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] -= countsRef[i];
            yErrors[i] = std::round(sqrt(std::max(0.0, counts[i]))); 
        }
      } else {
          std::cerr << "⚠️  Baseline missing or incompatible: not applied to " << filePath << "\n";
          for (size_t i = 0; i < counts.size(); ++i)
          yErrors[i] = std::round(sqrt(std::max(0.0, counts[i])));
      }


      std::tie(DCR, DCR_Error) = Analysis::calculateCounts(thresholds, counts, 0.5 * pe);
      auto [CT_value, CT_error] = Analysis::calculateCounts(thresholds, counts, 1.5 * pe);
      std::tie(CT, CT_Error) = Analysis::calculateCTDCRRatio(DCR, DCR_Error, CT_value, CT_error);

      ov = voltage;
      infoTree->Fill();
    } else {
      std::cerr << "⚠️  No peak found for file " << filePath << "\n";
    }

    // Styling and graph storage
    Int_t dynamicColor = ColorPaletteManager::getColorFromOvervoltage(voltage);
    GraphicsUtils::setGraphStyle(gRawOrig, 20, dynamicColor, 1);
    GraphicsUtils::setGraphStyle(gDerivOrig, 20, dynamicColor, 1);

    TString dirName = fname;
    dirName.ReplaceAll(".txt", "");
    TDirectory* subDir = outFile->mkdir(dirName.Data());
    subDir->cd();

    Utils::writeGraphTTree(gRawOrig, "gRaw");
    Utils::writeGraphTTree(gDerivOrig, "gDeriv");

    // Create individual canvases with labels
    auto cRaw = GraphicsUtils::createCanvas("cRaw_" + dirName, "Raw counts", 800, 600, true);
    TMultiGraph* mgSingle = new TMultiGraph();
    mgSingle->Add(gRawOrig, "LPE");
    GraphicsUtils::setAxisAndDraw(mgSingle, "Threshold (mV)", "Counts", true);
    if (!sensorName.IsNull() && !sensorName.IsWhitespace()) {
      TString label = TString::Format("Sensor: %s  |  V = %.1f V", sensorName.Data(), voltage);
      GraphicsUtils::addLabelToCanvas(cRaw, label, 0.15, 0.82, 0.55, 0.88, 0, 12, 1, 42);
    }
    cRaw->Write();

    auto cDeriv = GraphicsUtils::createCanvas("cDeriv_" + dirName, "Derivative", 800, 600, true);
    TMultiGraph* mgDeriv = new TMultiGraph();
    mgDeriv->Add(gDerivOrig, "LPE");
    GraphicsUtils::setAxisAndDraw(mgDeriv, "Threshold (mV)", "First Derivative", true);
    if (!sensorName.IsNull() && !sensorName.IsWhitespace()) {
      TString label = TString::Format("Sensor: %s  |  @%.1f V", sensorName.Data(), voltage);
      GraphicsUtils::addLabelToCanvas(cDeriv, label, 0.15, 0.82, 0.55, 0.88, 0, 12, 1, 42);
    }
    TString label2 = TString::Format("p.e.: %.2f #pm %.2f (mV) |  #sigma: %.2f #pm %.2f (mV) |  #chi^{2}/ndf: %.2f", pe, peErr, std, stdErr, reducedChi2);
    GraphicsUtils::addLabelToCanvas(cDeriv, label2, 0.15, 0.75, 0.55, 0.81, 0, 12, 1, 42);

    cDeriv->Write();

    outFile->cd();

    TString legendLabel;
    legendLabel.Form("@ %.1f V", voltage);
    gRawOrig->SetTitle(legendLabel.Data());
    mgOrig->Add(gRawOrig, "LPE");
    legendEntries.push_back({voltage, legendLabel, dynamicColor, gRawOrig});
  }

  // Final combined canvas
  auto cOrig = GraphicsUtils::createCanvas("cOrig", "All curves", 800, 600, true);
  GraphicsUtils::setAxisAndDraw(mgOrig, "Threshold (mV)", "Counts", true);

  std::sort(legendEntries.begin(), legendEntries.end(),
            [](const LegendEntry& a, const LegendEntry& b) {
              return a.voltage < b.voltage;
            });

  auto legOrig = GraphicsUtils::createLegend(0.7, 0.7, 0.89, 0.89);
  for (const auto& entry : legendEntries)
    legOrig->AddEntry(entry.graph, entry.label.Data(), "pe");
  legOrig->Draw();

  if (!sensorName.IsNull() && !sensorName.IsWhitespace()) {
    TString label = "Sensor: " + sensorName;
    GraphicsUtils::addLabelToCanvas(cOrig, label, 0.15, 0.82, 0.45, 0.88, 0, 12, 1, 42);
  }

  cOrig->Write();
  outFile->Write();
  outFile->Close();
}
