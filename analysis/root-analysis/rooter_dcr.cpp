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


#include "root-library/utils.h"
#include "root-library/fit.h"

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
// ----------------------------------------------------------------------------
// @brief Structure to hold legend entry data
// ----------------------------------------------------------------------------
struct LegendEntry {
  double voltage;
  TString label;
  int color;
  TGraphErrors* graph;
};

// ----------------------------------------------------------------------------
// @brief Main function for DCR analysis
// ----------------------------------------------------------------------------
void rooter_dcr(const TString& outputFileName,
         const std::vector<std::tuple<std::string, std::string, std::string>>& input,
         const TString& sensorName = "",
         double sampling_rate = 750,       // in MHz
         int n_wf = 30000) {
  
  TMultiGraph* mgAll = new TMultiGraph(); // per contenere tutti i gRawOrig
  std::vector<LegendEntry> legendEntries;

  TFile* outFile = new TFile(outputFileName, "RECREATE");
  outFile->mkdir("Summary");
  outFile->cd("Summary");

  TTree* infoTree = new TTree("Info", "Information on overvoltages [current cycle]");
  double ov, pe, std, peErr, stdErr, DCR, CT;
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

  // Caricamento baseline reference
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

  int colorIndex = 1; // per colorare i grafici

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

    if (baselineAvailable && thresholds.size() == thresholdsRef.size()) {
      for (size_t i = 0; i < counts.size(); ++i) {
        counts[i] -= countsRef[i];
        yErrors[i] = std::round(sqrt(std::max(0.0, counts[i])));
      }
    } else {
      std::cerr << "⚠️  Baseline missing or incompatible: not applied to " << filePath << "\n";
    }
      
    TGraphErrors* gRawOrig = Utils::createGraph(thresholds, counts, xErrors, yErrors);

    auto [sMedie, derivate, derivErrors] = Utils::computeDerivative(thresholds, counts);
    TGraphErrors* gDerivOrig = Utils::createGraph(sMedie, derivate, std::vector<double>(sMedie.size(), 0.0), derivErrors);
      
    int peakIndex = Utils::findMaximum(sMedie, derivate, 4.5);
    Fitting::gauss_fit(gDerivOrig, peakIndex, pe, std, peErr, stdErr);
    
    double x1, y1, x2, y2;
    double CT_value = 0, CT_error = 0;
    int n = gRawOrig->GetN();
    if (pe > 2.0) {
      double DCR_counts = 0;
      // pe*0.5 → DCR
      for (int i = 0; i < n-1; ++i) {
        gRawOrig->GetPoint(i, x1, y1);
        gRawOrig->GetPoint(i+1, x2, y2);
        double thr = pe * 0.5;
        if (x1 <= thr && thr <= x2) {
          DCR_counts = y1 + (y2 - y1) * (thr - x1) / (x2 - x1);
          DCR_Error = std::sqrt(std::max(0.0, DCR_counts));
          DCR = DCR_counts / (n_wf * 1024 * (1.0 / sampling_rate));
          DCR_Error /= (n_wf * 1024 * (1.0 / sampling_rate));
          DCR *= 1000.0;
          DCR_Error *= 1000.0;
          break;
        }
      }

      // pe*1.5 → CT_value
      for (int i = 0; i < n-1; ++i) {
        gRawOrig->GetPoint(i, x1, y1);
        gRawOrig->GetPoint(i+1, x2, y2);
        double thr = pe * 1.5;
        if (x1 <= thr && thr <= x2) {
          CT_value = y1 + (y2 - y1) * (thr - x1) / (x2 - x1);
          CT_error = std::sqrt(std::max(0.0, CT_value));
          break;
        }
      }

      CT = (DCR_counts > 0) ? CT_value / DCR_counts : 0;
      CT_Error = (DCR_counts > 0) ? std::sqrt((CT * (1 - CT)) / DCR_counts) : 0;
    } else {
      CT = 0;
      CT_Error = 0;
      DCR = 0;
      DCR_Error = 0;
    }

    ov = voltage;
    infoTree->Fill();
    
    // Styling
    gRawOrig->SetLineColor(colorIndex);
    gRawOrig->SetMarkerColor(colorIndex);
    gRawOrig->SetMarkerStyle(20);
    
    mgAll->Add(gRawOrig, "LP");

    LegendEntry entry;
    entry.voltage = voltage;
    entry.label = Form("%.1f V", voltage);
    entry.color = colorIndex;
    entry.graph = gRawOrig;
    legendEntries.push_back(entry);
    
    colorIndex++;

    // Salvataggio dei grafici singoli
    TString dirName = fname;
    dirName.ReplaceAll(".txt", "");
    TDirectory* subDir = outFile->mkdir(dirName.Data());
    subDir->cd();
    gRawOrig->SetName("gRaw");
    gRawOrig->Write();
    gDerivOrig->SetName("gDerivative");
    gDerivOrig->Write();
    outFile->cd();
  }

  // Creazione canvas con tutti i grafici
  outFile->cd("Summary");
  TCanvas* cAll = new TCanvas("cAll_gRawOrig", "All gRawOrig vs Voltage", 1200, 800);
  mgAll->SetTitle("Raw counts vs Threshold;Threshold;Counts");
  mgAll->Draw("A");
  
  TLegend* leg = new TLegend(0.75, 0.7, 0.9, 0.9);
  for (const auto& entry : legendEntries) {
    leg->AddEntry(entry.graph, entry.label, "lp");
  }
  leg->Draw();

  cAll->Write(); // salva la canvas nel file
  
  outFile->Write();
  outFile->Close();
}
