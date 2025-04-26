/******************************************************************************
 * @file dcr_plot.cpp
 * @brief Script for analyzing and plotting dark count rate (DCR), cross-talk (CT),
 *        and single photoelectron (p.e.) signal mean from ROOT data.
 *
 * This script performs the following operations:
 *   - Opens the previously generated ROOT file ("output.root").
 *   - Reads the "Info" TTree from the "Summary" directory.
 *   - Extracts overvoltage, mean signal, DCR, and CT data with their errors.
 *   - Calculates the dark count rate normalized to the acquisition window.
 *   - Constructs TGraphErrors for the mean signal, DCR, and CT.
 *   - Performs a linear fit on the mean signal vs overvoltage plot.
 *   - Draws graphs with error bands and fit results.
 *   - Generates additional plots for DCR and CT versus overvoltage.
 *
 ******************************************************************************/

#include "../../root/analysis.h"
#include "../../root/graphics.h"
#include "../../root/utils.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TTree.h"
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

// Define static constants
const double Analysis::wf = 80e3;
const double Analysis::timeWindow = 1.36e-6;

void dcr_plot() {
  // Open input ROOT file
  std::unique_ptr<TFile> inputFile(TFile::Open("output.root", "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Errore nell'apertura di output.root" << std::endl;
    return;
  }

  // Access "Summary" directory and retrieve TTree
  inputFile->cd("Summary");
  TTree *tree = static_cast<TTree *>(gDirectory->Get("Info"));
  if (!tree) {
    std::cerr << "TTree 'Info' non trovato nella directory Summary!" << std::endl;
    return;
  }

  // Define variables for reading TTree branches
  double ov, mean, meanErr, stdDev, DCR, DCR_Err, CT, CT_Err;
  tree->SetBranchAddress("overvoltage", &ov);
  tree->SetBranchAddress("pe", &mean);
  tree->SetBranchAddress("peError", &meanErr);
  tree->SetBranchAddress("stdDev", &stdDev);
  tree->SetBranchAddress("DCR", &DCR);
  tree->SetBranchAddress("DCR_Error", &DCR_Err);
  tree->SetBranchAddress("CT", &CT);
  tree->SetBranchAddress("CT_Error", &CT_Err);

  const int nEntries = tree->GetEntries();

  // Prepare vectors to hold data for plotting
  std::vector<double> ovV, meanV, meanErrV;
  std::vector<double> dcrV, dcrErrV;
  std::vector<double> ctV, ctErrV;
  ovV.reserve(nEntries);
  meanV.reserve(nEntries);
  meanErrV.reserve(nEntries);
  dcrV.reserve(nEntries);
  dcrErrV.reserve(nEntries);
  ctV.reserve(nEntries);
  ctErrV.reserve(nEntries);

  std::vector<BandPoint> bandData;
  bandData.reserve(nEntries);

  // Loop over tree entries to extract and compute data
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    auto darkPair = Analysis::calculateDarkRateAndError(
        DCR, DCR_Err, Analysis::wf, Analysis::timeWindow);
    double darkRate = darkPair.first;
    double darkRateErr = darkPair.second;

    ovV.push_back(ov);
    meanV.push_back(mean);
    meanErrV.push_back(meanErr);
    dcrV.push_back(darkRate);
    dcrErrV.push_back(darkRateErr);
    ctV.push_back(CT);
    ctErrV.push_back(CT_Err);

    bandData.push_back({ov, mean + stdDev, mean - stdDev});
  }

  // Create graphs
  std::vector<double> zeroErr(nEntries, 0.0);
  TGraphErrors *gMean = Utils::createGraph(ovV, meanV, zeroErr, meanErrV);
  TGraphErrors *gDCR = Utils::createGraph(ovV, dcrV, zeroErr, dcrErrV);
  TGraphErrors *gCT = Utils::createGraph(ovV, ctV, zeroErr, ctErrV);

  GraphicsUtils::setGraphStyle(gMean, 20, kBlue);
  GraphicsUtils::setGraphStyle(gDCR, 20, kRed);
  GraphicsUtils::setGraphStyle(gCT, 20, kGreen + 2);

  // Create and sort band graph
  std::sort(bandData.begin(), bandData.end(), [](auto &a, auto &b) {
    return a.x < b.x;
  });
  TGraph *gBand = GraphicsUtils::createBandGraph(bandData);

  // Draw main graph with band
  GraphicsUtils::drawGraphWithBand(gMean, bandData,
                                   "Overvoltage vs Mean + Band",
                                   "Overvoltage (V)", "Single p.e. (mV)", "c1");

  // Perform linear fit
  double minX = *std::min_element(ovV.begin(), ovV.end());
  double maxX = *std::max_element(ovV.begin(), ovV.end());
  TF1 *fit = new TF1("fit", "pol1", minX, maxX);
  gMean->Fit(fit, "RQ");

  double chi2 = fit->GetChisquare();
  int ndf = fit->GetNDF();
  double chi2red = (ndf != 0) ? chi2 / ndf : 0.0;

  // Prepare and draw fit information
  std::vector<std::string> fitResults = {
      "Fit: y = a*x + b",
      Form("a = %.1f #pm %.1f [mV/ov]", fit->GetParameter(1), fit->GetParError(1)),
      Form("b = %.1f #pm %.1f [mV]", fit->GetParameter(0), fit->GetParError(0)),
      Form("#chi^{2}/ndf = %.2f", chi2red)
  };
  TPaveText *pave = GraphicsUtils::createPaveText(0.1, 0.7, 0.35, 0.9, fitResults);
  pave->Draw();

  // Create and draw legend
  TLegend *leg = GraphicsUtils::createLegend(0.6, 0.7, 0.9, 0.9);
  GraphicsUtils::addEntryToLegend(leg, gMean, "Single p.e.", "lep");
  GraphicsUtils::addEntryToLegend(leg, gBand, "Single p.e. #pm #sigma", "f");
  leg->Draw();

  // Helper lambda for additional plots
  auto makeCanvasPlot = [&](TGraphErrors *g, const char *name,
                            const char *title, const char *yTitle) {
    TCanvas *c = GraphicsUtils::createCanvas(name, title, 800, 600);
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(g, "PE");
    mg2->SetTitle(Form("%s;Overvoltage (V);%s", title, yTitle));
    GraphicsUtils::setAxisAndDraw(mg2, "Overvoltage (V)", yTitle);
  };

  // Create DCR and CT plots
  makeCanvasPlot(gDCR, "c2", "Overvoltage vs DCR", "DCR (kHz)");
  makeCanvasPlot(gCT, "c3", "Overvoltage vs CT", "CT probability (%)");

  inputFile->Close();
  std::cout << "Grafici creati con successo!" << std::endl;
}
