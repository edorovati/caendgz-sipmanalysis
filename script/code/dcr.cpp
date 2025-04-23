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
#include "../../root/graphics.h"
#include "../../root/utils.h"
#include "../../root/analysis.h"

struct LegendEntry {
  double overvoltage;
  TString label;
  int color;
  TGraphErrors *graph;
};
std::vector<LegendEntry> legendEntries;

void dcr() {
  TFile *outFile = new TFile("output.root", "RECREATE");

  TMultiGraph *mgOrig = new TMultiGraph();

  std::vector<double> allOvervoltages;

  outFile->mkdir("Summary");
  outFile->cd("Summary");
  TTree *infoTree =
      new TTree("Info", "Informazioni sugli overvoltages [current cycle]");
  double ov, pe, std, peErr, stdErr, DCR, CT, reducedChi2;
  double DCR_Error, CT_Error;

  infoTree->Branch("overvoltage", &ov, "overvoltage/D");
  infoTree->Branch("pe", &pe, "pe/D");
  infoTree->Branch("stdDev", &std, "stdDev/D");
  infoTree->Branch("peError", &peErr, "peError/D");
  infoTree->Branch("stdDevError", &stdErr, "stdDevError/D");
  infoTree->Branch("DCR", &DCR, "DCR/D");
  infoTree->Branch("CT", &CT, "CT/D");
  infoTree->Branch("DCR_Error", &DCR_Error, "DCR_Error/D");
  infoTree->Branch("CT_Error", &CT_Error, "CT_Error/D");
  std::vector<int> customColors = {kBlue,   kRed,  kGreen + 2, kMagenta,
                                   kOrange, kCyan, kViolet,    kGray};
  int colorIndex = 0;

  TSystemDirectory dir("Scan", "Scan");
  TList *files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "Nessun file trovato nella cartella Scan!" << std::endl;
    return;
  }

  TSystemFile *fileItem = nullptr;
  TIter next(files);
  while ((fileItem = (TSystemFile *)next())) {
    TString fname = fileItem->GetName();
    if (fileItem->IsDirectory() || fname.BeginsWith(".") || fname == "0OV.txt")
      continue;
    if (!fname.Contains("OV") || !fname.EndsWith(".txt"))
      continue;

    std::vector<double> soglie;
    std::vector<double> conteggi;

    if (!Utils::readDataFile(TString("Scan/") + fname, soglie, conteggi))
      continue;
    int N = soglie.size();
    if (N == 0)
      continue;

    double overvoltage = Utils::extractOvervoltage(fname);
    if (overvoltage != 0.0) {
      allOvervoltages.push_back(overvoltage);
    }

    std::vector<double> yErrors;
    for (size_t i = 0; i < conteggi.size(); ++i) {
      yErrors.push_back(sqrt(conteggi[i]));
    }

    std::vector<double> xErrors(soglie.size(), 0.25);
    TGraphErrors *gRawOrig =
        Utils::createGraph(soglie, conteggi, xErrors, yErrors);

    auto [sMedie, derivate, derivateErr] =
        Utils::computeDerivativeWithErrors(soglie, conteggi, yErrors, false);
    std::vector<double> xErrors_deriv(sMedie.size(), 0.0);
    TGraphErrors *gDerivOrig =
        Utils::createGraph(sMedie, derivate, xErrors_deriv, derivateErr);

    int peakIndex = Utils::findMaximum(sMedie, derivate, 6.0);

    if (peakIndex != -1) {
      std::tie(pe, std, peErr, stdErr, reducedChi2) = Analysis::fitGaussian(
          gDerivOrig, sMedie[peakIndex] - 4.0, sMedie[peakIndex] + 4.0,
          derivate[peakIndex], sMedie[peakIndex], 1.0);
      std::tie(DCR, DCR_Error) =
          Analysis::calculateCounts(soglie, conteggi, 0.5 * pe);
      auto [CT_value, CT_error] =
          Analysis::calculateCounts(soglie, conteggi, 1.5 * pe);
      std::tie(CT, CT_Error) =
          Analysis::calculateCTDCRRatio(DCR, DCR_Error, CT_value, CT_error);
      ov = overvoltage;
      infoTree->Fill();
    }

    GraphicsUtils::setGraphStyle(
        gRawOrig, 20, customColors[colorIndex % customColors.size()], 1);

    ++colorIndex;

    GraphicsUtils::setGraphStyle(gDerivOrig, 20, gRawOrig->GetLineColor(), 1);

    TString dirName = fname;
    dirName.ReplaceAll(".txt", "");
    TDirectory *subDir = outFile->mkdir(dirName.Data());
    subDir->cd();

    Utils::writeGraphTTree(gRawOrig, "gRaw");
    Utils::writeGraphTTree(gDerivOrig, "gDeriv");
    outFile->cd();

    TString legendLabel;
    legendLabel.Form("@ %.1f OV", overvoltage);
    gRawOrig->SetTitle(legendLabel.Data());
    mgOrig->Add(gRawOrig, "LPE");
    LegendEntry entry;
    entry.overvoltage = overvoltage;
    entry.label = legendLabel;
    entry.color = gRawOrig->GetLineColor();
    entry.graph = gRawOrig;
    legendEntries.push_back(entry);
  }

  auto cOrig = GraphicsUtils::createCanvas("cOrig", "Tutte le curve originali",
                                           800, 600, true);
  GraphicsUtils::setAxisAndDraw(mgOrig, "Threshold (mV)", "Counts", true);
  std::sort(legendEntries.begin(), legendEntries.end(),
            [](const LegendEntry &a, const LegendEntry &b) {
              return a.overvoltage < b.overvoltage;
            });

  auto legOrig = GraphicsUtils::createLegend(0.7, 0.7, 0.89, 0.89);
  for (const auto &entry : legendEntries) {
    legOrig->AddEntry(entry.graph, entry.label.Data(), "pe");
  }
  legOrig->Draw();
  cOrig->Write();
  outFile->Write();
  outFile->Close();
}
