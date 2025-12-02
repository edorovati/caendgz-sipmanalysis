#include <TFile.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include "../root-library/graphics.h"
#include "../root-library/gauss_stuff.h"

void analyze_fixed_threshold(const TString& filename, const TString& channel_name, bool doFit = true) {
   // --- Apri file ROOT ---
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore apertura file: " << filename << std::endl;
        return;
    }

    // ---- Accedi alla cartella fixed_threshold ----
    f->cd("constant_fraction_discrimination");
    TString hist_name = channel_name + "_cfd_dt_vs_amp";
    TH2F* h2d = (TH2F*)gDirectory->Get(hist_name);
    if (!h2d) {
        std::cerr << "Non trovo istogramma " << hist_name << std::endl;
        return;
    }
    double total_events = h2d->GetEntries();

    // ---- Accedi alla cartella calibration e recupera il TGraph ----
    f->cd("calibration");
       TString graph_name = channel_name + "_mean_vs_pe";
       TGraphErrors* gr_cal = (TGraphErrors*)gDirectory->Get(graph_name);
       if (!gr_cal) {
           std::cerr << "Non trovo TGraph " << graph_name << std::endl;
           f->Close();
           return;
       }

       double slope = 0;
       double intercept = 0;

       if (doFit) {
           // ---- Fit lineare al TGraph ----
           TF1* linearFit = new TF1("linearFit", "pol1", 1, 10);
           gr_cal->Fit(linearFit, "R"); // "Q" = quiet
           slope     = linearFit->GetParameter(1);
           intercept = linearFit->GetParameter(0);
       } else {


           // ---- Prendi solo il primo punto del TGraph ----
           if (gr_cal->GetN() > 1) {
               double x1, y1;
               gr_cal->GetPoint(1, x1, y1);  // il secondo punto ha indice 1
               slope = y1;
               intercept = 0;
           } else {
               std::cerr << "TGraph vuoto, non posso prendere il primo punto." << std::endl;
           }
       }


    // ---- Definizione intervalli p.e. come in probability ----
    std::vector<std::pair<double,double>> ranges;
    for (int i = 1; i < 6; ++i) {
        double start = (i - 0.5) * slope + intercept; // NOTA: centrato su i p.e.
        double end   = (i + 0.5) * slope + intercept;
        if (end > h2d->GetXaxis()->GetXmax()) break;
        if (start < h2d->GetXaxis()->GetXmin()) continue;
        ranges.push_back({start, end});
    }

    // ---- Disegna istogramma con linee verticali ----
    TCanvas* c1 = Graphics::CreateCanvas("c1", "Fixed Threshold Analysis", 1200, 800);
    Graphics::FormatAxis(h2d->GetXaxis(), h2d->GetYaxis(),
                         "#Deltat_{th-fixed} (ns)", "counts",
                         0.03, 0.03, 62, true, 0, 50, true, 30, 50);
    Graphics::Format2DHisto(h2d, kRainbow, "COLZ", 0.85, 0.04, 0.5, false);

    for (auto& r : ranges) {
        Graphics::CreateLine(r.first, 35, r.first, 45)->Draw("SAME");
        Graphics::CreateLine(r.second, 35, r.second, 45)->Draw("SAME");
    }

    // ---- Conta quante proiezioni sono "utili" ----
    int n_proj = 0;
    for (size_t i = 0; i < ranges.size(); ++i) {
        int bin_low  = h2d->GetXaxis()->FindBin(ranges[i].first);
        int bin_high = h2d->GetXaxis()->FindBin(ranges[i].second);
        TH1D* proj = h2d->ProjectionY(Form("proj_check_%zu", i), bin_low, bin_high);
        double frac = proj->Integral() / total_events * 100.0;
        delete proj;
        if (frac < 0.1) break;
        n_proj++;
    }

    // ---- Layout proporzionale ----
    int cols = TMath::CeilNint(TMath::Sqrt(n_proj));
    int rows = TMath::CeilNint((double)n_proj / cols);
    TCanvas* c2 = Graphics::CreateCanvas("c2", "Proiezioni", 1500, 800);
    c2->Divide(cols, rows);

    // ---- Loop proiezioni e fit ----
    const Style_t markerStyles[] = {45, 46, 47, 28, 34, 40, 41, 48};
    const int nMarkers = sizeof(markerStyles) / sizeof(markerStyles[0]);

    const Color_t colors[] = {kGreen+2, kBlue+2, kMagenta+2, kRed+2, kOrange+7, kViolet+2};
    const int nColors = sizeof(colors) / sizeof(colors[0]);

    std::vector<double> pe_index, sigma_ps, sigma_err_ps;
    int pad_index = 1;
    for (int i = 0; i < n_proj; ++i) {
        int bin_low  = h2d->GetXaxis()->FindBin(ranges[i].first);
        int bin_high = h2d->GetXaxis()->FindBin(ranges[i].second);
        TH1D* proj = h2d->ProjectionY(Form("fixed-th_%dpe", i+1), bin_low, bin_high);

        c2->cd(pad_index++);
        Graphics::SetupPad(false, true, true);

        Color_t col = colors[i % nColors];
        Style_t marker = markerStyles[i % nMarkers];

        Graphics::FormatAxis(proj->GetXaxis(), proj->GetYaxis(),
                             "#Deltat_{CFD} (ns)", "counts",
                             0.03, 0.03, 62,
                             true, 30, 50, false);

        //TF1* fitFunc = GaussStuff::FitQGaussPlusConst(proj, Form("fit_%dpe", i));
        TF1* fitFunc = GaussStuff::FitQGauss(proj, Form("fit_%dpe", i));
        Graphics::SetupHistogramStyle(proj, true, 111111, col, col, marker, 1, "HISTE", true, col, 3354);
        fitFunc->SetLineColor(col);
        fitFunc->Draw("SAME");

        double sigma_ns  = fitFunc->GetParameter(2);
        double sigma_err_ns = fitFunc->GetParError(2);
        pe_index.push_back(i+1);
        sigma_ps.push_back(sigma_ns * 1000.0);
        sigma_err_ps.push_back(sigma_err_ns * 1000.0);
    }

    // ---- Sigma vs indice p.e. ----
    TCanvas* c3 = Graphics::CreateCanvas("c3", "Sigma vs PE", 600, 600, false, false, true);
    // stampa i valori di sigma_ps / 1000
    std::cout << "sigma_ps (in ns): ";
    for (size_t i = 0; i < sigma_ps.size(); ++i) {
        double val_ns = sigma_ps[i] / 1000.0;
        double err_ns = sigma_err_ps[i] / 1000.0;
        std::cout << val_ns << " ± " << err_ns << " ";
    }
    std::cout << std::endl;

    // Creazione TGraphErrors (gli errori restano in ps se vuoi tenerli così,
    // altrimenti puoi fare la conversione qui direttamente)
    TGraphErrors* gr = new TGraphErrors(
        pe_index.size(),
        pe_index.data(),
        sigma_ps.data(),       // valori in ps, se vuoi in ns devi dividere per 1000
        nullptr,
        sigma_err_ps.data()    // errori in ps, se vuoi in ns devi dividere per 1000
    );


    Graphics::CustomizeGraph(gr, 47, kBlue+1, "APE", kBlue+2);
    Graphics::FormatAxis(gr->GetXaxis(), gr->GetYaxis(),
                         "#photo-electrons", "#sigma_{t} (ps)");
    gr->Draw("APE");

    c1->Update();
    c2->Update();
    c3->Update();
}
