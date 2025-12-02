#include <TFile.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include "../root-library/graphics.h"
#include "../root-library/utils.h"
#include "../root-library/analysis.h"

void crossing(const TString& filename, const TString& channel_name, bool doFit = true) {
    // --- Apri file ROOT ---
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore apertura file: " << filename << std::endl;
        return;
    }

    // ---- Accedi alla cartella fixed_threshold ----
    f->cd("fixed_threshold");

    TString hist_name_first  = channel_name + "_fixed_th_diff_laser";        // first crossing
    TString hist_name_signal = channel_name + "_fixed_th_diff_laser_clean"; // signal crossing

    TH2F* h2d_first  = (TH2F*)gDirectory->Get(hist_name_first);
    TH2F* h2d_signal = (TH2F*)gDirectory->Get(hist_name_signal);

    if (!h2d_first || !h2d_signal) {
        std::cerr << "Non trovo istogrammi richiesti!" << std::endl;
        return;
    }

    // ---- Accedi alla cartella calibration ----
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
        // Fit lineare
        TF1* linearFit = new TF1("linearFit", "pol1", 1, 10);
        gr_cal->Fit(linearFit, "RQ"); // "Q" = quiet
        slope     = linearFit->GetParameter(1);
        intercept = linearFit->GetParameter(0);
        std::cout << slope << " " << intercept << "\n";
    } else {
        if (gr_cal->GetN() > 1) {
            double x0, y0;
            gr_cal->GetPoint(1, x0, y0);
            slope = y0;
            intercept = 0;
        } else {
            std::cerr << "TGraph vuoto, non posso prendere il primo punto." << std::endl;
        }
    }

    // ---- Recupera n_waveforms ----
    f->cd();
    TTree* infoTree = (TTree*)f->Get("info");
    Int_t n_wf = 0;
    if (infoTree) {
        infoTree->SetBranchAddress("n_waveforms", &n_wf);
        infoTree->GetEntry(0);
        std::cout << "[INFO] Numero di waveforms: " << n_wf << std::endl;  // <<< stampa qui

    } else {
        std::cerr << "Attenzione: infoTTree non trovato, n_waveforms = 0" << std::endl;
    }

    // ---- Soglia a 0.5 p.e. ----
    double low_cut  = 0.5*slope + intercept;

    // Prendi da low_cut fino all’ultimo bin disponibile
    int bin_low_first  = h2d_first->GetXaxis()->FindBin(low_cut);
    int bin_high_first = h2d_first->GetNbinsX();

    int bin_low_signal  = h2d_signal->GetXaxis()->FindBin(low_cut);
    int bin_high_signal = h2d_signal->GetNbinsX();

    // ---- Proiezioni ----
    TH1D* proj_first  = h2d_first->ProjectionY("proj_first",  bin_low_first,  bin_high_first);
    TH1D* proj_signal = h2d_signal->ProjectionY("proj_signal", bin_low_signal, bin_high_signal);

    // ---- Trova bin massimo e calcola S e B ----
    int b_first  = Utils::FindMaxBin(proj_first);
    int b_signal = Utils::FindMaxBin(proj_signal);

    Analysis::SignalBackground sb_first  = Analysis::CalculateSignalAndBackground(proj_first,  b_first);
    Analysis::SignalBackground sb_signal = Analysis::CalculateSignalAndBackground(proj_signal, b_signal);

    // ---- Probabilità di avere segnale ----
    auto [prob_first,  prob_first_err]  = Analysis::CalculateProbability(sb_first.diff,  n_wf);
    auto [prob_signal, prob_signal_err] = Analysis::CalculateProbability(sb_signal.diff, n_wf);

    std::cout << "=== Probabilità di segnale ===" << std::endl;
    std::cout << "First crossing : " << prob_first*100.0  << " ± " << prob_first_err*100.0  << " %" << std::endl;
    std::cout << "Signal crossing: " << prob_signal*100.0 << " ± " << prob_signal_err*100.0 << " %" << std::endl;
    std::cout << "=============================" << std::endl;

    // ---- Disegna confronto ----
    TCanvas* c_compare = Graphics::CreateCanvas("c_compare", "First vs Signal Crossing", 1000, 700);

    Graphics::FormatAxis(proj_first->GetXaxis(), proj_first->GetYaxis(),
                         "#Deltat (ns)", "counts", 0.03, 0.03, 62, true, 30, 50);

    Graphics::SetupHistogramStyle(proj_first, false, 1111,
                                  kBlue+2, kBlue+2, 47, 2, "E1", false);

    Graphics::SetupHistogramStyle(proj_signal, false, 1111,
                                  kRed+1, kRed+1, 46, 2, "E1 SAME", false);

    TLegend* leg = Graphics::CreateLegend(0.15, 0.85, 0.4, 0.7);
    Graphics::AddLegendEntry(leg, proj_first,  "first crossing",  "lpe");
    Graphics::AddLegendEntry(leg, proj_signal, "signal crossing", "lpe");
    leg->Draw();

    TString infoText;
    infoText += Form("P_{first}  = %.2f #pm %.2f %%\n",  prob_first*100.0,  prob_first_err*100.0);
    infoText += Form("P_{signal} = %.2f #pm %.2f %%\n", prob_signal*100.0, prob_signal_err*100.0);

    TPaveText* infoBox = Graphics::CreateInfoBox(0.55, 0.75, 0.9, 0.9,
                                                 infoText, kWhite, kBlack, 0.04, 1);
    infoBox->Draw();

    c_compare->Update();
}
