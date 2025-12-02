#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "../root-library/graphics.h"

// ========================
// Funzione di debug: stampa il contenuto di una directory ricorsivamente
// ========================
void PrintDirectory(TDirectory* dir, int level = 0) {
    if (!dir) return;
    TString indent(level * 2, ' ');
    std::cout << indent << "Directory: " << dir->GetPath() << std::endl;

    TIter nextkey(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        std::cout << indent << "  " << key->GetName()
                  << " (" << key->GetClassName() << ")" << std::endl;
        if (strcmp(key->GetClassName(), "TDirectoryFile") == 0) {
            TDirectory* subdir = (TDirectory*)dir->Get(key->GetName());
            PrintDirectory(subdir, level + 1);
        }
    }
}

// ========================
// Helper: Fit salita tra 80%-20% del massimo
// ========================
TF1* FitRisingEdge(TProfile* prof, double maxPos, double maxVal) {
    double thr80 = 0.8 * maxVal;
    double thr20 = 0.2 * maxVal;

    int bin80 = -1, bin20 = -1;

    // Cerca il bin di 80% partendo dall'inizio (o vicino al massimo)
    for (int b = 1; b <= prof->GetNbinsX(); ++b) {
        double y = prof->GetBinContent(b);
        if (y >= thr80) {
            bin80 = b;
            break;
        }
    }

    if (bin80 < 0) {
        std::cerr << "[DEBUG] Non trovo bin80" << std::endl;
        return nullptr;
    }

    // Cerca bin20 andando all'indietro partendo da bin80
    for (int b = bin80; b >= 1; --b) {
        double y = prof->GetBinContent(b);
        if (y <= thr20) {
            bin20 = b;
            break;
        }
    }

    if (bin20 < 0) {
        std::cerr << "[DEBUG] Non trovo bin20" << std::endl;
        return nullptr;
    }

    double x1 = prof->GetBinCenter(bin20);
    double x2 = prof->GetBinCenter(bin80);

    std::cout << "[DEBUG] Max: " << maxVal << " at t=" << maxPos
              << ", t20=" << x1 << ", t80=" << x2 << std::endl;

    TF1* linFit = new TF1("linFit", "pol1", x1, x2);
    prof->Fit(linFit, "IMREQ0");
    return linFit;
}


// ========================
// Helper: Fit discesa con exp dal massimo fino al 5%
// ========================
TF1* FitFallingEdge(TProfile* prof, double maxPos, double maxVal) {
    double thr5 = 0.05 * maxVal;
    int bin5 = -1;
    for (int b = prof->FindBin(maxPos); b <= prof->GetNbinsX(); ++b) {
        if (prof->GetBinContent(b) <= thr5) { bin5 = b; break; }
    }
    if (bin5 < 0) return nullptr;

    double x1 = maxPos;
    double x2 = prof->GetBinCenter(bin5);

    TF1* expFit = new TF1("expFit", "expo", x1, x2);
    prof->Fit(expFit, "IMREQ0");
    return expFit;
}

// ========================
// Funzione principale di analisi
// ========================
void analyze_full(const TString& filename, const TString& channel_name, bool doFit = true) {
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERRORE] Non riesco ad aprire il file: " << filename << std::endl;
        return;
    }

    std::cout << "===== Contenuto del file ROOT =====" << std::endl;
    PrintDirectory(f);
    // ====================
    // --- PROFILES + CLASSIFICATION ---
    // ====================
    std::map<int,double> baseline_sigma;
    TCanvas* c_prof = Graphics::CreateCanvas("c_prof","Profiles overlay",1200,800);
    c_prof->Divide(3,2);
    int ipad=1;

    for (int npe=1;npe<=5;++npe) {
        // Profili
        TProfile* prof = (TProfile*)f->Get("profile_classification_first_crossing/" + channel_name + Form("_class_%d_profile",npe));
        if (!prof) {
            std::cerr << "[ERRORE] Non trovo profilo npe=" << npe << std::endl;
            continue;
        }
        
        int binmax = prof->GetMaximumBin();
        double maxVal = prof->GetBinContent(binmax);
        double maxPos = prof->GetBinCenter(binmax);

        TF1* fRise = FitRisingEdge(prof,maxPos,maxVal);
        TF1* fFall = FitFallingEdge(prof,maxPos,maxVal);

        // Istogramma 2D
        TH2F* h2dclass = (TH2F*)f->Get("classification_first_crossing/" + channel_name + Form("_class_%d",npe));
        if (!h2dclass) std::cerr << "[ERRORE] Non trovo TH2F class npe=" << npe << std::endl;

        c_prof->cd(ipad++);
        if (h2dclass) h2dclass->Draw("COLZ");
        Graphics::FormatAxis(h2dclass->GetXaxis(), h2dclass->GetYaxis(),
                             "ref time - hit time (ns)", "amplitude (mV)",
                             0.03, 0.03, 62, false, -40, 250, false, 30, 50);
        Graphics::Format2DHisto(h2dclass, kCividis, "COLZ", 0.85, 0.04, 0.5);

        prof->SetLineColor(kBlack);
        prof->SetLineStyle(6);
        prof->Draw("SAME");
        if (fRise) fRise->Draw("SAME");
        if (fFall) fFall->Draw("SAME");
        Graphics::SetupPad(false, false, false, true);
    }

    // ====================
    // --- BASELINE ---
    // ====================
    TCanvas* c_base = Graphics::CreateCanvas("c_base","Baseline fits",1200,800);
    c_base->Divide(3,2);
    ipad=1;
    const Style_t markerStyles[] = {45, 46, 47, 28, 34, 40, 41, 48};
    const int nMarkers = sizeof(markerStyles) / sizeof(markerStyles[0]);

    const Color_t colors[] = {kGreen+2, kBlue+2, kMagenta+2, kRed+2, kOrange+7, kViolet+2};
    const int nColors = sizeof(colors) / sizeof(colors[0]);

    for (int npe=1;npe<=5;++npe) {
        TH1D* hbase = (TH1D*)f->Get("baseline/" + channel_name + Form("_class_%d_amp_baseline",npe));
        if (!hbase) {
            std::cerr << "[ERRORE] Non trovo baseline npe=" << npe << std::endl;
            continue;
        }

        c_base->cd(ipad++);
        Graphics::SetupPad(false, false, true);
        Color_t col = colors[npe % nColors];
        Style_t marker = markerStyles[npe % nMarkers];
        Graphics::FormatAxis(hbase->GetXaxis(), hbase->GetYaxis(),
                             "amplitude (mV)", "counts",
                             0.03, 0.03, 62,
                             true, -10, 10, false);
        Graphics::SetupHistogramStyle(hbase, true, 111111, col, col, marker, 1, "HISTE", true, col, 3354, 0.35, 0.6, 0.9, 0.45, 0.9);
        TF1* gfit = new TF1("gfit","gaus", hbase->GetXaxis()->GetXmin(), hbase->GetXaxis()->GetXmax());
        gfit->SetNpx(1000);
        hbase->Fit(gfit,"FQ");
        gfit->SetLineColor(col);
        gfit->Draw("SAME");
        baseline_sigma[npe] = gfit->GetParameter(2);
    }
    
    
    // ====================
    // --- FIXED THRESHOLD con calibrazione ---
    // ====================

    // --- Step 1: recupera calibrazione ---
    f->cd("calibration");
    TString graph_name = channel_name + "_mean_vs_pe";
    TGraphErrors* gr_cal = (TGraphErrors*)gDirectory->Get(graph_name);
    if (!gr_cal) {
        std::cerr << "[ERRORE] Non trovo TGraph " << graph_name << std::endl;
        return;
    }

    double slope = 0, intercept = 0;
    if (doFit) {
        // Fit lineare del TGraph
        TF1* linearFit = new TF1("linearFit", "pol1", 1, 10);
        gr_cal->Fit(linearFit, "RQ"); // R = range fit, Q = quiet
        slope     = linearFit->GetParameter(1);
        intercept = linearFit->GetParameter(0);
        std::cout << "[DEBUG] Calibrazione con fit: slope=" << slope
                  << " intercept=" << intercept << std::endl;
    } else {
        // Prendi la differenza tra i primi due punti come slope
        if (gr_cal->GetN() > 1) {
            double x1,y1, x2,y2;
            gr_cal->GetPoint(0,x1,y1);
            gr_cal->GetPoint(1,x2,y2);
            slope = y2 - y1;
            intercept = y1 - slope*x1;
            std::cout << "[DEBUG] Calibrazione da punti: slope=" << slope
                      << " intercept=" << intercept << std::endl;
        }
    }

    // --- Step 2: usa calibrazione per definire i ranges ---
    TH2F* h2d_fixed = (TH2F*)f->Get("fixed_threshold/" + channel_name + "_fixed_th_diff_laser");
    if (!h2d_fixed) {
        std::cerr << "[ERRORE] Non trovo TH2F fixed_threshold per " << channel_name << std::endl;
        return;
    }

    std::vector<std::pair<double,double>> ranges;
    for (int i = 1; i <= 5; ++i) {
        double start = (i - 0.5)*slope + intercept;
        double end   = (i + 0.5)*slope + intercept;
        if (end > h2d_fixed->GetXaxis()->GetXmax()) break;
        if (start < h2d_fixed->GetXaxis()->GetXmin()) continue;
        ranges.push_back({start, end});
        std::cout << "[DEBUG] Range p.e.=" << i << " : " << start << " -> " << end << std::endl;
    }

    // --- Step 3: proietta e fitta ---
    std::vector<double> pe_index, sigma_ps, sigma_err_ps;
    TCanvas* c_fix = Graphics::CreateCanvas("c_fix","Fixed threshold projections",1200,800);
    c_fix->Divide(3,2);

    for (size_t i=0; i<=ranges.size(); ++i) {
        int binL = h2d_fixed->GetXaxis()->FindBin(ranges[i].first);
        int binH = h2d_fixed->GetXaxis()->FindBin(ranges[i].second);
        TH1D* proj = h2d_fixed->ProjectionY(Form("proj_%zu",i), binL, binH);

        c_fix->cd(i+1);
        TF1* fit = GaussStuff::FitQGauss(proj, Form("fit_%zu",i));
        if (!fit) continue;

        double sigma_ns  = fit->GetParameter(2);
        double sigma_err = fit->GetParError(2);

        pe_index.push_back(i+1);
        sigma_ps.push_back(sigma_ns*1000.0);   // ps
        sigma_err_ps.push_back(sigma_err*1000.0);
    }

    
    // ====================
    // --- CALCOLO RAPPORTO BASELINE/SLOPE ---
    // ====================
    std::vector<double> pe_vals, ratio, ratio_err;

    for (int npe=1;npe<=5;++npe) {
        // --- Profili: slope della salita ---
        TProfile* prof = (TProfile*)f->Get("profile_classification_first_crossing/" + channel_name + Form("_class_%d_profile",npe));
        if (!prof) continue;

        int binmax = prof->GetMaximumBin();
        double maxVal = prof->GetBinContent(binmax);
        double maxPos = prof->GetBinCenter(binmax);

        TF1* fRise = FitRisingEdge(prof,maxPos,maxVal);
        if (!fRise) continue;

        double slope = fRise->GetParameter(1);
        double slope_err = fRise->GetParError(1);

        // --- Baseline: RMS ---
        TH1D* hbase = (TH1D*)f->Get("baseline/" + channel_name + Form("_class_%d_amp_baseline",npe));
        if (!hbase) continue;

        TF1* gfit = new TF1("gfit","gaus");
        hbase->Fit(gfit,"FQ");
        double sigma_base = gfit->GetParameter(2);

        // --- Rapporto e propagazione errore ---
        double val = sigma_base / slope;

        // Prendi anche l'errore sulla baseline
        double sigma_base_err = gfit->GetParError(2);

        // Propagazione corretta: errore combinato di due variabili indipendenti
        double err = val * std::sqrt(
            std::pow(sigma_base_err / sigma_base, 2) +
            std::pow(slope_err / slope, 2)
        );

        pe_vals.push_back(npe);
        ratio.push_back(val * 1000);
        ratio_err.push_back(err * 1000);

    }

    // ====================
    // --- Disegna entrambe le curve sulla stessa canvas ---
    // ====================
    TCanvas* c_compare = Graphics::CreateCanvas("c_compare","Resolution comparison",800,600);

    // --- Grafico Fixed Threshold ---
    TGraphErrors* gr_fixed = new TGraphErrors(pe_index.size(), pe_index.data(),
                                              sigma_ps.data(), nullptr, sigma_err_ps.data());
    Graphics::FormatAxis(gr_fixed->GetXaxis(), gr_fixed->GetYaxis(),
                         "photo-electrons", "#sigma_{t} (ps)",
                         0.03, 0.03, 62, false, 0, 0, true, 0, 120);
    Graphics::CustomizeGraph(gr_fixed, 47, kBlue-3, "APE", kBlue-3);
    gr_fixed->Draw("APE");  // Draw axis + points + error bars

    // --- Grafico Baseline / Slope ---
    TGraphErrors* gr_ratio = new TGraphErrors(pe_vals.size(), pe_vals.data(),
                                              ratio.data(), nullptr, ratio_err.data());
    Graphics::CustomizeGraph(gr_ratio, 46, kRed-3, "PE SAME", kRed-3);
    gr_ratio->Draw("P SAME");  // Draw points on same pad

    // --- Fit con parametri stimati dai dati ---
    // Parametri iniziali per gr_fixed
    double a_guess_fixed = sigma_ps.front() * std::sqrt(pe_index.front());
    double b_guess_fixed = sigma_ps.back();
    double xmin_fixed = pe_index.front();
    double xmax_fixed = pe_index.back();

    TF1* f_time_res_fixed = new TF1("f_time_res_fixed", "sqrt([0]*[0]/x + [1]*[1])", xmin_fixed, xmax_fixed);
    f_time_res_fixed->SetParNames("a","b");
    f_time_res_fixed->SetParameters(a_guess_fixed, b_guess_fixed);
    gr_fixed->Fit(f_time_res_fixed, "RQ"); // R=range, Q=quiet
    f_time_res_fixed->SetLineColor(kBlue-3);
    f_time_res_fixed->Draw("SAME");

    // Parametri iniziali per gr_ratio
    double a_guess_ratio = ratio.front() * std::sqrt(pe_vals.front());
    double b_guess_ratio = ratio.back();
    double xmin_ratio = pe_vals.front();
    double xmax_ratio = pe_vals.back();

    TF1* f_time_res_ratio = new TF1("f_time_res_ratio", "sqrt([0]*[0]/x + [1]*[1])", xmin_ratio, xmax_ratio);
    f_time_res_ratio->SetParNames("a","b");
    f_time_res_ratio->SetParameters(a_guess_ratio, b_guess_ratio);
    gr_ratio->Fit(f_time_res_ratio, "RQ");
    f_time_res_ratio->SetLineColor(kRed-3);
    f_time_res_ratio->Draw("SAME");

    // --- Legenda ---
    auto leg = Graphics::CreateLegend(0.6, 0.7, 0.88, 0.88);
    Graphics::AddLegendEntry(leg, gr_fixed, "Fixed threshold", "lp");
    Graphics::AddLegendEntry(leg, gr_ratio, "Baseline / Slope", "lp");
    leg->Draw();

}
