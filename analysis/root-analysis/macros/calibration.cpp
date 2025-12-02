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

void calibration(const TString& filename, const TString& ch_name) {
    TFile* f = TFile::Open(filename, "UPDATE");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore apertura file: " << filename << std::endl;
        return;
    }

    TDirectory* dir_amp = f->GetDirectory("amp_max-distribution");
    if (!dir_amp) {
        std::cerr << "Directory amp_max-distribution non trovata nel file!" << std::endl;
        f->Close();
        return;
    }

    TDirectory* dir_cal = f->GetDirectory("calibration");
    if (!dir_cal) {
        std::cerr << "Directory calibration non trovata nel file!" << std::endl;
        f->Close();
        return;
    }

    dir_amp->cd();
    

    // cerca istogramma che contiene sia ch_name che "_diff"
    TH1* hist = nullptr;
    TIter next(dir_amp->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TString keyName = key->GetName();
        if (keyName.Contains(ch_name) && keyName.Contains("_signal")) {
            hist = dynamic_cast<TH1*>(key->ReadObj());
            break;
        }
    }

    if (!hist) {
        std::cerr << "Non trovo istogramma che contiene '"
                  << ch_name << "' and '_diff'" << std::endl;
        f->Close();
        return;
    }

    

    TH1F* histF = dynamic_cast<TH1F*>(hist);
    if (!histF) {
        std::cerr << "Istogramma non è di tipo TH1F, impossibile fare il fit!" << std::endl;
        f->Close();
        return;
    }

    TString fit_name = ch_name + "-multi-gauss";
    TF1* fit = GaussStuff::multi_gauss(histF, 5, fit_name);

    TCanvas* c_amp = Graphics::CreateCanvas("c_amp", "Amp Max Distribution", 1000, 700, false, true, true);
    Graphics::FormatAxis(hist->GetXaxis(), hist->GetYaxis(), "amplitude (mV)", "counts", 0.03, 0.03, 62, true, 0, 50);
    Graphics::SetupHistogramStyle(hist, true, 111111, kMagenta-2, kBlue+2, 46, 2, "HISTE", true, kMagenta-2, 3149, 0.35);
    fit->SetLineColor(kRed);
    fit->SetLineStyle(1);
    fit->SetLineWidth(4);
    fit->Draw("SAME");
    c_amp->Update();
    // Numero di picchi fittati
    const int n_peaks = 6;

    // Vettori per i valori
    std::vector<double> pe_values;   // fotoelettroni (1, 2, ..., n_peaks)
    std::vector<double> mean_values; // medie estratte dal fit
    std::vector<double> mean_errors; // errori sulle medie (dai parametri fit)
    std::vector<double> pe_errors(n_peaks, 0); // errori sugli x (zero se non disponibili)

    for (int i = 0; i < n_peaks; ++i) {
        double pe = i;  // assumiamo picchi da 1 a 5 pe
        pe_values.push_back(pe);

        int mean_idx = 1 + 3*i;  // indice parametro media gaussiana i-esima
        double mean = fit->GetParameter(mean_idx);
        double mean_err = fit->GetParError(mean_idx);

        mean_values.push_back(mean);
        mean_errors.push_back(mean_err);
    }

    // Costruisci TGraphErrors
    TGraphErrors* graph = new TGraphErrors(n_peaks,
                                          pe_values.data(),
                                          mean_values.data(),
                                          pe_errors.data(),
                                          mean_errors.data());
    graph->SetName(ch_name + "_mean_vs_pe");
    graph->SetTitle("Mean amplitude vs photo-electrons");

    TF1* linear_fit = new TF1("linear_fit", "[0] + [1]*x", 1, pe_values.back());
    linear_fit->SetParNames("q", "m");  // q = intercept, m = slope
    linear_fit->SetLineColor(kMagenta+2);
    linear_fit->SetLineStyle(6);
    linear_fit->SetLineWidth(2);
    graph->Fit(linear_fit, "IMREQ");
    double q = linear_fit->GetParameter(0); // intercetta
    double m = linear_fit->GetParameter(1); // slope


    // Disegna il grafico
    TCanvas* c_mean = Graphics::CreateCanvas("c_mean", "Mean vs Photo-electrons", 800, 600, false, false, true);
    Graphics::FormatAxis(graph->GetXaxis(), graph->GetYaxis(), "#photo-electrons", "amplitude (mV)", 0.03, 0.03);
    Graphics::CustomizeGraph(graph, 47, kMagenta-2, "APE", kMagenta-2);
    graph->Draw("APE");

    c_mean->Update();

    c_amp->cd();  // torna sulla canvas dell'istogramma

    int n_lines = n_peaks+10;
    for (int i = 0; i < n_lines; ++i) {
        double x_line = q + (i + 0.5) * m;
        TLine* line = new TLine(x_line, 0, x_line, histF->GetMaximum());
        line->SetLineColor(kOrange+7);
        line->SetLineStyle(1);
        line->SetLineWidth(2);
        line->Draw("SAME");
    }
    gPad->Update();          // aggiorna il pad corrente
    TPaveStats* stats = (TPaveStats*)histF->GetListOfFunctions()->FindObject("stats");
    if (stats) stats->Draw();  // ridisegna il box

    c_amp->Modified();
    c_amp->Update();

    
    // Calcolo della soglia: 0.5*slope + intercept
    double threshold = q + 0.5 * m;

    // Integrale totale dell'istogramma
    double total_integral = histF->Integral();

    // Trova il bin corrispondente alla soglia
    int bin_threshold = histF->FindBin(threshold);

    // Integrale dalla soglia fino alla fine
    double integral_above = histF->Integral(bin_threshold, histF->GetNbinsX());

    // Rapporto
    double fraction = integral_above / total_integral;

    // Errore statistico binomiale
    int n_entries = histF->GetEntries();  // numero totale di eventi (non normalizzato)
    double fraction_err = std::sqrt(fraction * (1.0 - fraction) / n_entries);

    // Stampa a terminale
    std::cout << "Soglia (0.5*slope + intercept) = " << threshold << " mV" << std::endl;
    std::cout << "Integrale totale = " << total_integral << std::endl;
    std::cout << "Integrale sopra soglia = " << integral_above << std::endl;
    std::cout << "Frazione sopra soglia = " << fraction
              << " ± " << fraction_err << std::endl;


}
