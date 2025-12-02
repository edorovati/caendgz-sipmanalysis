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


// Funzione principale
void probability(const TString& filename, const TString& channel_name, bool doFit = true) {
    // --- Apri file ROOT ---
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore apertura file: " << filename << std::endl;
        return;
    }

    // ---- Accedi alla cartella fixed_threshold ----
    f->cd("fixed_threshold");
    TString hist_name = channel_name + "_fixed_th_diff_laser";
    TH2F* h2d = (TH2F*)gDirectory->Get(hist_name);
    if (!h2d) {
        std::cerr << "Non trovo istogramma " << hist_name << std::endl;
        return;
    }

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
           TF1* linearFit = new TF1("linearFit", "pol1");
           gr_cal->Fit(linearFit, "Q"); // "Q" = quiet
           slope     = linearFit->GetParameter(1);
           intercept = linearFit->GetParameter(0);
       } else {
           // ---- Prendi solo il primo punto del TGraph ----
           if (gr_cal->GetN() > 1) {
               double x0, y0;
               gr_cal->GetPoint(1, x0, y0);
               slope = y0;
               intercept = 0;
           } else {
               std::cerr << "TGraph vuoto, non posso prendere il primo punto." << std::endl;
           }
       }


    // ---- Recupera n_waveforms dal TTree infoTTree ----
    f->cd();
    TTree* infoTree = (TTree*)f->Get("info");
    Int_t n_wf = 0;
    if (infoTree) {
        infoTree->SetBranchAddress("n_waveforms", &n_wf);
        infoTree->GetEntry(0);
    } else {
        std::cerr << "Attenzione: infoTTree non trovato, n_waveforms = 0" << std::endl;
    }

    // ---- Calcola i tagli ----
    double low_cut  = 0.5*slope + intercept;
    double high_cut = 1.5*slope + intercept;
    int bin_lowcut  = h2d->GetXaxis()->FindBin(low_cut);
    int bin_highcut = h2d->GetXaxis()->FindBin(high_cut);
    int bin_max     = h2d->GetXaxis()->GetNbins();

    // ---- Proiezioni ----
    TH1D* hist_above_low  = h2d->ProjectionY("events > 0.5 p.e.", bin_lowcut, bin_max);
    TH1D* hist_above_high = h2d->ProjectionY("events > 1.5 p.e.", bin_highcut, bin_max);

    // Trova il bin massimo per ogni istogramma
    int b_low = Utils::FindMaxBin(hist_above_low);
    int b_high = Utils::FindMaxBin(hist_above_high);

    // Calcola S e B per entrambi
    Analysis::SignalBackground sb_low = Analysis::CalculateSignalAndBackground(hist_above_low, b_low);
    Analysis::SignalBackground sb_high = Analysis::CalculateSignalAndBackground(hist_above_high, b_high);

    // Calcola crosstalk e errore
    auto [ct, ct_err] = Analysis::CalculateCrosstalk(sb_high.diff,sb_low.diff);
    
    // Tempo del bin massimo
    double t_peak = hist_above_low->GetXaxis()->GetBinCenter(b_low);

    // Tempo massimo per la finestra dark = picco - 2 ns
    double t_end_dark = t_peak - 2.0;  // ns

    // Trova il primo bin con eventi
   
    // Bin di start (a -30 ns invece che dal primo bin con eventi)
    int bin_start_dark = hist_above_low->GetXaxis()->FindBin(-30.0);

    // Se per qualche motivo è 0 o fuori range, fallback
    if (bin_start_dark <= 0) {
        std::cerr << "Warning: bin di start a -30 ns non valido, fallback a 1!\n";
        bin_start_dark = 1;
    }

    // Bin di stop
    int bin_end_dark = hist_above_low->GetXaxis()->FindBin(t_end_dark);

    // Conta eventi nel range
    double dark_counts = Utils::CountEvents(hist_above_low, bin_start_dark, bin_end_dark);

    // Calcolo della finestra temporale
    double time_min = hist_above_low->GetXaxis()->GetBinLowEdge(bin_start_dark);
    double time_window = t_end_dark + fabs(time_min);

    // Calcola DCR
    auto [dcr, dcr_err] = Analysis::CalculateDCR(dark_counts, time_window * 1e-9, n_wf);

    // Stampa di tutte le variabili rilevanti
    std::cout << "=== DEBUG Dark Counts ===" << std::endl;
    std::cout << "t_peak       : " << t_peak << " ns" << std::endl;
    std::cout << "t_end_dark   : " << t_end_dark << " ns" << std::endl;
    std::cout << "bin_start_dark : " << bin_start_dark << std::endl;
    std::cout << "bin_end_dark   : " << bin_end_dark << std::endl;
    std::cout << "time_min     : " << time_min << " ns" << std::endl;
    std::cout << "time_window  : " << time_window << " ns" << std::endl;
    std::cout << "dark_counts  : " << dark_counts << std::endl;
    std::cout << "DCR          : " << dcr << " Hz" << std::endl;
    std::cout << "DCR error    : " << dcr_err << " Hz" << std::endl;
    std::cout << "=========================" << std::endl;

    
    // Calcolo probabilità di segnale per hist_above_low (diff S-B) / n_wf
    auto [prob_low, prob_err_low] = Analysis::CalculateProbability(sb_low.diff, n_wf);
    auto [prob_high, prob_err_high] = Analysis::CalculateProbability(sb_high.diff, n_wf);

   
    // ---- Canvas unica con entrambi gli istogrammi ----
    TCanvas* c_thresholds = Graphics::CreateCanvas("c_thresholds", "Eventi sopra soglie", 800, 600, false, true, true);
    Graphics::FormatAxis(hist_above_low->GetXaxis(), hist_above_low->GetYaxis(),
                         "#Deltat", "counts", 0.03, 0.03, 62, true, 30, 50);

    Graphics::SetupHistogramStyle(hist_above_low, false, 1111,
                                   kBlue+2, kBlue+2, 47, 2, "E1", false);
    
    /*Graphics::SetupHistogramStyle(hist_above_high, false, 1111,
                                   kRed+1, kRed+1, 46, 2, "E1 SAME", false);*/
   
    TLegend* leg = Graphics::CreateLegend(0.1, 0.9, 0.3, 0.7);
    Graphics::AddLegendEntry(leg, hist_above_low, "> 0.5 p.e.", "lpe");
    //Graphics::AddLegendEntry(leg, hist_above_high, "> 1.5 p.e.", "lpe");
    leg->Draw();
    TString infoText;
    infoText += Form("CT: %.1f %% #pm %.1f %%\n", ct, ct_err);
    infoText += Form("DCR: %.0f kHz #pm %.0f kHz\n", dcr, dcr_err);
    infoText += Form("P_{signal}: %.2f #pm %.2f %%\n", prob_low*100, prob_err_low*100);
    
    TPaveText* infoBox = Graphics::CreateInfoBox(0.6, 0.6, 0.9, 0.9, infoText, kWhite, kBlack, 0.035, 1);
    infoBox->Draw();

    std::vector<double> npe_values;
    std::vector<double> prob_values;
    std::vector<double> prob_err_values;

    double sum_prob_nonzero = 0.0;
    // ---- Loop per intervalli multipli di slope ----
    std::vector<std::pair<double,double>> ranges;
    for (int i = 0; i < 10; ++i) {
        double start = (0.5 + i) * slope + intercept;
        double end   = (1.5 + i) * slope + intercept;
        if (end > h2d->GetXaxis()->GetXmax()) break;
        ranges.push_back({start, end});
    }

    // Primo ciclo: conto quanti intervalli hanno eventi
    int n_valid = 0;
    for (auto& r : ranges) {
        int bin_start = h2d->GetXaxis()->FindBin(r.first);
        int bin_end   = h2d->GetXaxis()->FindBin(r.second);
        TH1D* proj = h2d->ProjectionY(Form("%.2f-%.2f p.e.", r.first, r.second), bin_start, bin_end);

        if (proj->Integral() == 0) {
            break;  // se zero, esco subito (come richiesto)
        }
        n_valid++;
    }
    c_thresholds->Update();
    
    // Secondo ciclo: disegno solo i validi
    // Array di colori ben distinti
    const Color_t colors[] = {
        kBlue+2,    // Blu
        kMagenta+1, // Magenta
        kRed+1,     // Rosso
        kOrange+7,  // Arancio
        kGreen+2,   // Verde
        kViolet+1,  // Viola
        kCyan+1,    // Ciano
        kPink+9,    // Rosa
        kTeal+3,    // Verde acqua
        kAzure+4    // Blu chiaro
    };

    // Array di marker diversi
    const Style_t markerStyles[] = {45, 46, 47, 28, 34, 40, 41, 48, 30, 33};

    int cols = TMath::CeilNint(TMath::Sqrt(n_valid));
    int rows = TMath::CeilNint((double)n_valid / cols);
    TCanvas* c_ranges = Graphics::CreateCanvas("c_ranges", "Eventi per intervallo", 1500, 800);
    c_ranges->Divide(cols, rows);

    int pad_idx = 1;
    for (size_t i = 0; i < ranges.size(); ++i) {
        auto& r = ranges[i];
        int bin_start = h2d->GetXaxis()->FindBin(r.first);
        int bin_end   = h2d->GetXaxis()->FindBin(r.second);
        TH1D* proj = h2d->ProjectionY(Form("%zu p.e.", i+1), bin_start, bin_end);

        if (proj->Integral() == 0) break;

        int b_proj = Utils::FindMaxBin(proj);
        Analysis::SignalBackground sb_proj = Analysis::CalculateSignalAndBackground(proj, b_proj);
        auto [prob, prob_err] = Analysis::CalculateProbability(sb_proj.diff, n_wf);

        c_ranges->cd(pad_idx++);
        Graphics::FormatAxis(proj->GetXaxis(), proj->GetYaxis(),
                             "#Deltat", "counts", 0.03, 0.03, 62, true, 30, 50);

        // Uso colore e marker diversi per ogni i
        int colorIndex = i % (sizeof(colors)/sizeof(colors[0]));
        int markerIndex = i % (sizeof(markerStyles)/sizeof(markerStyles[0]));

        Graphics::SetupHistogramStyle(
            proj,
            true,               // fill?
            1111,               // line style
            colors[colorIndex], // line color
            colors[colorIndex], // marker color
            markerStyles[markerIndex], // marker style
            2,                  // marker size
            "HISTE1",           // draw option
            true,               // fill area?
            colors[colorIndex], // fill color
            3354,               // fill style (pattern carino)
            0.35, 0.6, 0.9,     // legend box x1,y1,x2
            0.6, 0.9            // legend box y1,y2
        );

        Graphics::SetupPad(false, true, true);

        TString infoText;
        infoText += Form("P_{%zu p.e.}: %.2f #pm %.2f %%", i+1, prob*100.0, prob_err*100.0);

        TPaveText* infoBox = Graphics::CreateInfoBox(0.1, 0.9, 0.4, 0.8,
                                                     infoText, kWhite, kBlack, 0.05, 1);
        infoBox->Draw();
        double prob_percent = prob * 100.0;
        double prob_err_percent = prob_err * 100.0;

        npe_values.push_back(i+1);
        prob_values.push_back(prob_percent);
        prob_err_values.push_back(prob_err_percent);
        sum_prob_nonzero += prob_percent;
    }

    c_ranges->Update();
    
    // --- Calcolo somma errore probabilità non-zero in percentuale ---
    double sum_prob_nonzero_err_sq = 0.0;
    for (size_t i = 0; i < prob_err_values.size(); ++i) {
        sum_prob_nonzero_err_sq += prob_err_values[i] * prob_err_values[i];
    }
    double sum_prob_nonzero_err = std::sqrt(sum_prob_nonzero_err_sq);

    // ---- Calcolo P(0 p.e.) e µ con errore propagato ----
    double P0_percent = 100.0 - sum_prob_nonzero;       // in percentuale
    double P0 = P0_percent / 100.0;                      // frazione (0..1)

    double P0_percent_err = sum_prob_nonzero_err;       // errore in percentuale
    double P0_err = P0_percent_err / 100.0;             // errore frazionario

    // Protezione contro P0 <= 0 (log non definito)
    if (P0 <= 0.0) {
        std::cerr << "Errore: P0 <= 0, impossibile calcolare mu" << std::endl;
        P0 = 1e-6;
    }

    double mu_from_count = -std::log(P0);

    // Propagazione errore su mu: sigma_mu = sigma_P0 / P0
    double mu_from_count_err = P0_err / P0;

    // --- Canvas ---
    TCanvas* c_prob_vs_pe = Graphics::CreateCanvas("c_prob_vs_pe", "Probability vs Photo-electrons", 1000, 700,false, true, true);

    // Grafico con errori
    TGraphErrors* g_prob = new TGraphErrors(npe_values.size() + 1);

    // Punto a 0 p.e.
    g_prob->SetPoint(0, 0, P0_percent);
    g_prob->SetPointError(0, 0.0, P0_percent_err);

    // Altri punti (shiftati di +1 nell’indice)
    for (size_t i = 0; i < npe_values.size(); ++i) {
        g_prob->SetPoint(i+1, npe_values[i], prob_values[i]);
        g_prob->SetPointError(i+1, 0.0, prob_err_values[i]);
    }

    // Personalizza assi e stile grafico
    Graphics::CustomizeGraph(g_prob, 47, kOrange+7, "APE", kOrange+7, 1, 1);

    // Poisson teorica dal conto
    TF1* f_poisson_theory = new TF1("f_poisson_theory",
        "[0]*TMath::PoissonI(x, [1])", 0, npe_values.size()+1);
    f_poisson_theory->SetParameters(100.0, mu_from_count);
    f_poisson_theory->SetLineColor(kBlack);
    f_poisson_theory->SetLineStyle(6);
    f_poisson_theory->SetLineWidth(4);
    f_poisson_theory->Draw("SAME");

    // Fit Poisson ai dati
    TF1* f_poisson_fit = new TF1("f_poisson_fit",
        "[0]*TMath::PoissonI(x, [1])", 0, npe_values.size()+1);
    f_poisson_fit->SetParameters(100.0, 1.0);

    g_prob->Fit(f_poisson_fit, "Q");
    double mu_from_fit     = f_poisson_fit->GetParameter(1);
    double mu_from_fit_err = f_poisson_fit->GetParError(1);

    f_poisson_fit->SetLineColor(kRed);
    f_poisson_fit->SetLineStyle(6);
    f_poisson_fit->SetLineWidth(4);

    // Legenda
    TLegend* leg_prob = Graphics::CreateLegend(0.1, 0.25, 0.35, 0.40);
    Graphics::AddLegendEntry(leg_prob, g_prob, "p.e. distribution", "lep");
    Graphics::AddLegendEntry(leg_prob, f_poisson_theory,"poisson estimation", "l");
    Graphics::AddLegendEntry(leg_prob, f_poisson_fit, "poisson fit", "l");
    leg_prob->Draw();

    // Info box
    TString infoText_prob;
    infoText_prob += Form("P(0): #mu = %.2f #pm %.2f\n", mu_from_count, mu_from_count_err);
    infoText_prob += Form("P(fit):   #mu = %.2f #pm %.2f\n", mu_from_fit, mu_from_fit_err);

    TPaveText* infoBox_prob = Graphics::CreateInfoBox(0.1, 0.1, 0.35, 0.25,
                                                      infoText_prob, kWhite, kBlack, 0.035, 1);
    infoBox_prob->Draw();

    Graphics::FormatAxis(
        g_prob->GetXaxis(),
        g_prob->GetYaxis(),
        "# photo-electrons",
        "probability (%)", 0.03,0.03,62, false, 0,0,true, 0, 100,42);

}
