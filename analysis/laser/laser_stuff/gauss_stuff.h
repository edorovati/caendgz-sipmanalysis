#ifndef GAUSS_STUFF_H
#define GAUSS_STUFF_H

#include <cmath>
#include <vector>
#include <utility>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>

class GaussStuff {
public:
    static double meanPeak;



    static TF1* build_multi_gauss_fit_alternative(
        TH1F* hist,
        double amplitude_threshold_fraction = 0.005, // es. 0.5%
        double default_sigma = 0.5,
        double fit_window_mv = 2.0,
        int tail_bins = 10,
        int max_peaks = 10
    ) {
        const int nbins = hist->GetNbinsX();
        const double bin_width = hist->GetBinWidth(1);
        const double x_min = hist->GetXaxis()->GetXmin();
        const double x_max = hist->GetXaxis()->GetXmax();
        const double hist_max = hist->GetMaximum();
        const double total_entries = hist->GetEntries();

        // Smooth preliminare con media mobile (finestra 3 bin)
        std::vector<double> smooth_y(nbins+2, 0);
        for (int b = 1; b <= nbins; ++b) {
            double sum = 0;
            int count = 0;
            for (int w = -1; w <= 1; ++w) {
                int bb = b + w;
                if (bb >= 1 && bb <= nbins) {
                    sum += hist->GetBinContent(bb);
                    count++;
                }
            }
            smooth_y[b] = sum / count;
        }

        // Trova picchi locali (più permissivo), in numero massimo max_peaks
        std::vector<int> peak_bins;
        for (int b = 2; b < nbins - 2 && (int)peak_bins.size() < max_peaks; ++b) {
            double y1 = smooth_y[b - 1];
            double y2 = smooth_y[b];
            double y3 = smooth_y[b + 1];
            if (y2 > y1 && y2 > y3 && y2 > 0.01 * hist_max) {
                peak_bins.push_back(b);
                b += int(fit_window_mv / bin_width); // salta finestre vicine
            }
        }

        // Se non trovo abbastanza picchi, aggiungo posizioni a scaglioni
        while ((int)peak_bins.size() < max_peaks) {
            if (peak_bins.empty()) {
                for (int i = 0; i < max_peaks; ++i) {
                    double pos = x_min + (i + 0.5) * (x_max - x_min) / max_peaks;
                    int bin = hist->GetXaxis()->FindBin(pos);
                    peak_bins.push_back(bin);
                }
                break;
            } else {
                int last_bin = peak_bins.back();
                int step = int(fit_window_mv / bin_width);
                int next_bin = last_bin + step;
                if (next_bin > nbins) break;
                peak_bins.push_back(next_bin);
            }
        }

        // Costruisco la formula dinamicamente e aggiungo solo le gaussiane valide
        TString formula;
        int npeaks = 0;

        for (int i = 0; i < (int)peak_bins.size(); ++i) {
            double mu = hist->GetBinCenter(peak_bins[i]);
            double amp = hist->GetBinContent(peak_bins[i]);
            if (amp < 1.0) amp = 1.0; // evitare amp=0

            // Verifico se il bin ha entries > soglia relativa rispetto a total_entries
            double bin_content = hist->GetBinContent(peak_bins[i]);
            if (bin_content < amplitude_threshold_fraction * total_entries) {
                // Soglia non superata: fermo l'aggiunta di ulteriori gaussiane
                break;
            }

            if (npeaks > 0) formula += " + ";
            formula += Form("[%d]*exp(-0.5*((x-[%d])/[ %d ])^2)", 3*npeaks, 3*npeaks+1, 3*npeaks+2);
            npeaks++;
        }

        if (npeaks == 0) return nullptr; // Nessun picco valido trovato

        TF1* f_combined = new TF1("multi_gaus_fit_until_thresh", formula.Data(), x_min, x_max);
        f_combined->SetNpx(1000);
        f_combined->SetLineColor(kMagenta + 2);

        // Inizializzo parametri solo per i picchi usati
        for (int i = 0; i < npeaks; ++i) {
            double amp = hist->GetBinContent(peak_bins[i]);
            if (amp < 1.0) amp = 1.0;

            double mu = hist->GetBinCenter(peak_bins[i]);

            f_combined->SetParameter(3 * i + 0, amp);
            f_combined->SetParameter(3 * i + 1, mu);
            f_combined->SetParameter(3 * i + 2, default_sigma);

            f_combined->SetParName(3 * i + 0, Form("A_{%d p.e.}", i + 1));
            f_combined->SetParName(3 * i + 1, Form("#mu_{%d p.e.}", i + 1));
            f_combined->SetParName(3 * i + 2, Form("#sigma_{%d p.e.}", i + 1));

            f_combined->SetParLimits(3 * i + 0, 0, amp * 10);
            f_combined->SetParLimits(3 * i + 1, mu - fit_window_mv * 2, mu + fit_window_mv * 2);
            f_combined->SetParLimits(3 * i + 2, 0.1 * bin_width, 10.0);
        }

        // Definisco intervallo fit dinamico basato sull'ultimo picco incluso
        double fit_start = hist->GetXaxis()->GetBinLowEdge(1);
        int bin_end = hist->GetXaxis()->FindBin(f_combined->GetParameter(3 * (npeaks - 1) + 1));
        bin_end = std::min(bin_end + tail_bins, nbins);
        double fit_end = hist->GetXaxis()->GetBinUpEdge(bin_end);

        hist->Fit(f_combined, "IMREQ", "", fit_start, fit_end);

        return f_combined;
    }

    
    
};
#endif  // GAUSS_STUFF_H
