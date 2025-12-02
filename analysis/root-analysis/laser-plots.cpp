#include <TString.h>
#include <iostream>
#include "macros/fixed-th.cpp"
#include "macros/calibration.cpp"
#include "macros/probability.cpp"
#include "macros/crossing.cpp"
#include "macros/jitter.cpp"

void run_analysis(const char* mode, const char* filename, const char* channel_name, bool doFit = true) {
    gSystem->Load("libanalysis.so");
    TString mode_str(mode);
    TString file_str(filename);
    TString channel_str(channel_name);

    if (mode_str == "fixed_threshold-analysis") {
        std::cout << "Eseguo analisi fixed_threshold sul file " << file_str.Data()
                  << " per canale " << channel_str.Data() << std::endl;

        analyze_fixed_threshold(file_str, channel_str, doFit);
    }
    if (mode_str == "crossing") {
        std::cout << "Eseguo analisi fixed_threshold sul file " << file_str.Data()
                  << " per canale " << channel_str.Data() << std::endl;

        crossing(file_str, channel_str, doFit);
    }
    if (mode_str == "calibration") {
        std::cout << "Eseguo analisi fixed_threshold sul file " << file_str.Data()
                  << " per canale " << channel_str.Data() << std::endl;

        calibration(file_str, channel_str);
    }
    if (mode_str == "probability") {
        std::cout << "Eseguo analisi fixed_threshold sul file " << file_str.Data()
                  << " per canale " << channel_str.Data() << std::endl;

        probability(file_str, channel_str, doFit);
    }
    if (mode_str == "jitter") {
        std::cout << "Eseguo analisi fixed_threshold sul file " << file_str.Data()
                  << " per canale " << channel_str.Data() << std::endl;

        analyze_full(file_str, channel_str, doFit);
    } else {
        std::cerr << "Modalità sconosciuta: " << mode << std::endl;
        std::cerr << "Modalità supportata: fixed_threshold-analysis" << std::endl;
    }
}

void run_analysis_macro() {
    run_analysis("fixed_threshold-analysis", "arturo.root", "ch1");
}
