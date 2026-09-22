void batch_and_merge() {
    const int nFiles = 100;
    std::vector<TString> files;

    // Primo passo: batch timing
    for (int i = 0; i < nFiles; ++i) {
        TString input_file  = Form("dati/Highstat2/rooted_%d.root", i);
        TString output_file = Form("dati/Dati-hist/root%d.root", i);

        std::cout << "Processing " << input_file << " -> " << output_file << std::endl;
        timing(input_file, output_file);

        files.push_back(output_file); // salva il file per il merge
    }

    std::cout << "Batch timing finished for " << nFiles << " files." << std::endl;

    // Secondo passo: merge
    TString merged_file = "dati/mergedbsi-sptr.root";
    merge_histo(files, merged_file);
}

void do_merge() {
    std::vector<TString> files;

    const int nFiles = 100; // il numero di file che hai prodotto
    for (int i = 0; i < nFiles; ++i) {
        files.push_back(Form("dati/data-hist-hama/root%d.root", i));
    }

    TString output_file = "dati/data-hist-hama/merged.root";

    merge_histo(files, output_file);
}
