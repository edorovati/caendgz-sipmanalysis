void run_timing_analysis()
{
    std::vector<std::tuple<TString, TString, TString>> file_list;

    // Reference file
    file_list.emplace_back("scan_vbias_hama_gain/rooted/rooted_vbias_50.0_run_0.root", "reference", "50.0");

    // Data files from 51.0 to 56.8 in 0.2 steps
    for (float vb = 51.0; vb <= 57.0; vb += 0.2) {
        TString vb_str = Form("%.1f", vb);
        TString fname = Form("scan_vbias_hama_gain/rooted/rooted_vbias_%s_run_0.root", vb_str.Data());
        file_list.emplace_back(fname, "data", vb_str);
    }

    
    // Baseline window for each channel (edit channel names if needed)
    std::map<TString, std::pair<float, float>> baseline_windows;
    baseline_windows["ch1"] = std::make_pair(0.0, 30.0);

    // Call main function in starter.cpp
    starter(file_list, baseline_windows, "output_analysis.root");
}
