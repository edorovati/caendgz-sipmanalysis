void run_dcr_analysis()
{
    // Lista dei file: { percorso, vbias, tipo }
    std::vector<std::tuple<std::string, std::string, std::string>> file_list;

    // Qui aggiungi direttamente i file che vuoi processare
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_52.5.txt", "52.5", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_52.txt", "52", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_53.5.txt", "53.5", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_53.txt", "53", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_54.5.txt", "54.5", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_54.txt", "54", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_55.txt", "55", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_56.txt", "56", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_57.txt", "57", "data");
    file_list.emplace_back("dati/processed-hama/sn41-A1/dati/get_transitions_58.txt", "58", "data");
    
    // Nome del file di output ROOT
    TString outputFileName = "dcr_output.root";

    // Chiama la funzione principale definita nel tuo codice
    dcr(outputFileName, file_list, "custom-hama"); // "" se non vuoi specificare sensorName
}
