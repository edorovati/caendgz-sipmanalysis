#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TString.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <TError.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

/* ------------------------------------------------------------------
   Utility: pulisce funzioni/fit dagli istogrammi
------------------------------------------------------------------- */
void CleanHistogram(TH1* h)
{
    if (!h) return;
    if (h->GetListOfFunctions()) {
        h->GetListOfFunctions()->Delete();
        h->GetListOfFunctions()->Clear();
    }
}

/* ------------------------------------------------------------------
   Overload per sommare TH1F e TH2F
------------------------------------------------------------------- */
void AddHistogram(TH1F*& dest, TH1F* src)
{
    if (!src) return;
    CleanHistogram(src);
    if (!dest) {
        dest = (TH1F*)src->Clone();
        CleanHistogram(dest);
        dest->SetDirectory(0);
    } else {
        dest->Add(src);
    }
}

void AddHistogram(TH2F*& dest, TH2F* src)
{
    if (!src) return;
    CleanHistogram(src);
    if (!dest) {
        dest = (TH2F*)src->Clone();
        CleanHistogram(dest);
        dest->SetDirectory(0);
    } else {
        dest->Add(src);
    }
}

/* ------------------------------------------------------------------
   Ricorsiva: processa una directory e tutte le subdirectory
------------------------------------------------------------------- */
void ProcessDirectory(
    TDirectory* indir,
    TDirectory* outdir,
    std::map<std::string, TH1F*>& h1_map,
    std::map<std::string, TH2F*>& h2_map,
    std::map<std::string, std::vector<TGraph*>>& g_map,
    const std::string& prefix = "")
{
    if (!indir || !outdir) return;

    TIter next(indir->GetListOfKeys());
    TKey* key{};
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        std::string name = prefix + key->GetName();

        if (obj->InheritsFrom("TDirectory")) {
            // --- entra in subdirectory ---
            TDirectory* subIn  = (TDirectory*)obj;
            TDirectory* subOut = outdir->GetDirectory(key->GetName());
            if (!subOut) {
                subOut = outdir->mkdir(key->GetName());
            }
            ProcessDirectory(subIn, subOut, h1_map, h2_map, g_map, name + "/");
        }
        else if (obj->InheritsFrom("TH1F")) {
            AddHistogram(h1_map[name], (TH1F*)obj);
        }
        else if (obj->InheritsFrom("TH2F")) {
            AddHistogram(h2_map[name], (TH2F*)obj);
        }
        else if (obj->InheritsFrom("TGraph")) {
            g_map[name].push_back((TGraph*)obj->Clone());
        }
        // altro (TF1, ecc.) scartato
    }
}

/* ------------------------------------------------------------------
   Media dei TGraph: produce un TGraphErrors con media e errore
------------------------------------------------------------------- */
TGraphErrors* AverageGraphs(const std::vector<TGraph*>& graphs, const std::string& name)
{
    if (graphs.empty()) return nullptr;

    int nPoints = graphs[0]->GetN();
    auto gmean = new TGraphErrors(nPoints);
    gmean->SetName(name.c_str());
    gmean->SetTitle(name.c_str());

    for (int i = 0; i < nPoints; ++i) {
        double x0, y0;
        graphs[0]->GetPoint(i, x0, y0);

        std::vector<double> ys;
        ys.reserve(graphs.size());
        ys.push_back(y0);

        for (size_t j = 1; j < graphs.size(); ++j) {
            double xj, yj;
            graphs[j]->GetPoint(i, xj, yj);
            if (std::fabs(xj - x0) > 1e-9) {
                std::cerr << "Warning: x-values differ between graphs at point " << i
                          << " (" << x0 << " vs " << xj << ")\n";
            }
            ys.push_back(yj);
        }

        // calcola media
        double sum = 0;
        for (auto v : ys) sum += v;
        double mean = sum / ys.size();

        // deviazione standard campionaria
        double sqsum = 0;
        for (auto v : ys) sqsum += (v - mean) * (v - mean);
        double stdev = (ys.size() > 1) ? std::sqrt(sqsum / (ys.size() - 1)) : 0.0;

        // errore sulla media
        double err = (ys.size() > 0) ? stdev / std::sqrt((double)ys.size()) : 0.0;

        gmean->SetPoint(i, x0, mean);
        gmean->SetPointError(i, 0.0, err);
    }

    return gmean;
}

/* ------------------------------------------------------------------
   Funzione principale
------------------------------------------------------------------- */
void merge_histo(const std::vector<TString>& input_files,
                        const TString&       output_file)
{
    std::map<std::string, TH1F*> h1_map;
    std::map<std::string, TH2F*> h2_map;
    std::map<std::string, std::vector<TGraph*>> g_map;

    Long64_t total_waveforms = 0;  // accumulatore per n_waveforms

    // Loop sui file di input
    for (const auto& fname : input_files) {
        Int_t oldLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;
        TFile* file = TFile::Open(fname);
        gErrorIgnoreLevel = oldLevel;

        if (!file || file->IsZombie()) {
            std::cerr << "Errore nell’aprire " << fname << std::endl;
            continue;
        }

        // --- processa istogrammi e grafici ---
        ProcessDirectory(file, file, h1_map, h2_map, g_map);

       // --- leggi TTree Info ---
        TTree* tinfo = (TTree*)file->Get("info");
        if (tinfo) {
            Int_t n_wf = 0;  // <-- usa Int_t perché il branch originale è Int_t
            tinfo->SetBranchAddress("n_waveforms", &n_wf);
            if (tinfo->GetEntry(0) > 0) {
                total_waveforms += n_wf;  // accumulatore Long64_t
            }
        } else {
            std::cerr << "Warning: file " << fname << " non contiene TTree 'Info'\n";
        }


        file->Close();
    }

    // Scrittura file di output
    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Errore nella creazione di " << output_file << std::endl;
        return;
    }

    // Istogrammi
    for (auto& [fullname, h1] : h1_map) {
        if (!h1) continue;
        std::string path = fullname.substr(0, fullname.find_last_of('/'));
        std::string name = fullname.substr(fullname.find_last_of('/') + 1);
        TDirectory* dir = fout;
        if (!path.empty()) dir = fout->mkdir(path.c_str(), "", true);
        dir->cd();
        h1->Write(name.c_str());
    }
    for (auto& [fullname, h2] : h2_map) {
        if (!h2) continue;
        std::string path = fullname.substr(0, fullname.find_last_of('/'));
        std::string name = fullname.substr(fullname.find_last_of('/') + 1);
        TDirectory* dir = fout;
        if (!path.empty()) dir = fout->mkdir(path.c_str(), "", true);
        dir->cd();
        h2->Write(name.c_str());
    }

    // Grafici
    for (auto& [fullname, vec] : g_map) {
        if (vec.empty()) continue;
        TGraphErrors* gmean = AverageGraphs(vec, fullname);
        if (!gmean) continue;

        std::string path = fullname.substr(0, fullname.find_last_of('/'));
        std::string name = fullname.substr(fullname.find_last_of('/') + 1);
        TDirectory* dir = fout;
        if (!path.empty()) dir = fout->mkdir(path.c_str(), "", true);
        dir->cd();
        gmean->Write(name.c_str());
    }

    // TTree Info con somma n_waveforms
    fout->cd();
    TTree* tinfo_out = new TTree("info", "Merged info tree");
    Int_t n_waveforms_out = (Int_t) total_waveforms;  // truncato se > 2 miliardi
    tinfo_out->Branch("n_waveforms", &n_waveforms_out, "n_waveforms/I");
    tinfo_out->Fill();
    tinfo_out->Write();


    fout->Close();
    std::cout << "File '" << output_file
              << "' creato (istogrammi sommati, grafici mediati, n_waveforms sommati)."
              << std::endl;
}
