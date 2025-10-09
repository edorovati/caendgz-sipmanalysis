#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TString.h>
#include <TError.h>     // gErrorIgnoreLevel
#include <iostream>
#include <vector>
#include <map>

/* ------------------------------------------------------------------
   1)  Funzione di utilità: cancella le funzioni/fit associati
------------------------------------------------------------------- */
void CleanHistogram(TH1* h)
{
    if (!h) return;
    if (h->GetListOfFunctions()) {
        h->GetListOfFunctions()->Delete();  // distrugge i TF1/TFormula
        h->GetListOfFunctions()->Clear();
    }
}

/* ------------------------------------------------------------------
   2)  Overload per sommare TH1F e TH2F
------------------------------------------------------------------- */
void AddHistogram(TH1F*& dest, TH1F* src)
{
    if (!src) return;
    CleanHistogram(src);                    // elimina fit dal sorgente
    if (!dest) {
        dest = (TH1F*)src->Clone();
        CleanHistogram(dest);               // assicura anche il clone pulito
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
   3)  Funzione principale
       Esempio:
       .L merge_histo_noFit.cpp
       merge_histos_noFit({"f1.root","f2.root"},"merged.root");
------------------------------------------------------------------ */
void merge_histo(const std::vector<TString>& input_files,
                        const TString&              output_file)
{
    std::map<TString, TH1F*> h1_map;
    std::map<TString, TH2F*> h2_map;

    // --‑‑ loop sui file di input --------------------------------------------------
    for (const auto& fname : input_files) {

        // opzionale: disattiva messaggi “Error” temporaneamente
        Int_t oldLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;          // ignora Warning+Error finché leggiamo

        TFile* file = TFile::Open(fname);
        gErrorIgnoreLevel = oldLevel;        // ripristina livello messaggi

        if (!file || file->IsZombie()) {
            std::cerr << "Errore nell’aprire " << fname << std::endl;
            continue;
        }

        TIter next(file->GetListOfKeys());
        TKey* key{};
        while ((key = (TKey*)next())) {
            TObject* obj  = key->ReadObj();
            TString   name = key->GetName();

            if      (obj->InheritsFrom("TH1F"))
                AddHistogram(h1_map[name], (TH1F*)obj);
            else if (obj->InheritsFrom("TH2F"))
                AddHistogram(h2_map[name], (TH2F*)obj);
            // qualunque altra cosa viene scartata (TF1, TGraph, ecc.)
        }
        file->Close();
    }

    // --‑‑ scrittura file di output -----------------------------------------------
    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Errore nella creazione di " << output_file << std::endl;
        return;
    }

    for (auto& [name, h1] : h1_map) if (h1) h1->Write(name);
    for (auto& [name, h2] : h2_map) if (h2) h2->Write(name);

    fout->Close();
    std::cout << "File '" << output_file
              << "' creato (istogrammi sommati, fit ignorati)." << std::endl;
}
