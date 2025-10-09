#ifndef DATA_EXTRACTOR_H
#define DATA_EXTRACTOR_H

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <map>
#include <iostream>
#include <TKey.h>
class DataExtractor {
public:
    explicit DataExtractor(TFile* file)
        : file_(file) {}

    bool extract(std::map<TString, TTree*>& tree_map,
                 std::map<TString, Double_t*>& time_map,
                 std::map<TString, Double_t*>& amp_map) const
    {
        if (!file_ || file_->IsZombie()) {
            std::cerr << "Errore: file ROOT non valido" << std::endl;
            return false;
        }

        // Estrae tutti gli alberi nel file ROOT
        TIter next(file_->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            TTree* tree = dynamic_cast<TTree*>(obj);
            if (!tree) continue;

            TString name = tree->GetName();
            tree_map[name] = tree;
            time_map[name] = Utils::setupBranch_dgz(tree, "time");
            amp_map[name]  = Utils::setupBranch_dgz(tree, "amplitude");
        }
        return true;
    }

private:
    TFile* file_;
};

#endif // DATA_EXTRACTOR_H
