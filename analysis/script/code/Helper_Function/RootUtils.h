/******************************************************************************
 * @file RootUtils.h
 * @brief Declaration of the RootUtils class: static helper functions for loading
 *        TTrees from ROOT files and drawing explanatory legends on canvases.
 *
 * This header provides:
 *   - RootUtils::loadTreeFromFile     : Load a TTree from a ROOT file
 *   - RootUtils::drawExplanationLegend: Draw a legend explaining marker and band
 ******************************************************************************/
#ifndef ROOTUTILS_H
#define ROOTUTILS_H

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TBox.h"
#include "TLegend.h"
#include <iostream>

/******************************************************************************
 *                               RootUtils Class
 * ----------------------------------------------------------------------------
 * @class RootUtils
 * @brief Collection of static utility methods for ROOT I/O and plotting helpers.
 *
 * This class includes the following public static functions:
 * ----------------------------------------------------------------------------
 * | **Public Functions**               | **Description**                          |
 * ----------------------------------------------------------------------------
 * | **loadTreeFromFile**               | Opens a ROOT file, navigates to a       |
 * |                                    | specified directory, and retrieves a     |
 * |                                    | TTree by name.                          |
 * ----------------------------------------------------------------------------
 * | **drawExplanationLegend**          | Draws a TPaveText legend showing        |
 * |                                    | markers and confidence bands            |
 * ----------------------------------------------------------------------------
 ******************************************************************************/
class RootUtils {
public:
    //--------------------------------------------------------------------------------
    /// @brief Load a TTree from a ROOT file and return success status.
    /// @param filePath      Path to the ROOT file.
    /// @param treeName      Name of the TTree to load.
    /// @param outTree       Reference to TTree* to receive the pointer.
    /// @param directoryName Name of subdirectory (default "Summary").
    /// @return True if tree is successfully loaded; false otherwise.
    //--------------------------------------------------------------------------------
    static bool loadTreeFromFile(const TString& filePath,
                                 const char* treeName,
                                 TTree*& outTree,
                                 const char* directoryName = "Summary") {
        TFile* file = TFile::Open(filePath);
        if (!file || file->IsZombie()) {
            std::cerr << "WARNING: Couldn't open file: " << filePath << "\n";
            return false;
        }

        if (file->GetListOfKeys()->Contains(directoryName)) {
            TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(directoryName));
            if (dir) {
                dir->cd();
            } else {
                std::cerr << "WARNING: Directory '" << directoryName
                          << "' is not a valid TDirectory in file: " << filePath << "\n";
                file->Close();
                return false;
            }
        } else {
            std::cout << "INFO: Directory '" << directoryName
                      << "' not found, reading from file root.\n";
            file->cd();
        }

        outTree = dynamic_cast<TTree*>(gDirectory->Get(treeName));
        if (!outTree) {
            std::cerr << "WARNING: Tree '" << treeName
                      << "' not found in file: " << filePath << "\n";
            file->Close();
            return false;
        }

        return true;
    }

    //--------------------------------------------------------------------------------
    /// @brief Draw an explanatory legend showing single p.e. marker and band.
    /// @param x1,y1,x2,y2 NDC coordinates for the legend box.
    //--------------------------------------------------------------------------------
    static void drawExplanationLegend(double x1 = 0.15,
                                      double y1 = 0.6,
                                      double x2 = 0.5,
                                      double y2 = 0.7) {
        TGraph* markerGraph = new TGraph();
        markerGraph->SetMarkerStyle(20);
        markerGraph->SetMarkerColor(kBlack);
        markerGraph->SetMarkerSize(1);

        TBox* bandBox = new TBox(0, 0, 1, 1);
        bandBox->SetFillColorAlpha(kGray, 0.3);
        bandBox->SetLineColor(kGray);

        TLegend* legend = new TLegend(x1, y1, x2, y2);
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetFillColor(kWhite);
        legend->SetTextFont(62);
        legend->SetTextSize(0.035);
        legend->SetLineColor(kBlack);
        legend->AddEntry(markerGraph, "single p.e. (mV)", "p");
        legend->AddEntry(bandBox, "single p.e. +/- sigma (mV)", "f");
        legend->Draw();
    }
};

#endif // ROOTUTILS_H
