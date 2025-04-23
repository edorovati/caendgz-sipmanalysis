#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <vector>

struct BandPoint {
    double x;
    double yTop;
    double yBot;
};

class GraphicsUtils {
public:

    static void setGraphStyle(TGraphErrors* graph, int markerStyle, int color, int lineStyle = 1) {
        graph->SetMarkerStyle(markerStyle);
        graph->SetLineColor(color);
        graph->SetMarkerColor(color);
        graph->SetLineStyle(lineStyle);
    }
    
    static TCanvas* createCanvas(const char* name, const char* title, int width, int height, bool logY = false) {
        TCanvas* canvas = new TCanvas(name, title, width, height);
        if (logY) {
            canvas->SetLogy();
        }
        canvas->SetGrid();
        return canvas;
    }
    static void setAxisAndDraw(TMultiGraph* mg, const char* xTitle, const char* yTitle, bool drawAPE = true) {
        mg->SetTitle(Form("%s; %s; %s", mg->GetTitle(), xTitle, yTitle));
        mg->GetXaxis()->SetTitleFont(62);
        mg->GetYaxis()->SetTitleFont(62);
    if (drawAPE) {
            mg->Draw("APE");
        } else {
            mg->Draw("ALE");
        }
    }

    static TLegend* createLegend(double x1, double y1, double x2, double y2, bool noFill = true) {
        TLegend* legend = new TLegend(x1, y1, x2, y2);
        if (noFill) {
            legend->SetFillStyle(0);
            legend->SetLineStyle(0);
            legend->SetBorderSize(0);
        }
        legend->SetTextFont(62);
        legend->SetTextSize(0.03);
        return legend;
    }

    static TGraph* createBandGraph(const std::vector<BandPoint>& bandData) {
        const int bandN = 2 * bandData.size() + 1;
        TGraph* gBand = new TGraph(bandN);
        
        for (int i = 0; i < bandData.size(); ++i) {
            gBand->SetPoint(i, bandData[i].x, bandData[i].yTop);
        }

        for (int i = 0; i < bandData.size(); ++i) {
            auto& pt = bandData[bandData.size() - 1 - i];
            gBand->SetPoint(bandData.size() + i, pt.x, pt.yBot);
        }

        gBand->SetPoint(bandN - 1, bandData[0].x, bandData[0].yTop);

        gBand->SetFillColorAlpha(kBlue, 0.2);
        gBand->SetLineColor(kBlue);
        gBand->SetLineStyle(2);

        return gBand;
    }
    static void drawGraphWithBand(TGraphErrors* graph, const std::vector<BandPoint>& bandData,
                                      const char* title, const char* xTitle, const char* yTitle, const char* canvasName) {
            
            TCanvas* canvas = createCanvas(canvasName, title, 800, 600);
            TGraph* gBand = createBandGraph(bandData);
            TMultiGraph* mg = new TMultiGraph();
            mg->Add(gBand, "F");
            mg->Add(graph, "PE");
            setAxisAndDraw(mg, xTitle, yTitle);
            canvas->Draw();
        }
    
    static void addEntryToLegend(TLegend* legend, TGraph* graph, const char* label, const char* type) {
            if (strcmp(type, "lep") == 0) {
                legend->AddEntry(graph, label, "lep");
            } else if (strcmp(type, "f") == 0) {
                legend->AddEntry(graph, label, "f");
            } else {
                std::cerr << "Tipo di entry non valido: " << type << std::endl;
            }
        }
    
        static TPaveText* createPaveText(double x1, double y1, double x2, double y2, const std::vector<std::string>& lines) {
            TPaveText* pave = new TPaveText(x1, y1, x2, y2, "NDC");
            pave->SetFillColor(0);
            pave->SetTextFont(42);
            pave->SetTextAlign(12);
            for (const auto& line : lines) {
                pave->AddText(line.c_str());
            }

            return pave;
        }
};

#endif 
