#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <TCanvas.h>
#include <TH1.h>
#include <TLine.h>
#include <TAttAxis.h>
#include <TString.h>

class Graphics {
public:
    // ====== Crea Canvas ======
    static TCanvas* CreateCanvas(
        const TString& name = "c1",
        const TString& title = "",
        int width = 1200,
        int height = 800,
        bool logx = false,
        bool logy = false,
        bool grid = false
    ) {
        TCanvas* c = new TCanvas(name, title, width, height);
        c->SetLogx(logx);
        c->SetLogy(logy);
        if (grid) {
            c->SetGrid();
        }
        return c;
    }

    static void SetupPad(bool logx = false, bool logy = false, bool grid = false, bool logz = false) {
        if (logx) gPad->SetLogx();
        else gPad->SetLogx(0);

        if (logy) gPad->SetLogy();
        else gPad->SetLogy(0);

        if (grid) {
            gPad->SetGridx(1);
            gPad->SetGridy(1);
        } else {
            gPad->SetGridx(0);
            gPad->SetGridy(0);
        }

        
        if (logz) gPad->SetLogz();
        else gPad->SetLogz(0);
    }

    
    // formatta Assi per tutti ======
    static void FormatAxesGraph(
        TGraph* gr,
        const TString& xTitle,
        const TString& yTitle,
        double titleSize = 0.03,
        double labelSize = 0.03,
        int font = 62,              // 42=Helvetica, 62=Helvetica bold
        bool setXLimits = false,
        double xMin = 0,
        double xMax = 0,
        bool setYLimits = false,
        double yMin = 0,
        double yMax = 0
    ) {
        gr->GetXaxis()->SetTitle(xTitle);
        gr->GetYaxis()->SetTitle(yTitle);

        gr->GetXaxis()->SetTitleSize(titleSize);
        gr->GetYaxis()->SetTitleSize(titleSize);

        gr->GetXaxis()->SetLabelSize(labelSize);
        gr->GetYaxis()->SetLabelSize(labelSize);

        gr->GetXaxis()->SetTitleFont(font);
        gr->GetYaxis()->SetTitleFont(font);

        gr->GetXaxis()->SetLabelFont(font);
        gr->GetYaxis()->SetLabelFont(font);

        if (setXLimits) {
            gr->GetXaxis()->SetRangeUser(xMin, xMax);
        }

        if (setYLimits) {
            gr->GetYaxis()->SetRangeUser(yMin, yMax);
        }
    }

    static void FormatAxis(
        TAxis* xAxis,
        TAxis* yAxis,
        const TString& xTitle,
        const TString& yTitle,
        double titleSize = 0.03,
        double labelSize = 0.03,
        int titlefont = 62,
        bool setXLimits = false,
        double xMin = 0,
        double xMax = 0,
        bool setYLimits = false,
        double yMin = 0,
        double yMax = 0, int labelfont=42
    ) {
        xAxis->SetTitle(xTitle);
        yAxis->SetTitle(yTitle);

        xAxis->SetTitleSize(titleSize);
        yAxis->SetTitleSize(titleSize);

        xAxis->SetLabelSize(labelSize);
        yAxis->SetLabelSize(labelSize);

        xAxis->SetTitleFont(titlefont);
        yAxis->SetTitleFont(titlefont);
        
        xAxis->SetLabelFont(labelfont);
        yAxis->SetLabelFont(labelfont);

        if (setXLimits) xAxis->SetRangeUser(xMin, xMax);
        if (setYLimits) yAxis->SetRangeUser(yMin, yMax);
    }

    // ====== Crea Linea ======
    static TLine* CreateLine(
        double x1, double y1, double x2, double y2,
        int color = kRed,
        int style = 2,
        int width = 2
    ) {
        TLine* line = new TLine(x1, y1, x2, y2);
        line->SetLineColor(color);
        line->SetLineStyle(style);
        line->SetLineWidth(width);
        return line;
    }
    
    static void Format2DHisto(TH2* h2d,
                              int palette = kViridis,
                              const TString& drawOption = "COLZ",
                              double colorBarPos = 0.9,
                              double colorBarWidth = 0.03,
                              double colorBarHeight = 0.4,
                              bool showStats = false,
                              int statOptions = 1111,
                              double statPosX = 0.6,
                              double statPosY = 0.75,
                              double statWidth = 0.25,
                              double statHeight = 0.15) {

        gStyle->SetPalette(palette);
        h2d->SetStats(showStats);

        // Draw histogram
        h2d->Draw(drawOption);

        // Palette bar
        TPaletteAxis* paletteAxis = (TPaletteAxis*) h2d->GetListOfFunctions()->FindObject("palette");
        if (paletteAxis) {
            paletteAxis->SetX1NDC(colorBarPos);
            paletteAxis->SetX2NDC(colorBarPos + colorBarWidth);
            paletteAxis->SetY1NDC(0.5 - colorBarHeight / 2);
            paletteAxis->SetY2NDC(0.5 + colorBarHeight / 2);
            paletteAxis->Draw();
        }

        // Stat box options and position
        if (showStats) {
            TPaveStats* statBox = (TPaveStats*) h2d->FindObject("stats");
            if (statBox) {
                statBox->SetX1NDC(statPosX);
                statBox->SetX2NDC(statPosX + statWidth);
                statBox->SetY1NDC(statPosY);
                statBox->SetY2NDC(statPosY + statHeight);
                statBox->SetOptStat(statOptions);
                statBox->Draw();
            }
        }
    }
    
    static void SetupHistogramStyle(
        TH1* hist,
        bool showStats = false,
        int statOptions = 1111,            // default: entries, mean, RMS, under/overflow
        Color_t lineColor = kBlue,
        Color_t markerColor = kBlue,
        Style_t markerStyle = 20,
        int lineWidth = 1,
        const TString& drawOption = "",    // es. "" o "E"
        bool enableFillAlpha = false,      // abilita SetFillAlpha (default: off)
        Color_t fillColor = kBlue - 10,    // colore di riempimento (default: tonalità più chiara di blue)
        Style_t fillStyle = 1001,          // stile riempimento (default: pieno)
        double fillAlpha = 0.4,            // trasparenza, se abilitata (0=trasparente, 1=opaco)
        double statX1NDC = 0.6,            // posizione stat box (sinistra)
        double statX2NDC = 0.9,            // posizione stat box (destra)
        double statY1NDC = 0.6,            // posizione stat box (basso)
        double statY2NDC = 0.9             // posizione stat box (alto)
    ) {
        // Stile linee e marker
        hist->SetLineColor(lineColor);
        hist->SetMarkerColor(markerColor);
        hist->SetMarkerStyle(markerStyle);
        hist->SetLineWidth(lineWidth);

        // Imposta stile riempimento (fill)
        hist->SetFillColor(fillColor);
        hist->SetFillStyle(fillStyle);
        if (enableFillAlpha) {
            hist->SetFillColorAlpha(fillColor, fillAlpha);
        }

        // Stat box
        hist->SetStats(showStats);

        // Disegna istogramma
        hist->Draw(drawOption);

        if (showStats) {
            gPad->Update();  // Aggiorna il pad per far comparire la stat box
            TPaveStats* stats = (TPaveStats*)hist->FindObject("stats");
            if (stats) {
                stats->SetOptStat(statOptions);
                gStyle->SetOptFit(statOptions);  // Mostra parametri fit nel box

                // Imposta posizione stat box da parametri
                stats->SetX1NDC(statX1NDC);
                stats->SetX2NDC(statX2NDC);
                stats->SetY1NDC(statY1NDC);
                stats->SetY2NDC(statY2NDC);
                stats->Draw();
            }
        }
    }


    static void CustomizeGraph(TGraphErrors* graph,
                        int markerStyle,
                        int markerColor,
                        const char* drawOption,
                        int lineColor,
                        int lineStyle=1, double markerSize = 1)
    {
        if (!graph) return;

        graph->SetMarkerStyle(markerStyle);
        graph->SetMarkerColor(markerColor);
        graph->SetMarkerSize(markerSize);
        graph->SetLineColor(lineColor);
        graph->SetLineStyle(lineStyle);
        graph->SetLineWidth(2); // opzionale

        graph->Draw(drawOption);
    }

    static TLegend* CreateLegend(double x1, double y1, double x2, double y2,
                                 int fillColor = 0, double borderSize = 1, int font = 62, int textAlign = 22)
    {
        TLegend* leg = new TLegend(x1, y1, x2, y2);
        leg->SetBorderSize(borderSize);
        leg->SetFillColor(fillColor);
        leg->SetTextFont(font);
        leg->SetTextAlign(textAlign);
        
        return leg;
    }

    static void AddLegendEntry(TLegend* legend,
                               TObject* obj,
                               const TString& label,
                               const TString& option = "l") // default: sinistra, centrato verticalmente
    {
        if (!legend || !obj) return;

        // Crea l'entry con TLatex bold
        TLegendEntry* entry = legend->AddEntry(obj, Form("%s", label.Data()), option);
    }

    
    static TPaveText* CreateInfoBox(double x1, double y1, double x2, double y2, const TString& text,
                             Color_t fillColor = kWhite, Color_t textColor = kBlack,
                             double textSize = 0.03, int borderSize = 1)
    {
        TPaveText* pave = new TPaveText(x1, y1, x2, y2, "NDC"); // coordinate in Normalized Device Coordinates (0-1)
        pave->SetFillColor(fillColor);
        pave->SetTextColor(textColor);
        pave->SetTextSize(textSize);
        pave->SetBorderSize(borderSize);
        pave->SetShadowColor(0);

        // Separa la stringa multilinea in righe e le aggiunge una a una
        TObjArray* lines = text.Tokenize("\n");
        for (int i = 0; i < lines->GetEntriesFast(); ++i) {
            TString line = ((TObjString*)lines->At(i))->GetString();
            pave->AddText(line);
        }
        delete lines;

        return pave;
    }



};

#endif
