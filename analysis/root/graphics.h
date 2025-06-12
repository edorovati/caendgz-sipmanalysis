/******************************************************************************
 * @file GraphicsUtils.h
 * @brief Declaration of the GraphicsUtils class: static utilities for styling
 *        and drawing ROOT canvases, graphs, legends, and confidence bands.
 ******************************************************************************/

#ifndef GRAPHICSUTILS_H
#define GRAPHICSUTILS_H

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <vector>
#include <string>

/**
 * @struct BandPoint
 * @brief Data point for constructing confidence bands.
 *
 * Contains an x-coordinate and upper/lower y-values.
 */
struct BandPoint {
    double x;     ///< x-coordinate
    double yTop;  ///< upper band boundary at x
    double yBot;  ///< lower band boundary at x
};

/******************************************************************************
 *                              GraphicsUtils Class
 * ----------------------------------------------------------------------------
 * @class GraphicsUtils
 * @brief Provides a suite of static methods to simplify creation and styling
"
 *        of ROOT graphical objects, including canvases, graphs, legends, and bands.
"
 *
 * This class includes the following public static functions:
 * ----------------------------------------------------------------------------
 * | **Public Functions**            | **Description**                                       |
 * ----------------------------------------------------------------------------
 * | **createCanvas**                | Instantiate and configure a TCanvas                   |
 * | **setGraphStyle**               | Apply marker, line, and color style to TGraphErrors   |
 * | **setAxesBoldTitles**           | Set bold axis titles and optional axis limits         |
 * | **createBandGraph**             | Build a vertical confidence band (filled)             |
 * | **createHorizontalBandGraph**   | Build a horizontal confidence band                    |
 * | **drawGraphWithBand**           | Draw TGraphErrors overlaid with a confidence band     |
 * | **createDummyGraph**            | Create a dummy graph to fix axes and ranges           |
 * | **createLegend**                | Create a TLegend with optional transparent styling    |
 * | **addEntryToLegend**            | Add a graph entry to a legend, handling micro-units   |
 * | **addLabelToCanvas**            | Place a TPaveText annotation on a canvas              |
 * | **setAxisAndDraw**              | Set axis titles on TMultiGraph and draw it            |
 * | **createPaveText**              | Create a TPaveText box with multiple text lines       |
 ******************************************************************************/
class GraphicsUtils {
public:
    /**
     * @brief Create and configure a new TCanvas.
     * @param name    Internal canvas name.
     * @param title   Display title.
     * @param width   Width in pixels.
     * @param height  Height in pixels.
     * @param logY    If true, enable logarithmic Y-axis.
     * @return Pointer to created TCanvas.
     */
    static TCanvas* createCanvas(const char* name, const char* title,
                                 int width, int height, bool logY = false) {
        TCanvas* canvas = new TCanvas(name, title, width, height);
        if (logY) canvas->SetLogy();
        canvas->SetGrid();
        return canvas;
    }

    /**
     * @brief Apply marker style, line style, and color to a TGraphErrors.
     * @param graph        Graph to style.
     * @param markerStyle  ROOT marker style code.
     * @param color        Color index for both marker and line.
     * @param lineStyle    Line style code (default 1).
     */
    static void setGraphStyle(TGraphErrors* graph, int markerStyle,
                              int color, int lineStyle = 1) {
        graph->SetMarkerStyle(markerStyle);
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetLineStyle(lineStyle);
    }

    /**
     * @brief Set bold axis titles and optional axis limits on a graph.
     * @param graph     Graph to adjust.
     * @param xTitle    X-axis title.
     * @param yTitle    Y-axis title.
     * @param setLimits If true, apply the provided limits.
     * @param xMin      Minimum X (if setLimits).
     * @param xMax      Maximum X (if setLimits).
     * @param yMin      Minimum Y (if setLimits).
     * @param yMax      Maximum Y (if setLimits).
     */
    static void setAxesBoldTitles(TGraph* graph,
                                  const TString& xTitle,
                                  const TString& yTitle,
                                  bool setLimits = false,
                                  double xMin = 0.0,
                                  double xMax = 0.0,
                                  double yMin = 0.0,
                                  double yMax = 0.0) {
        graph->GetXaxis()->SetTitle(xTitle);
        graph->GetXaxis()->SetTitleFont(62);
        graph->GetXaxis()->SetTitleSize(0.03);

        graph->GetYaxis()->SetTitle(yTitle);
        graph->GetYaxis()->SetTitleFont(62);
        graph->GetYaxis()->SetTitleSize(0.03);

        if (setLimits) {
            graph->GetXaxis()->SetLimits(xMin, xMax);
            graph->SetMinimum(yMin);
            graph->SetMaximum(yMax);
        }
    }

    /**
     * @brief Create a filled band graph from top/bottom values.
     * @param bandData   Vector of BandPoint {x,yTop,yBot}.
     * @param fillColor  Fill color index.
     * @param fillAlpha  Fill transparency.
     * @param lineColor  Outline color index.
     * @param lineStyle  Outline line style.
     * @return Pointer to the TGraph band.
     */
    static TGraph* createBandGraph(const std::vector<BandPoint>& bandData,
                                   Color_t fillColor = kBlue,
                                   Double_t fillAlpha = 0.2,
                                   Color_t lineColor = kBlue,
                                   Style_t lineStyle = 2) {
        int n = bandData.size();
        int tot = 2*n + 1;
        TGraph* g = new TGraph(tot);
        // upper
        for (int i=0; i<n; ++i) g->SetPoint(i, bandData[i].x, bandData[i].yTop);
        // lower
        for (int i=0; i<n; ++i) g->SetPoint(n+i,
            bandData[n-1-i].x, bandData[n-1-i].yBot);
        // close
        g->SetPoint(tot-1, bandData[0].x, bandData[0].yTop);
        g->SetFillColorAlpha(fillColor, fillAlpha);
        g->SetLineColor(lineColor);
        g->SetLineStyle(lineStyle);
        return g;
    }

    /**
     * @brief Create a horizontal band graph from swapped bandpoints.
     * @param bandData   Vector of BandPoint {x,yTop,yBot}.
     * @param fillColor  Fill color.
     * @param fillAlpha  Transparency.
     * @param lineColor  Outline color.
     * @param lineStyle  Outline style.
     * @return Pointer to TGraph.
     */
    static TGraph* createHorizontalBandGraph(const std::vector<BandPoint>& bandData,
                                             Color_t fillColor = kBlue,
                                             Double_t fillAlpha = 0.2,
                                             Color_t lineColor = kBlue,
                                             Style_t lineStyle = 2) {
        int n = bandData.size();
        int tot = 2*n + 1;
        TGraph* g = new TGraph(tot);
        for (int i=0; i<n; ++i) g->SetPoint(i, bandData[i].yTop, bandData[i].x);
        for (int i=0; i<n; ++i) g->SetPoint(n+i,
            bandData[n-1-i].yBot, bandData[n-1-i].x);
        g->SetPoint(tot-1, bandData[0].yTop, bandData[0].x);
        g->SetFillColorAlpha(fillColor, fillAlpha);
        g->SetLineColor(lineColor);
        g->SetLineStyle(lineStyle);
        return g;
    }

    /**
     * @brief Draw a TGraphErrors with overlaid band in a new canvas.
     * @param graph       Data graph (with errors).
     * @param bandData    Band points for confidence region.
     * @param title       Canvas title.
     * @param xTitle      X-axis label.
     * @param yTitle      Y-axis label.
     * @param canvasName  Internal name for canvas.
     */
    static void drawGraphWithBand(TGraphErrors* graph,
                                  const std::vector<BandPoint>& bandData,
                                  const char* title,
                                  const char* xTitle,
                                  const char* yTitle,
                                  const char* canvasName) {
        TCanvas* c = createCanvas(canvasName, title, 800, 600);
        TGraph* band = createBandGraph(bandData);
        TMultiGraph* mg = new TMultiGraph();
        mg->Add(band, "F");
        mg->Add(graph, "PE");
        setAxisAndDraw(mg, xTitle, yTitle);
        c->Draw();
    }

    /**
     * @brief Create a legend box.
     * @param x1     Lower-left x (NDC).
     * @param y1     Lower-left y (NDC).
     * @param x2     Upper-right x (NDC).
     * @param y2     Upper-right y (NDC).
     * @param noFill If true, disable fill and border.
     * @return Pointer to TLegend.
     */
    static TLegend* createLegend(double x1, double y1,
                                 double x2, double y2,
                                 bool noFill = true) {
        TLegend* legend = new TLegend(x1,y1,x2,y2);
        if (noFill) {
            legend->SetFillStyle(0);
            legend->SetLineStyle(0);
            legend->SetBorderSize(0);
        }
        legend->SetTextFont(62);
        legend->SetTextSize(0.03);
        return legend;
    }

    /**
     * @brief Add a graph entry to a legend, handling micro-unit replacement.
     * @param legend  Target legend.
     * @param graph   Graph pointer.
     * @param label   Entry text ("um"→"µm").
     * @param option  Draw option ("lep" or "f").
     */
    static void addEntryToLegend(TLegend* legend,
                                 TGraph* graph,
                                 const char* label,
                                 const char* option) {
        TString txt(label);
        if (txt.Contains("um")) txt.ReplaceAll("um","#mum");
        legend->AddEntry(graph, txt, option);
    }

    /**
     * @brief Add a TPaveText label on a canvas.
     * @param canvas     Canvas to modify.
     * @param text       Raw label text ("um"→"µm").
     * @param xmin,ymin  NDC lower-left.
     * @param xmax,ymax  NDC upper-right.
     * @param fillColor  Background fill.
     * @param textAlign  Text alignment.
     * @param borderSize Border size.
     * @param textFont   Font code.
     */
    static void addLabelToCanvas(TCanvas* canvas,
                                 const TString& text,
                                 double xmin=0.15, double ymin=0.82,
                                 double xmax=0.45, double ymax=0.88,
                                 int fillColor=0,
                                 int textAlign=12,
                                 int borderSize=1,
                                 int textFont=42) {
        canvas->cd();
        TString lbl(text);
        if (lbl.Contains("um")) lbl.ReplaceAll("um","#mum");
        TPaveText* pt = new TPaveText(xmin,ymin,xmax,ymax,"NDC");
        pt->SetFillColor(fillColor);
        pt->SetTextAlign(textAlign);
        pt->SetBorderSize(borderSize);
        pt->SetTextFont(textFont);
        pt->AddText(lbl);
        pt->Draw();
    }

    /**
     * @brief Draw a TMultiGraph with bold axis titles.
     * @param mg        Pointer to TMultiGraph.
     * @param xTitle    X-axis title.
     * @param yTitle    Y-axis title.
     * @param drawAPE   If true use "APE", else "ALE".
     */
    static void setAxisAndDraw(TMultiGraph* mg,
                               const char* xTitle,
                               const char* yTitle,
                               bool drawAPE=true) {
        mg->SetTitle(Form(";%s;%s",xTitle,yTitle));
        mg->GetXaxis()->SetTitleFont(62);
        mg->GetYaxis()->SetTitleFont(62);
        mg->Draw(drawAPE?"APE":"ALE");
    }

    /**
     * @brief Create a TPaveText box containing multiple lines of text.
     * @param x1,y1   NDC lower-left.
     * @param x2,y2   NDC upper-right.
     * @param lines   Vector of text lines.
     * @return TPaveText pointer.
     */
    static TPaveText* createPaveText(double x1, double y1,
                                     double x2, double y2,
                                     const std::vector<std::string>& lines) {
        TPaveText* pave = new TPaveText(x1,y1,x2,y2,"NDC");
        pave->SetFillColor(0);
        pave->SetTextFont(42);
        pave->SetTextAlign(12);
        for (auto& l : lines) pave->AddText(l.c_str());
        return pave;
    }
};

#endif // GRAPHICSUTILS_H
