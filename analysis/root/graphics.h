/******************************************************************************
 * @file GraphicsUtils.h
 * @brief Declaration of the GraphicsUtils class, offering a set of static
 *        utility functions for styling and drawing ROOT graphs, canvases,
 *        legends, and confidence bands.
 *
 * This header defines:
 *   - GraphicsUtils: convenience routines to apply styles, create canvases,
 *     set axis labels, build band graphs, and manage legends and text boxes.
 *
 ******************************************************************************/

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <vector>


/******************************************************************************
 *                              GraphicsUtils Class
 * ----------------------------------------------------------------------------
 * @class GraphicsUtils
 * @brief Provides a suite of static methods to simplify creation and styling
 *        of ROOT graphical objects.
 *
 * Public Functions:
 * ----------------------------------------------------------------------------
 * | Name                  | Description                                            |
 * ----------------------------------------------------------------------------
 * | setGraphStyle         | Apply marker, line, and color style to a TGraphErrors |
 * ----------------------------------------------------------------------------
 * | createCanvas          | Instantiate and configure a TCanvas (size, grid, log) |
 * ----------------------------------------------------------------------------
 * | setAxisAndDraw        | Set axis titles and draw a TMultiGraph with options    |
 * ----------------------------------------------------------------------------
 * | createLegend          | Create a TLegend with optional transparent background |
 * ----------------------------------------------------------------------------
 * | createBandGraph       | Build a filled band (TGraph) from top/bottom y-values |
 * ----------------------------------------------------------------------------
 * | drawGraphWithBand     | Draw TGraphErrors overlaid with a confidence band     |
 * ----------------------------------------------------------------------------
 * | addEntryToLegend      | Add a graph entry to a TLegend with specified style    |
 * ----------------------------------------------------------------------------
 * | createPaveText        | Create a TPaveText box populated with multiple lines   |
 * ----------------------------------------------------------------------------
 ******************************************************************************/
struct BandPoint {
    double x;
    double yTop;
    double yBot;
};

class GraphicsUtils {
public:
// ------------------------------------------------------------------------------------ //
/// @brief  Apply marker, line and color style to a TGraphErrors.
/// @param  graph        Pointer to the TGraphErrors to style.
/// @param  markerStyle  ROOT marker style code.
/// @param  color        Color index for both marker and line.
/// @param  lineStyle    Line style code (default = 1).
// ------------------------------------------------------------------------------------ //
    static void setGraphStyle(TGraphErrors* graph, int markerStyle, int color, int lineStyle = 1) {
        graph->SetMarkerStyle(markerStyle);
        graph->SetLineColor(color);
        graph->SetMarkerColor(color);
        graph->SetLineStyle(lineStyle);
    }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Create and configure a new TCanvas.
/// @param  name   Internal name of the canvas.
/// @param  title  Display title of the canvas.
/// @param  width  Canvas width in pixels.
/// @param  height Canvas height in pixels.
/// @param  logY   If true, enable logarithmic y-axis.
/// @return Pointer to the newly created TCanvas.
// ------------------------------------------------------------------------------------ //
    static TCanvas* createCanvas(const char* name, const char* title, int width, int height, bool logY = false) {
        TCanvas* canvas = new TCanvas(name, title, width, height);
        if (logY) {
            canvas->SetLogy();
        }
        canvas->SetGrid();
        return canvas;
    }
// ------------------------------------------------------------------------------------ //
    
// ------------------------------------------------------------------------------------//
/// @brief  Set axis titles and draw a TMultiGraph.
/// @param  mg       Pointer to the TMultiGraph to draw.
/// @param  xTitle   Label for the x-axis.
/// @param  yTitle   Label for the y-axis.
/// @param  drawAPE  If true, draw with "APE" option; otherwise "ALE".
// ------------------------------------------------------------------------------------//
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
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Create a legend box.
/// @param  x1     Lower-left x (NDC).
/// @param  y1     Lower-left y (NDC).
/// @param  x2     Upper-right x (NDC).
/// @param  y2     Upper-right y (NDC).
/// @param  noFill If true, disable fill, border, and lines.
/// @return Pointer to the newly created TLegend.
// ------------------------------------------------------------------------------------ //
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
// ------------------------------------------------------------------------------------ //
    
    
// ------------------------------------------------------------------------------------ //
/// @brief  Build a filled “band” graph from top/bottom y-values.
/// @param  bandData  Vector of BandPoint {x, yTop, yBot}.
/// @return Pointer to the TGraph representing the band.
// ------------------------------------------------------------------------------------ //
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
// ------------------------------------------------------------------------------------ //
    
// ------------------------------------------------------------------------------------ //
/// @brief  Draw a TGraphErrors with an overlaid band on its own canvas.
/// @param  graph       Pointer to the TGraphErrors to draw.
/// @param  bandData    Vector of BandPoint defining the band.
/// @param  title       Canvas title.
/// @param  xTitle      Label for the x-axis.
/// @param  yTitle      Label for the y-axis.
/// @param  canvasName  Internal name of the canvas.
// ------------------------------------------------------------------------------------ //
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
// ------------------------------------------------------------------------------------ //
    
// ------------------------------------------------------------------------------------ //
/// @brief  Add an entry to a legend for a TGraph.
/// @param  legend  Pointer to the TLegend.
/// @param  graph   Pointer to the TGraph to add.
/// @param  label   Text label for the legend entry.
/// @param  type    Entry type ("lep" or "f").
// ------------------------------------------------------------------------------------ //
    static void addEntryToLegend(TLegend* legend, TGraph* graph, const char* label, const char* type) {
            if (strcmp(type, "lep") == 0) {
                legend->AddEntry(graph, label, "lep");
            } else if (strcmp(type, "f") == 0) {
                legend->AddEntry(graph, label, "f");
            } else {
                std::cerr << "Tipo di entry non valido: " << type << std::endl;
            }
        }
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
/// @brief  Create a TPaveText box with multiple lines.
/// @param  x1     Lower-left x (NDC).
/// @param  y1     Lower-left y (NDC).
/// @param  x2     Upper-right x (NDC).
/// @param  y2     Upper-right y (NDC).
/// @param  lines  Vector of text lines to add.
/// @return Pointer to the newly created TPaveText.
// ------------------------------------------------------------------------------------ //
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
// ------------------------------------------------------------------------------------ //
    
};

#endif 
