ROOT-Based Analysis Utilities
====================================
A suite of utility classes designed to support signal, statistical,
and graphical analysis for waveform-based lab activities.


Classes Overview
====================================

1. Analysis        - Signal fitting, threshold counting, CT/DCR ratio, dark rate.
2. GraphicsUtils   - Styling and drawing utilities for ROOT graphs.
3. Utils           - File I/O, interpolation, graph creation, derivatives, etc.
4. ColorPaletteManager - Color mappings based on several conditions.

Analysis Class (analysis.h)
====================================
Static methods for statistical operations on waveform or histogram data.

▸ `fitGaussian(TGraph*, fitMin, fitMax, amp, mean, sigma)`
   → Fits a Gaussian to a graph. Returns (mean, sigma, meanErr, sigmaErr, χ²/ndf)

▸ `calculateCounts(x, y, threshold)`
   → Interpolates number of counts at given threshold. Returns (counts, error)

▸ `calculateCTDCRRatio(DCR, DCRerr, CT, CTerr)`
   → Returns cross-talk to DCR ratio (%) and associated propagated uncertainty

▸ `calculateDarkRateAndError(DCR, DCRerr, wf, timeWindow)`
   → Calculates the dark count rate (kHz) from raw DCR and time window


GraphicsUtils Class (graphics.h)
====================================
ROOT graph utilities for styling, drawing, and canvas/legend configuration.

▸ `setGraphStyle(graph, markerStyle, color, lineStyle=1)`
   → Set marker/line style and color on a TGraphErrors

▸ `createCanvas(name, title, width, height, logY=false)`
   → Creates a TCanvas with optional log scale

▸ `setAxisAndDraw(mg, xTitle, yTitle, drawAPE=true)`
   → Sets axis titles and draws a TMultiGraph ("APE" or "ALE")

▸ `createLegend(x1, y1, x2, y2, noFill=true)`
   → Configures a legend (optionally transparent)

▸ `createBandGraph(vector<BandPoint>)`
   → Builds a confidence band (TGraph) from yTop/yBot vectors

▸ `drawGraphWithBand(graph, bandData, title, xTitle, yTitle, canvasName)`
   → Plots a TGraphErrors overlaid with a filled band

▸ `addEntryToLegend(legend, graph, label, type)`
   → Adds a graph to a legend with marker ("lep") or filled ("f")

▸ `createPaveText(x1, y1, x2, y2, lines)`
   → Creates a TPaveText with multiple text lines


Utils Class (utils.h)
====================================
A set of helper functions for numerical processing and ROOT I/O.

▸ `extractOvervoltage(filename)`
   → Extracts a numeric overvoltage value (supports "5,0OV.txt" format)

▸ `interpolate(x, x1, y1, x2, y2, yErr1, yErr2)`
   → Linear interpolation with error propagation between two points

▸ `readDataFile(filepath, xVec, yVec)`
   → Loads column data from file (thresholds, counts, etc.)

▸ `createGraph(xVals, yVals, xErrs, yErrs)`
   → Builds a TGraphErrors from data and error vectors

▸ `findMaximum(x, y, threshold)`
   → Returns the index of the maximum y-value beyond a threshold

▸ `findLocalMaxima(y)`
   → Finds the indices of local maxima in a vector of data

▸ `computeDerivativeWithErrors(x, y, yErr, correlation)`
   → Computes numerical derivative of y vs x, with error propagation

▸ `writeGraphTTree(graph, name)`
   → Writes a TGraphErrors into the current ROOT file


ColorPaletteManager Class (colorpalettemanager.h)
================================================
A utility class providing color mappings based on different cases

▸ `getColorFromOvervoltage(overvoltage)`
   → Returns a color corresponding to the given overvoltage value. 

