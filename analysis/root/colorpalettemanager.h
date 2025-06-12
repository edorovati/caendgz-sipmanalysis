/******************************************************************************
 * @file ColorPaletteManager.h
 * @brief Declaration of the ColorPaletteManager class: static utilities to map
 *        overvoltage values to ROOT colors (RGB).
 *
 * This header provides:
 *   - ColorPaletteManager::getColorFromOvervoltage : Retrieve ROOT color by ov.
 *   - ColorPaletteManager::getColorByIndex         : Retrieve color cycling indices
 ******************************************************************************/
#ifndef COLORPALETTEMANAGER_H
#define COLORPALETTEMANAGER_H

#include <TColor.h>
#include <vector>
#include <string>
#include <map>

/******************************************************************************
 *                          ColorPaletteManager Class
 * ----------------------------------------------------------------------------
 * @class ColorPaletteManager
 * @brief Collection of static methods for assigning colors based on overvoltage
 *        levels or numeric indices, facilitating consistent graph styling.
 *
 * This class includes the following public static functions:
 * ----------------------------------------------------------------------------
 * | **Public Functions**             | **Description**                           |
 * ----------------------------------------------------------------------------
 * | **getColorFromOvervoltage**      | Map specific OV values to predefined RGB  |
 * ----------------------------------------------------------------------------
 * | **getColorByIndex**              | Cycle through a palette by integer index  |
 * ----------------------------------------------------------------------------
 ******************************************************************************/
class ColorPaletteManager {
public:
    //--------------------------------------------------------------------------------
    /// @brief Retrieve ROOT color for a given overvoltage (in volts).
    /// @param overvoltage Overvoltage value to map.
    /// @return ROOT color index (fallback black if unrecognized).
    //--------------------------------------------------------------------------------
    static Int_t getColorFromOvervoltage(double overvoltage) {
        if (std::abs(overvoltage - 52) < 1e-3) return kOrange + 7;
        if (std::abs(overvoltage - 52.5) < 1e-3) return kRed + 1;
        if (std::abs(overvoltage - 53) < 1e-3) return kGreen + 2;
        if (std::abs(overvoltage - 53.5) < 1e-3) return kBlue + 1;
        if (std::abs(overvoltage - 54) < 1e-3) return kViolet + 6;
        if (std::abs(overvoltage - 54.5) < 1e-3) return kPink + 7;
        if (std::abs(overvoltage - 55) < 1e-3) return kCyan - 6;
        if (std::abs(overvoltage - 56) < 1e-3) return kGray + 2;
        if (std::abs(overvoltage - 57) < 1e-3) return kMagenta - 3;
        if (std::abs(overvoltage - 58) < 1e-3) return kOrange - 3;
        if (std::abs(overvoltage - 5.0) < 1e-3) return kAzure + 5;
        if (std::abs(overvoltage - 5.5) < 1e-3) return kSpring + 9;
        if (std::abs(overvoltage - 6.0) < 1e-3) return kGreen + 1;
        if (std::abs(overvoltage - 6.5) < 1e-3) return kOrange + 3;
        if (std::abs(overvoltage - 7.0) < 1e-3) return kBlue + 3;
        if (std::abs(overvoltage - 7.5) < 1e-3) return kTeal + 4;
        if (std::abs(overvoltage - 8.0) < 1e-3) return kPink + 5;
        if (std::abs(overvoltage - 8.5) < 1e-3) return kSpring + 4;
        if (std::abs(overvoltage - 9.0) < 1e-3) return kAzure + 3;
        if (std::abs(overvoltage - 9.5) < 1e-3) return kOrange + 1;
        if (std::abs(overvoltage - 10.0) < 1e-3) return kViolet + 2;
        std::cerr << "[WARN] OV " << overvoltage << " unrecognized; fallback black." << std::endl;
        return kBlack;
    }

    //--------------------------------------------------------------------------------
    /// @brief Get a reusable ROOT color cycling through an internal palette.
    /// @param index Zero-based index to cycle.
    /// @return Newly allocated ROOT color index corresponding to RGB tuple.
    //--------------------------------------------------------------------------------
    static Int_t getColorByIndex(size_t index) {
        static const std::map<size_t, std::tuple<Int_t,Int_t,Int_t>> colorMap = {
            {0,{255,165,0}}, {1,{255,0,0}},   {2,{0,255,0}},  {3,{0,0,255}},
            {4,{238,130,238}},{5,{255,105,180}},{6,{0,255,255}},{7,{169,169,169}},
            {8,{255,0,255}},{9,{255,140,0}},{10,{240,255,255}},{11,{255,255,0}},
            {12,{0,128,128}},{13,{240,255,255}}
        };
        size_t idx = index % colorMap.size();
        auto [r,g,b] = colorMap.at(idx);
        Int_t col = TColor::GetFreeColorIndex();
        new TColor(col, r/255.0, g/255.0, b/255.0);
        return col;
    }
};

#endif // COLORPALETTEMANAGER_H
