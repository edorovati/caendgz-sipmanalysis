/******************************************************************************
 * @file ColorPaletteManager.h
 * @brief Declaration of the ColorPaletteManager class, which provides a set of
 *        static utility functions to retrieve colors based on specific overvoltage
 *        values. Each color is represented by its RGB value and associated with
 *        an overvoltage threshold.
 *
 * This header defines:
 *   - ColorPaletteManager: routines for obtaining colors corresponding to
 *            different overvoltage values.
 *
 ******************************************************************************/
#ifndef COLORPALETTEMANAGER_H
#define COLORPALETTEMANAGER_H

#include <TColor.h>
#include <vector>
#include <string>

/******************************************************************************
 *                                ColorPaletteManager Class
 * ----------------------------------------------------------------------------
 * @class ColorPaletteManager
 * @brief Collection of static functions for retrieving colors based on specific
 *        overvoltage values. These colors are pre-defined and associated with
 *        values ranging from 0.0 to 10.0, providing a visual representation
 *        of the overvoltage levels.
 *
 * Public Functions:
 * ----------------------------------------------------------------------------
 * | Name                          | Description                                           |
 * ----------------------------------------------------------------------------
 * | getColorFromOvervoltage        | Retrieve a color based on the overvoltage value.      |
 * ----------------------------------------------------------------------------
 ******************************************************************************/

class ColorPaletteManager {
public:
// ------------------------------------------------------------------------------------ //
/// @brief  Retrieve the color associated with a specific overvoltage value.
/// @param  overvoltage  The overvoltage value (in volts).
/// @return The corresponding color as an integer value, representing RGB color.
///         Returns a fallback color (black) if the overvoltage is not recognized.
// ------------------------------------------------------------------------------------ //
    static Int_t getColorFromOvervoltage(double overvoltage) {
        if (std::abs(overvoltage - 0.0) < 1e-3) return customGoldenrod();
        if (std::abs(overvoltage - 0.5) < 1e-3) return customCrimson();
        if (std::abs(overvoltage - 1.0) < 1e-3) return customForestGreen();
        if (std::abs(overvoltage - 1.5) < 1e-3) return customRoyalBlue();
        if (std::abs(overvoltage - 2.0) < 1e-3) return customOrchid();
        if (std::abs(overvoltage - 2.5) < 1e-3) return customCoral();
        if (std::abs(overvoltage - 3.0) < 1e-3) return customTurquoise();
        if (std::abs(overvoltage - 3.5) < 1e-3) return customSlateGray();
        if (std::abs(overvoltage - 4.0) < 1e-3) return customDeepPink();
        if (std::abs(overvoltage - 4.5) < 1e-3) return customChocolate();
        if (std::abs(overvoltage - 5.0) < 1e-3) return customDodgerBlue();
        if (std::abs(overvoltage - 5.5) < 1e-3) return customMediumVioletRed();
        if (std::abs(overvoltage - 6.0) < 1e-3) return customLimeGreen();
        if (std::abs(overvoltage - 6.5) < 1e-3) return customDarkOrange();
        if (std::abs(overvoltage - 7.0) < 1e-3) return customMediumSlateBlue();
        if (std::abs(overvoltage - 7.5) < 1e-3) return customDarkCyan();
        if (std::abs(overvoltage - 8.0) < 1e-3) return customHotPink();
        if (std::abs(overvoltage - 8.5) < 1e-3) return customOliveDrab();
        if (std::abs(overvoltage - 9.0) < 1e-3) return customSteelBlue();
        if (std::abs(overvoltage - 9.5) < 1e-3) return customSienna();
        if (std::abs(overvoltage - 10.0) < 1e-3) return customPlum();

        std::cerr << "[WARN] Overvoltage " << overvoltage << " unrecognised. Use fallback colour (black)." << std::endl;
        return kBlack;
    }
    
private:
    static Int_t customGoldenrod()      { return TColor::GetColor(218, 165, 32); }   // intense gold
    static Int_t customCrimson()        { return TColor::GetColor(220, 20, 60); }    // dark red
    static Int_t customForestGreen()    { return TColor::GetColor(34, 139, 34); }    // forest green
    static Int_t customRoyalBlue()      { return TColor::GetColor(65, 105, 225); }   // royal blue
    static Int_t customOrchid()         { return TColor::GetColor(218, 112, 214); }  // soft violet
    static Int_t customCoral()          { return TColor::GetColor(255, 127, 80); }   // orange/pink
    static Int_t customTurquoise()      { return TColor::GetColor(64, 224, 208); }   // turquoise
    static Int_t customSlateGray()      { return TColor::GetColor(112, 128, 144); }  // blue-grey
    static Int_t customDeepPink()       { return TColor::GetColor(255, 20, 147); }   // neon pink
    static Int_t customChocolate()      { return TColor::GetColor(210, 105, 30); }   // warm brown
    static Int_t customDodgerBlue()     { return TColor::GetColor(30, 144, 255); }   // bright blue
    static Int_t customMediumVioletRed(){ return TColor::GetColor(199, 21, 133); }   // violet red
    static Int_t customLimeGreen()      { return TColor::GetColor(50, 205, 50); }    // lime green
    static Int_t customDarkOrange()     { return TColor::GetColor(255, 140, 0); }    // dark orange
    static Int_t customMediumSlateBlue(){ return TColor::GetColor(123, 104, 238); }  // lavender blue
    static Int_t customDarkCyan()       { return TColor::GetColor(0, 139, 139); }    // dark cyan
    static Int_t customHotPink()        { return TColor::GetColor(255, 105, 180); }  // hot pink
    static Int_t customOliveDrab()      { return TColor::GetColor(107, 142, 35); }   // olive drab green
    static Int_t customSteelBlue()      { return TColor::GetColor(70, 130, 180); }   // steel blue
    static Int_t customSienna()         { return TColor::GetColor(160, 82, 45); }    // sienna brown
    static Int_t customPlum()           { return TColor::GetColor(221, 160, 221); }  // light plum

};

#endif  // COLORPALETTEMANAGER_H
