/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    writePCDict

Description
    Utility to generate a constant/phaseChange/phaseChangeDict dictionary file for
    phaseChangeFoam cases with full high-accuracy defaults and comments.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a constant/phaseChange/phaseChangeDict file for phaseChangeFoam simulations."
    );

    argList args(argc, argv);

    fileName casePath = args.path();
    fileName dictDir = casePath / "constant" / "phaseChange";
    fileName dictPath = dictDir / "phaseChangeDict";

    Info<< "Writing phaseChangeDict to " << dictPath << endl;

    mkDir(dictDir);

    OFstream os(dictPath);

    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << dictPath << " for writing"
            << exit(FatalError);
    }

    os  << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        << "| =========                 |                                                 |\n"
        << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        << "|  \\\\    /   O peration     | Version:  v2412                                 |\n"
        << "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
        << "|    \\\\/     M anipulation  |                                                 |\n"
        << "\\*---------------------------------------------------------------------------*/\n"
        << "FoamFile\n"
        << "{\n"
        << "    version     2.0;\n"
        << "    format      ascii;\n"
        << "    class       dictionary;\n"
        << "    location    \"constant/phaseChange\";\n"
        << "    object      phaseChangeDict;\n"
        << "}\n"
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
        << "active          true;\n\n"
        << "phaseChange\n"
        << "{\n"
        << "    active          true;\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 1. Solid Phase-Change Model & Direction Controls\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    phaseChangeMode EHCDarcy;   // Options: EHCDarcy (EHC + implicit Darcy drag), EHC (energy-only), enthalpyPorosity (alias), none\n"
        << "    direction       both;       // Options: both, forward (melting only), reverse (freezing only)\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 2. Forward / Melting Phase-Change Properties\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    forward\n"
        << "    {\n"
        << "        T_lowerBound    300.0;  // Solidus temperature T_solidus [K]\n"
        << "        T_upperBound    310.0;  // Liquidus temperature T_liquidus [K]\n"
        << "        latentHeat      100000.0; // Latent heat of fusion L [J/kg]\n"
        << "    }\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 3. Reverse / Freezing Phase-Change Properties (Thermal Hysteresis Mode)\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    reverse\n"
        << "    {\n"
        << "        T_lowerBound    295.0;  // Freezing lower bound temperature [K]\n"
        << "        T_upperBound    305.0;  // Freezing upper bound temperature [K]\n"
        << "        latentHeat      100000.0; // Latent heat of freezing L [J/kg]\n"
        << "    }\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 4. Stateful Thermal Hysteresis Tracking\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    hysteresis\n"
        << "    {\n"
        << "        active              true;  // Enable path-dependent thermal hysteresis\n"
        << "        reversalTolerance   1e-6;  // Temperature reversal detection tolerance [K]\n"
        << "    }\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 5. Phase Density & Thermophysical Properties\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    density\n"
        << "    {\n"
        << "        model                       linear;   // Options: linear, thermo\n"
        << "        rhoRef                      1000.0;   // Reference density [kg/m^3]\n"
        << "        rhoSolid                    1000.0;   // Solid-phase density [kg/m^3]\n"
        << "        rhoLiquid                   1000.0;   // Liquid-phase density [kg/m^3]\n"
        << "        allowNonConservativeDensity true;     // Allow unequal solid/liquid density\n"
        << "    }\n\n"
        << "    thermophysical\n"
        << "    {\n"
        << "        mode            custom; // Options: custom, thermo\n"
        << "        CpSolid         1980.0; // Solid specific heat capacity Cp [J/(kg K)]\n"
        << "        CpLiquid        2320.0; // Liquid specific heat capacity Cp [J/(kg K)]\n"
        << "        kSolid          1.0;    // Solid thermal conductivity k [W/(m K)]\n"
        << "        kLiquid         1.0;    // Liquid thermal conductivity k [W/(m K)]\n"
        << "    }\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 6. Enthalpy-Porosity (Mushy Zone Flow Resistance)\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    porosity\n"
        << "    {\n"
        << "        A_cu            1.0e5;  // Mushy zone Darcy constant [kg/(m^3 s)]\n"
        << "        eps             1.0e-3; // Small division tolerance to prevent zero-divide\n"
        << "    }\n\n\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    // 7. Fluid Phase-Change Model (Evaporation / Condensation)\n"
        << "    // -------------------------------------------------------------------------\n"
        << "    type            Lee;        // Options: Lee, constantSource, none\n"
        << "    liquid          liquid;     // Liquid species name in thermo composition\n"
        << "    vapor           vapor;      // Vapor species name in thermo composition\n"
        << "    C_evap          0.1;        // Evaporation rate frequency coefficient [1/s]\n"
        << "    C_cond          0.1;        // Condensation rate frequency coefficient [1/s]\n"
        << "    latentHeat      2.26e6;     // Latent heat of vaporization [J/kg]\n"
        << "    Tsat            373.15;     // Saturation temperature [K]\n"
        << "    pRef            101325;     // Reference pressure for saturation [Pa]\n"
        << "    enableTsatP     false;      // Enable Clausius-Clapeyron pressure-dependent Tsat(p)\n"
        << "}\n\n"
        << "// ************************************************************************* //\n";

    Info<< "Successfully generated " << dictPath << endl;

    return 0;
}

// ************************************************************************* //
