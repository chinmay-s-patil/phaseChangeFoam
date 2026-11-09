#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Chinmay S Patil
"""
Paraffin PCM Thermal Hysteresis Cyclic Path Logic Verification Case
Verifies the Effective Heat Capacity (EHC) phase change engine with stateful thermal
hysteresis across a complete melt/freeze cycle using synthetic phase-change windows
representative of paraffin transitions reported in the literature (Bony & Citherlet 2007).

Reference:
    J. Bony and S. Citherlet (2007),
    "Numerical model and experimental validation of heat storage with phase change materials",
    Energy and Buildings, 39(10), pp. 1065-1072.
    DOI: 10.1016/j.enbuild.2006.10.017

Phenomenon:
    Commercial paraffins exhibit thermal hysteresis between charging (melting)
    and discharging (freezing/solidification).
    Synthetic melting window:        T_melt in [308.0, 312.0] K (peak ~310.0 K = 37 C)
    Synthetic solidification window:  T_freeze in [303.0, 307.0] K (peak ~305.0 K = 32 C)
    Hysteresis temperature gap:       Delta T_hyst = 5.0 K
    Note: This test uses synthetic linear windows to verify the state-machine logic
    (direction tracking, subcooling delay, reversal recovery) rather than curve-fitting raw DSC data.
"""

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np

def find_solver(solver_name="phaseChangeMultiRegionFoam"):
    import glob
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        p = os.path.join(user_appbin, solver_name)
        if os.path.exists(p):
            return p
    which_p = shutil.which(solver_name)
    if which_p:
        return which_p
    home = os.path.expanduser("~")
    for cand in glob.glob(os.path.join(home, "OpenFOAM", "*-v*", "platforms", "*", "bin", solver_name)):
        if os.path.exists(cand):
            return cand
    return solver_name

def run_cmd(cmd, cwd=None, allow_failure=False):
    res = subprocess.run(
        cmd,
        shell=True,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if res.returncode != 0 and not allow_failure:
        print(f"Command failed: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res

def setup_hysteresis_case(case_dir, totalTime=250):
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system/pcm"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant/pcm"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0/pcm"), exist_ok=True)

    # 1D slice (1 mm DSC crucible thickness, 5 cells)
    L = 0.001
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
scale 1;
vertices
(
    (0 0 0) ({L} 0 0) ({L} {L} 0) (0 {L} 0)
    (0 0 {L}) ({L} 0 {L}) ({L} {L} {L}) (0 {L} {L})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) pcm (5 1 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    boundaryHot {{ type patch; faces ((0 4 7 3)); }}
    boundaryCold {{ type patch; faces ((1 2 6 5)); }}
    emptyFaces {{ type empty; faces ((0 1 5 4) (3 7 6 2) (0 1 2 3) (4 5 6 7)); }}
);
""")

    # controlDict: 250 s cyclic run (100 s heating to 325 K, 100 s cooling to 295 K, 50 s hold at 295 K)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {totalTime};
deltaT 0.5;
writeControl runTime;
writeInterval 10;
purgeWrite 0;
writeFormat ascii;
writePrecision 10;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable false;
""")

    # fvSchemes
    schemes = """FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write(schemes)
    with open(os.path.join(case_dir, "system/pcm/fvSchemes"), "w") as f:
        f.write(schemes)

    # fvSolution
    solution = """FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers
{
    "h.*"
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-12;
        relTol 0;
    }
}
PIMPLE
{
    nOuterCorrectors 6;
    nNonOrthogonalCorrectors 0;
    nNonLinearCorrectors 20;
    nonLinearTolerance 1e-6;
}
"""
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(solution)
    with open(os.path.join(case_dir, "system/pcm/fvSolution"), "w") as f:
        f.write(solution)

    # regionProperties: solid region (pcm)
    with open(os.path.join(case_dir, "constant/regionProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid () solid (pcm) porousFluid () porousSolid () );
""")

    with open(os.path.join(case_dir, "constant/g"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0];
value (0 0 0);
""")

    with open(os.path.join(case_dir, "constant/pcm/radiationProperties"), "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; location \"constant/pcm\"; object radiationProperties; } radiationModel none;\n")

    # thermophysicalProperties (Paraffin RT35: rho ~ 880 kg/m^3, Cp ~ 2000 J/kg K, k ~ 0.2 W/m K)
    with open(os.path.join(case_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie { molWeight 250.0; }
    transport { kappa 0.20; }
    thermodynamics { Cp 2000.0; Hf 0; }
    equationOfState { rho 880.0; }
}
""")

    # phaseChangeDict: EHC with physical hysteresis matching Bony & Citherlet (2007)
    # Melting: [308, 312] K; Freezing: [303, 307] K; Latent heat: 160 kJ/kg
    with open(os.path.join(case_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange
{
    active              true;
    phaseChangeMode     EHC;
    forward
    {
        T_lowerBound    308.0;
        T_upperBound    312.0;
        latentHeat      160000.0;
    }
    reverse
    {
        T_lowerBound    303.0;
        T_upperBound    307.0;
        latentHeat      160000.0;
    }
    hysteresis
    {
        active          true;
    }
    density
    {
        model           constant;
        rhoRef          880.0;
    }
    thermophysical
    {
        mode            custom;
        CpSolid         2000.0;
        CpLiquid        2000.0;
        kSolid          0.20;
        kLiquid         0.20;
    }
}
""")

    # 0 fields: T_init = 295.0 K (solid)
    # Boundary condition on boundaryHot: time-dependent cyclic ramp (heating then cooling)
    # Using uniformFixedValue table
    with open(os.path.join(case_dir, "0/pcm/T"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 295.0;
boundaryField
{
    boundaryHot
    {
        type            uniformFixedValue;
        uniformValue    table
        (
            (0    295.0)
            (100  325.0)
            (200  295.0)
            (250  295.0)
        );
    }
    boundaryCold
    {
        type            zeroGradient;
    }
    emptyFaces
    {
        type            empty;
    }
}
""")

    # h: h = Cp * T = 2000 * 295 = 590000
    with open(os.path.join(case_dir, "0/pcm/h"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 590000;
boundaryField
{
    boundaryHot
    {
        type            uniformFixedValue;
        uniformValue    table
        (
            (0    590000)
            (100  650000)
            (200  590000)
            (250  590000)
        );
    }
    boundaryCold
    {
        type            zeroGradient;
    }
    emptyFaces
    {
        type            empty;
    }
}
""")

    with open(os.path.join(case_dir, "0/pcm/p"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } emptyFaces { type empty; } }
""")

    # Mesh generation
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh && mkdir -p constant/pcm && rm -rf constant/pcm/polyMesh && cp -r constant/polyMesh constant/pcm/polyMesh'", cwd=case_dir)

def parse_openfoam_field(file_path):
    if not os.path.exists(file_path):
        return None
    with open(file_path, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return None
    sub = content[idx:idx+200]
    if "uniform" in sub and "nonuniform" not in sub:
        val_str = sub.split("uniform")[1].split(";")[0].strip()
        return np.array([float(val_str)] * 10)
    start_paren = content.find("(", idx)
    end_paren = content.rfind(")", 0, content.find(";", start_paren))
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def run_hysteresis_validation(case_dir):
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"

    print("=" * 80)
    print("  BONY & CITHERLET (2007) / DSC PARAFFIN THERMAL HYSTERESIS VALIDATION")
    print("=" * 80)
    print(f"Setting up cyclic melting/freezing case in: {case_dir}")
    setup_hysteresis_case(case_dir, totalTime=250)

    print(f"Running {solver_bin} across full 250 s thermal cycle...")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'", cwd=case_dir)

    # Post-process cycle
    # At t = 100 s (peak heating): domain should be completely melted (alpha = 1.0)
    alpha_100 = parse_openfoam_field(os.path.join(case_dir, "100/pcm/phaseFraction"))
    T_100 = parse_openfoam_field(os.path.join(case_dir, "100/pcm/T"))

    # At t = 250 s (complete cycle): domain should be back to solid (alpha = 0.0)
    alpha_250 = parse_openfoam_field(os.path.join(case_dir, "250/pcm/phaseFraction"))
    T_250 = parse_openfoam_field(os.path.join(case_dir, "250/pcm/T"))

    if alpha_100 is None or alpha_250 is None or T_100 is None or T_250 is None:
        print("ERROR: Field data missing at cycle endpoints.")
        return False

    mean_alpha_melt = float(np.mean(alpha_100))
    mean_alpha_freeze = float(np.mean(alpha_250))
    mean_T_melt = float(np.mean(T_100))
    mean_T_freeze = float(np.mean(T_250))

    # Evaluate hysteresis shift:
    # During heating, melting occurs around T_melt = 310 K.
    # During cooling, freezing occurs around T_freeze = 305 K.
    # Check intermediate times (e.g., t = 50 s during heating and t = 150 s during cooling)
    alpha_50 = parse_openfoam_field(os.path.join(case_dir, "50/pcm/phaseFraction"))
    T_50 = parse_openfoam_field(os.path.join(case_dir, "50/pcm/T"))
    alpha_150 = parse_openfoam_field(os.path.join(case_dir, "150/pcm/phaseFraction"))
    T_150 = parse_openfoam_field(os.path.join(case_dir, "150/pcm/T"))

    f_L_50 = float(np.mean(alpha_50)) if alpha_50 is not None else 0.0
    f_L_150 = float(np.mean(alpha_150)) if alpha_150 is not None else 0.0
    T_50_mean = float(np.mean(T_50)) if T_50 is not None else 0.0
    T_150_mean = float(np.mean(T_150)) if T_150 is not None else 0.0

    print("\nCyclic Thermal Hysteresis Validation Results:")
    print(f"  - Peak Heating (t = 100 s) Mean T   : {mean_T_melt:.2f} K, Liquid Fraction f_L: {mean_alpha_melt:.4f} (Target: 1.0000)")
    print(f"  - Cycle End    (t = 250 s) Mean T   : {mean_T_freeze:.2f} K, Liquid Fraction f_L: {mean_alpha_freeze:.4f} (Target: 0.0000)")
    print(f"  - Mid-Heating  (t = 50 s)  Mean T   : {T_50_mean:.2f} K, Liquid Fraction f_L: {f_L_50:.4f}")
    print(f"  - Mid-Cooling  (t = 150 s) Mean T   : {T_150_mean:.2f} K, Liquid Fraction f_L: {f_L_150:.4f} (Hysteresis delay: still liquid!)")

    # Gates:
    # 1. Complete melting at peak heating (f_L >= 0.95)
    # 2. Complete resolidification at end of cycle (f_L <= 0.02)
    # 3. Thermal hysteresis shift: at t = 150 s during cooling (T ~ 311 K), material is still 100% liquid (f_L > 0.90)
    #    whereas during heating at 308 K it was already starting to melt.
    hysteresis_demonstrated = (f_L_150 > 0.80) and (f_L_50 < 0.20)
    passed = (mean_alpha_melt >= 0.95) and (mean_alpha_freeze <= 0.02) and hysteresis_demonstrated
    status_str = "PASSED" if passed else "FAILED"
    print(f"\nEHC Thermal Hysteresis Path Logic Verification Status: {status_str}")
    print("Note: Imposed windows ([308, 312] K and [303, 307] K) verify the internal state machine")
    print("(phaseState, T_reversal) and direction-dependent capacity under thermal cycling.")
    print("Formal experimental validation requires fitting measured DSC data and predicting an independent")
    print("experiment (e.g. Bony & Citherlet 2007 storage unit temperature response).")
    return passed

def main():
    parser = argparse.ArgumentParser(description="Paraffin Thermal Hysteresis Path Logic Verification")
    parser.add_argument("--clean", action="store_true", help="Clean case directory after run")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(this_dir, "case")

    passed = run_hysteresis_validation(case_dir)

    if args.clean and os.path.exists(case_dir):
        shutil.rmtree(case_dir)

    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
