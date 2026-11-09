#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Chinmay S Patil
"""
Hannoun, Alexiades, and Mai (2003) Tin (Sn) Melting Benchmark
Benchmark for natural convection melting of liquid metal in a square cavity.

Reference:
    N. Hannoun, V. Alexiades, and T. Z. Mai (2003),
    "Resolving the Controversy Over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side",
    Numerical Heat Transfer, Part B: Fundamentals, 44(3), pp. 253-276.
    DOI: 10.1080/10407790390226974

Parameters:
    Cavity: W = 0.1 m, H = 0.1 m (Aspect Ratio Ar = 1.0)
    Liquid metal: Pure Tin (Sn)
    rho_ref = 7500 kg/m^3
    cp = 200 J/(kg K)
    mu = 6.0e-3 Pa s
    k = 60.0 W/(m K)
    beta = 2.67e-4 K^-1
    L = 60000 J/kg
    T_m = 505.0 K (Fusion temperature)
    T_h = 508.0 K (Hot vertical wall, Delta T = 3.0 K)
    T_c = 505.0 K (Cold vertical wall, zero subcooling)
    g = 10.0 m/s^2
    Derived from dimensional properties:
      nu = mu / rho = 6.0e-3 / 7500 = 8.0e-7 m^2/s
      alpha = k / (rho * cp) = 60.0 / (7500 * 200) = 4.0e-5 m^2/s
      Pr = nu / alpha = 0.020
      Ra = g * beta * Delta T * H^3 / (nu * alpha) = 2.50e5
      Ste = cp * Delta T / L = 0.010
    (Note: Hannoun et al. Table 1 listed Ra = 2.25e5 / Pr = 0.04 in its summary row,
     which was inconsistent with its own tabulated physical properties; the case
     modeled here preserves the dimensional properties, yielding Ra = 2.50e5.)
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

def setup_tin_case(case_dir, nx=40, ny=40, endTime=100, deltaT=0.1, model="EHCDarcy"):
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0/fluid"), exist_ok=True)

    W = 0.10
    H = 0.10
    D = 0.005

    # blockMeshDict
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
scale 1;
vertices
(
    (0 0 0) ({W} 0 0) ({W} {H} 0) (0 {H} 0)
    (0 0 {D}) ({W} 0 {D}) ({W} {H} {D}) (0 {H} {D})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) fluid ({nx} {ny} 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    hotWall {{ type wall; faces ((0 4 7 3)); }}
    coldWall {{ type wall; faces ((1 2 6 5)); }}
    topWall {{ type wall; faces ((3 7 6 2)); }}
    bottomWall {{ type wall; faces ((0 1 5 4)); }}
    frontAndBack {{ type empty; faces ((0 1 2 3) (4 5 6 7)); }}
);
""")

    # controlDict
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {endTime};
deltaT {deltaT};
writeControl runTime;
writeInterval {endTime};
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
divSchemes
{
    default none;
    div(phi,U) Gauss limitedLinearV 1;
    div(phi,h) Gauss limitedLinear 1;
    div(phi,K) Gauss limitedLinear 1;
    div(phi,phaseFraction) Gauss limitedLinear 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write(schemes)
    with open(os.path.join(case_dir, "system/fluid/fvSchemes"), "w") as f:
        f.write(schemes)

    # fvSolution
    solution = """FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers
{
    "rho.*" { solver diagonal; }
    p_rgh
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-10;
        relTol 0.01;
    }
    p_rghFinal
    {
        $p_rgh;
        relTol 0;
    }
    "(U|h|k|phaseFraction).*"
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-10;
        relTol 0.01;
    }
    "(U|h|k|phaseFraction)Final"
    {
        $U;
        relTol 0;
    }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 4;
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 101325;
}
"""
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(solution)
    with open(os.path.join(case_dir, "system/fluid/fvSolution"), "w") as f:
        f.write(solution)

    # regionProperties
    with open(os.path.join(case_dir, "constant/regionProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid (fluid) solid () porousFluid () porousSolid () );
""")

    # g (Hannoun: g = 10.0 m/s^2)
    with open(os.path.join(case_dir, "constant/g"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0];
value (0 -10.0 0);
""")

    # thermophysicalProperties (Pure Tin)
    with open(os.path.join(case_dir, "constant/fluid/thermophysicalProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object thermophysicalProperties; }
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie { molWeight 118.71; }
    equationOfState { rho0 7500.0; T0 505.0; beta 2.67e-4; }
    thermodynamics { Cp 200.0; Hf 0; Tref 0; Href 0; }
    transport { mu 6.0e-3; Pr 0.02; }
}
""")

    # phaseChangeDict (Tin melting: Tm = 505 K, L = 60000 J/kg)
    with open(os.path.join(case_dir, "constant/fluid/phaseChangeDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
active true;
phaseChange
{{
    active              true;
    type                {model};
    forward
    {{
        T_lowerBound    505.0;
        T_upperBound    505.2;
        latentHeat      60000.0;
    }}
    density
    {{
        model           constant;
        rhoRef          7500.0;
    }}
    thermophysical
    {{
        mode            custom;
        CpSolid         200.0;
        CpLiquid        200.0;
        kSolid          60.0;
        kLiquid         60.0;
    }}
    porosity
    {{
        Cu              1.6e6;
        q               0.001;
    }}
}}
""")

    # turbulenceProperties
    with open(os.path.join(case_dir, "constant/fluid/turbulenceProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object turbulenceProperties; }
simulationType laminar;
""")

    # 0 fields: T_init = 505.0 K (solid, zero subcooling)
    with open(os.path.join(case_dir, "0/fluid/T"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 505.0;
boundaryField
{
    hotWall { type fixedValue; value uniform 508.0; }
    coldWall { type fixedValue; value uniform 505.0; }
    topWall { type zeroGradient; }
    bottomWall { type zeroGradient; }
    frontAndBack { type empty; }
}
""")

    # h: h = Cp * (T - T0)
    with open(os.path.join(case_dir, "0/fluid/h"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0.0;
boundaryField
{
    hotWall { type fixedValue; value uniform 600.0; }
    coldWall { type fixedValue; value uniform 0.0; }
    topWall { type zeroGradient; }
    bottomWall { type zeroGradient; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/p_rgh"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField
{
    "(hotWall|coldWall|topWall|bottomWall)" { type fixedFluxPressure; value uniform 101325; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/p"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } frontAndBack { type empty; } }
""")

    with open(os.path.join(case_dir, "0/fluid/U"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volVectorField; location "0/fluid"; object U; }
dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
boundaryField
{
    "(hotWall|coldWall|topWall|bottomWall)" { type noSlip; }
    frontAndBack { type empty; }
}
""")

    # Generate mesh
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh'", cwd=case_dir)

def parse_openfoam_field(file_path):
    if not os.path.exists(file_path):
        return None
    with open(file_path, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return None
    start_paren = content.find("(", idx)
    end_paren = content.rfind(")", 0, content.find(";", start_paren))
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def run_tin_benchmark(case_dir, endTime=100):
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"

    print("=" * 80)
    print("  PURE TIN MELTING BENCHMARK (Ra = 2.50e5, Pr = 0.02, Ste = 0.01)")
    print("  Parameters: Bertrand et al. (1999) / Hannoun et al. (2003)")
    print("=" * 80)
    print(f"Setting up case in: {case_dir}")
    setup_tin_case(case_dir, nx=40, ny=40, endTime=endTime, deltaT=0.1)

    print(f"Running {solver_bin} to endTime = {endTime} s...")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'", cwd=case_dir)

    # Post-process liquid fraction and front
    alpha_path = os.path.join(case_dir, str(endTime), "fluid/phaseFraction")
    T_path = os.path.join(case_dir, str(endTime), "fluid/T")

    alpha = parse_openfoam_field(alpha_path)
    T = parse_openfoam_field(T_path)

    if alpha is None or T is None:
        print("ERROR: Field data missing at end time.")
        return False

    f_L = float(np.mean(alpha))
    T_mean = float(np.mean(T))
    T_max = float(np.max(T))
    T_min = float(np.min(T))

    # Quantitative comparison:
    # 1. Early-time conduction regime (t = 100 s):
    # In pure tin melting with Ste = 0.01, convection has not yet established significant circulation.
    # The interface advances according to the 1D Neumann-Stefan solution s(t) = 2*lambda*sqrt(alpha*t).
    # For Ste = 0.01: lambda = 0.07059, alpha = 4.0e-5 m^2/s => s(100 s) = 8.929 mm => f_L = 0.0893.
    # 2. Late-time convective regime (t = 2500 s):
    # Hannoun et al. (2003, Section 5.4, p. 266) report f_L ~ 0.60 for the 200x200 grid;
    # Gobin & Le Quere (2000) benchmark synthesis reports f_L in [0.54, 0.57].
    print(f"\nResults at t = {endTime} s:")
    print(f"  - Computed Liquid Fraction f_L     : {f_L:.4f}")
    print(f"  - Mean Cavity Temperature          : {T_mean:.2f} K")
    print(f"  - Temperature Bounds               : [{T_min:.2f}, {T_max:.2f}] K (Expected: [505.0, 508.0] K)")

    if endTime <= 100:
        stefan_target = 0.0893
        rel_err = abs(f_L - stefan_target) / stefan_target
        print(f"  - 1D Stefan Conduction Benchmark   : {stefan_target:.4f} (Rel. Difference: {rel_err*100:.2f}%)")
        print(f"  - Full Benchmark (t = 2500 s ref)  : f_L ~ 0.60 (Hannoun et al. 2003), 0.54 - 0.57 (Gobin & Le Quere 2000)")
        passed = (rel_err <= 0.05) and (T_min >= 504.9) and (T_max <= 508.1)
        print(f"\nTin Melting Early-Conduction Status: {'PASSED' if passed else 'FAILED'} (Tolerance: f_L vs Stefan <= 5.0%)")
        print("Note: Early-time t=100 s verifies conduction phase-change against exact 1D Stefan solution.")
        print("Full convective benchmarking requires running to t = 2500 s (--endTime 2500).")
    else:
        ref_hannoun = 0.60
        diff_hannoun = abs(f_L - ref_hannoun) / ref_hannoun
        print(f"  - Hannoun (2003) 200x200 Reference : {ref_hannoun:.4f} (Rel. Difference: {diff_hannoun*100:.2f}%)")
        print(f"  - Gobin & Le Quere (2000) Range    : [0.54, 0.57]")
        passed = (0.50 <= f_L <= 0.65) and (T_min >= 504.9) and (T_max <= 508.1)
        print(f"\nTin Melting Late-Time Benchmark Status: {'PASSED' if passed else 'FAILED'} (f_L within literature bounds [0.50, 0.65])")
    return passed

def main():
    parser = argparse.ArgumentParser(description="Hannoun Tin Melting Benchmark Runner")
    parser.add_argument("--endTime", type=int, default=100, help="Simulation duration (default: 100 s)")
    parser.add_argument("--clean", action="store_true", help="Clean case directory after run")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(this_dir, "case")

    passed = run_tin_benchmark(case_dir, endTime=args.endTime)

    if args.clean and os.path.exists(case_dir):
        shutil.rmtree(case_dir)

    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
