#!/usr/bin/env python3
import os
import sys
import subprocess
import math
import shutil
import re

def find_solver():
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        path = os.path.join(user_appbin, "phaseChangeMultiRegionFoam")
        if os.path.exists(path):
            return path
    path = shutil.which("phaseChangeMultiRegionFoam")
    if path:
        return path
    home = os.environ.get("HOME", "/home/lavender")
    matches = [
        os.path.join(home, "OpenFOAM", "lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam")
    ]
    for m in matches:
        if os.path.exists(m):
            return m
    return "phaseChangeMultiRegionFoam"

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing command: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res.stdout

def parse_openfoam_field(file_path, num_cells=100):
    if not os.path.exists(file_path):
        return []
    with open(file_path, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return []
    sub = content[idx:idx+200]
    if "uniform" in sub and "nonuniform" not in sub:
        val_str = sub.split("uniform")[1].split(";")[0].strip()
        val = float(val_str)
        return [val] * num_cells
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    block = content[start_paren+1:end_paren].strip()
    return [float(x) for x in block.split()]

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== 1D PCM Heating-Cooling Reversal Verification Test in {case_dir} ===")
    
    os.chdir(case_dir)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()

    # 1. Setup base 1D mesh
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/polyMesh")
    
    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices ((0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01));
blocks ( hex (0 1 2 3 4 5 6 7) pcm (100 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot { type patch; faces ((0 3 7 4)); }
    cold { type patch; faces ((1 5 6 2)); }
    bottom { type empty; faces ((0 1 5 4)); }
    top { type empty; faces ((3 2 6 7)); }
    back { type empty; faces ((0 1 2 3)); }
    front { type empty; faces ((4 5 6 7)); }
);
""")

    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 2000; deltaT 2;
writeControl runTime; writeInterval 100; purgeWrite 0; writeFormat ascii;
""")

    with open("system/fvSchemes", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")

    with open("system/fvSolution", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
PIMPLE { nNonOrthogonalCorrectors 0; nNonLinearCorrectors 25; nonLinearTolerance 1e-5; }
solvers { "h.*" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; } }
""")

    os.makedirs("system/pcm", exist_ok=True)
    with open("system/pcm/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/pcm/fvSolution", "w") as f: f.write(open("system/fvSolution").read())

    os.makedirs("constant/pcm", exist_ok=True)
    with open("constant/regionProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid () solid (pcm) porousFluid () porousSolid () );
""")
    with open("constant/g", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0]; value (0 0 0);
""")
    with open("constant/pcm/radiationProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object radiationProperties; }
radiationModel none;
""")
    with open("constant/pcm/thermophysicalProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture { specie { molWeight 100; } transport { kappa 1.0; } thermodynamics { Cp 2000; Hf 0; } equationOfState { rho 1000; } }
""")
    with open("constant/pcm/phaseChangeDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange
{
    active true; phaseChangeMode EHC;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    freezing { T_lowerBound 295.0; T_upperBound 305.0; latentHeat 100000.0; }
    hysteresis { active true; }
    density { model linear; rhoRef 1000.0; rhoSolid 1000.0; rhoLiquid 1000.0; }
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 1.0; kLiquid 1.0; }
}
""")

    os.makedirs("0/pcm", exist_ok=True)
    with open("0/pcm/T", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField {
    hot { type fixedValue; value uniform 350.0; }
    cold { type fixedValue; value uniform 280.0; }
    "(top|bottom|front|back)" { type empty; }
}
""")
    with open("0/pcm/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField {
    hot { type fixedValue; value uniform 700000; }
    cold { type fixedValue; value uniform 560000; }
    "(top|bottom|front|back)" { type empty; }
}
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } "(top|bottom|front|back)" { type empty; } }
""")

    # Run Phase 1: Heating (0 -> 2000s)
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/pcm/fvSolution 2>/dev/null || true")

    print("\n--- Running Phase 1: Heating (t = 0 -> 2000 s, T_hot = 350 K) ---")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_heat = parse_openfoam_field("2000/pcm/T")
    alpha_heat = parse_openfoam_field("2000/pcm/phaseFraction")
    mean_T_heat = sum(T_heat) / len(T_heat)
    mean_a_heat = sum(alpha_heat) / len(alpha_heat)
    print(f"Heating t=2000s: Mean T = {mean_T_heat:.2f} K, Mean alphaL = {mean_a_heat:.4f}")

    # Run Phase 2: Cooling Reversal (2000 -> 4000s) with T_hot switched to 270 K
    print("\n--- Running Phase 2: Cooling Reversal (t = 2000 -> 4000 s, T_hot = 270 K) ---")
    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime 4000; deltaT 2;
writeControl runTime; writeInterval 100; purgeWrite 0; writeFormat ascii;
""")

    # update boundary BC at t=2000 using exact regex targeting hot patch
    with open("2000/pcm/T", "r") as f:
        t_2000 = f.read()
    t_2000_cool = re.sub(r"(hot\s*\{\s*type\s+fixedValue;\s*value\s+uniform\s+)[0-9.]+(;\s*\})", r"\g<1>270.0\2", t_2000)
    assert t_2000 != t_2000_cool, "Regex substitution failed to update hot patch boundary value!"
    with open("2000/pcm/T", "w") as f:
        f.write(t_2000_cool)

    run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_cool = parse_openfoam_field("4000/pcm/T")
    alpha_cool = parse_openfoam_field("4000/pcm/phaseFraction")
    mean_T_cool = sum(T_cool) / len(T_cool)
    mean_a_cool = sum(alpha_cool) / len(alpha_cool)
    print(f"Cooling t=4000s: Mean T = {mean_T_cool:.2f} K, Mean alphaL = {mean_a_cool:.4f}")

    # Check bounds: liquid fraction must strictly be in [0, 1]
    min_a_cool = min(alpha_cool)
    max_a_cool = max(alpha_cool)
    alpha_bounded = (min_a_cool >= -1e-6) and (max_a_cool <= 1.000001)

    # Regression check against hysteresis OFF heating benchmark (mean_a = 0.2887, mean_T = 299.37 K)
    reg_a_passed = abs(mean_a_heat - 0.2887) < 0.002
    reg_T_passed = abs(mean_T_heat - 299.37) < 0.1

    # Assert cooling mean T at 4000s is 286.02 K (+/- 0.1 K)
    reg_T_cool_passed = abs(mean_T_cool - 286.02) < 0.10

    # Enthalpy balance
    rho = 1000.0
    Cp = 2000.0
    Lf = 100000.0
    T0 = 280.0
    dx = 0.1 / len(T_heat)
    H_heat = sum(rho * (Cp * (T - T0) + a * Lf) * dx for T, a in zip(T_heat, alpha_heat))
    H_cool = sum(rho * (Cp * (T - T0) + a * Lf) * dx for T, a in zip(T_cool, alpha_cool))
    delta_H_extracted = H_heat - H_cool
    delta_H_MJ = delta_H_extracted / 1e6

    # Extracted cooling enthalpy benchmark: 5.558 MJ/m^2 (expected 5.56 MJ/m^2 +/- 0.5%)
    reg_H_cool_passed = abs(delta_H_MJ - 5.558) / 5.558 < 0.005

    print("\n=======================================================")
    print("      HEATING & COOLING REVERSAL TEST RESULTS          ")
    print("=======================================================")
    print(f"Heating Peak (t=2000s) Mean T     : {mean_T_heat:.2f} K (Benchmark: 299.37 K, Pass={reg_T_passed})")
    print(f"Heating Peak (t=2000s) Mean alphaL : {mean_a_heat:.4f} (Benchmark: 0.2887, Pass={reg_a_passed})")
    print(f"Heating Peak (t=2000s) Domain H    : {H_heat:.2f} J/m^2")
    print(f"Cooling End  (t=4000s) Mean T     : {mean_T_cool:.2f} K (Benchmark: 286.02 K, Pass={reg_T_cool_passed})")
    print(f"Cooling End  (t=4000s) Mean alphaL : {mean_a_cool:.4f}")
    print(f"Cooling End  (t=4000s) Domain H    : {H_cool:.2f} J/m^2")
    print(f"Cooling Enthalpy Extracted         : {delta_H_extracted:.2f} J/m^2 ({delta_H_MJ:.4f} MJ/m^2, Pass={reg_H_cool_passed})")
    print(f"Liquid Fraction Bounded [0, 1]     : {alpha_bounded} (min={min_a_cool:.6f}, max={max_a_cool:.6f})")
    print("-------------------------------------------------------")

    if reg_a_passed and reg_T_passed and reg_T_cool_passed and reg_H_cool_passed and alpha_bounded:
        print("\nSTATUS: HEATING-COOLING REVERSAL VERIFICATION PASSED!")
        print("The solver committed alphaL_old statefully without edge shortcut artifacts, matching 286.02 K cooling mean T and 5.56 MJ/m^2 extracted enthalpy.")
    else:
        print("\nSTATUS: TEST FAILED")
        if not reg_a_passed or not reg_T_passed:
            print(f"Reason: Heating phase regression failed (Mean alphaL={mean_a_heat:.4f} vs 0.2887, Mean T={mean_T_heat:.2f} K vs 299.37 K)")
        if not reg_T_cool_passed:
            print(f"Reason: Cooling mean T at 4000s failed (Mean T={mean_T_cool:.2f} K vs 286.02 K)")
        if not reg_H_cool_passed:
            print(f"Reason: Cooling enthalpy extracted failed ({delta_H_MJ:.4f} MJ/m^2 vs 5.558 MJ/m^2)")
        if not alpha_bounded:
            print(f"Reason: Liquid fraction out of bounds [0, 1]: min={min_a_cool}, max={max_a_cool}")
        sys.exit(1)

if __name__ == "__main__":
    main()
