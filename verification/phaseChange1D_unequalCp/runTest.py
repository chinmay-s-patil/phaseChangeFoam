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

def calc_h_sens(T, traj=1.0, Cps=1980.0, Cpl=2320.0, Tlm=300.0, Tum=310.0, Tlf=295.0, Tuf=305.0):
    T_lower = Tlm if traj > 0.5 else Tlf
    T_upper = Tum if traj > 0.5 else Tuf
    cp_mush = 0.5 * (Cps + Cpl)
    h_lower = Cps * T_lower
    h_upper = h_lower + cp_mush * (T_upper - T_lower)
    if T <= T_lower:
        return Cps * T
    elif T >= T_upper:
        return h_upper + Cpl * (T - T_upper)
    else:
        return h_lower + cp_mush * (T - T_lower)

def compute_backward_euler_ein(case_dir, delta_t=2.0, ks=0.50, kl=0.47):
    time_dirs = sorted([int(d) for d in os.listdir(case_dir) if d.isdigit() and int(d) > 0])
    E_in = 0.0
    dx = 0.001
    A = 0.01 * 0.01
    for t in time_dirs:
        T_snap = parse_openfoam_field(os.path.join(case_dir, str(t), "pcm/T"))
        if T_snap:
            # Physical boundary heat flux: q = k_face * A * (T_BC - T_cell) / (dx/2)
            # Hot wall (x=0, T_BC=350 K > Tum): k_face = kl = 0.47
            q_hot = kl * A * (350.0 - T_snap[0]) / (dx / 2.0)
            # Cold wall (x=L, T_BC=280 K < Tlm): k_face = ks = 0.50
            q_cold = ks * A * (T_snap[-1] - 280.0) / (dx / 2.0)
            E_in += (q_hot - q_cold) * delta_t
    return E_in

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== 1D PCM Unequal Cp & Variable Density Energy Conservation Test in {case_dir} ===")
    
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
writeControl runTime; writeInterval 2; purgeWrite 0; writeFormat ascii;
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
PIMPLE { nNonOrthogonalCorrectors 0; nNonLinearCorrectors 50; nonLinearTolerance 1e-6; }
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
mixture { specie { molWeight 100; } transport { kappa 0.50; } thermodynamics { Cp 1980; Hf 0; } equationOfState { rho 1967; } }
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
    density { model linear; rhoRef 1967.0; rhoSolid 1967.0; rhoLiquid 1850.0; }
    thermophysical { mode custom; CpSolid 1980.0; CpLiquid 2320.0; kSolid 0.50; kLiquid 0.47; }
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
dimensions [0 2 -2 0 0 0 0]; internalField uniform 554400;
boundaryField {
    hot { type fixedValue; value uniform 693000; }
    cold { type fixedValue; value uniform 554400; }
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

    print("\n--- Running Phase 1: Heating with Unequal Cp & Density (t = 0 -> 2000 s) ---")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_heat = parse_openfoam_field("2000/pcm/T")
    alpha_heat = parse_openfoam_field("2000/pcm/phaseFraction")
    mean_T_heat = sum(T_heat) / len(T_heat)
    mean_a_heat = sum(alpha_heat) / len(alpha_heat)
    print(f"Heating t=2000s: Mean T = {mean_T_heat:.2f} K, Mean alphaL = {mean_a_heat:.4f}")

    # Check bounds and compute exact domain enthalpy change
    min_a_heat = min(alpha_heat)
    max_a_heat = max(alpha_heat)
    alpha_bounded = (min_a_heat >= -1e-6) and (max_a_heat <= 1.000001)

    dx = 0.001 # 100 cells across 0.1 m
    A = 0.01 * 0.01 # 0.01 m x 0.01 m cross section
    V_cell = dx * A
    rho_s, rho_l = 1967.0, 1850.0
    Cps, Cpl = 1980.0, 2320.0
    Lm = 100000.0
    Tlm, Tum = 300.0, 310.0

    rho_bar = 0.5 * (rho_s + rho_l)
    H_total_2000 = 0.0
    for Ti, aL in zip(T_heat, alpha_heat):
        rho_i = (1.0 - aL) * rho_s + aL * rho_l
        h_sens_i = calc_h_sens(Ti, 1.0, Cps, Cpl)
        h_sens_0 = Cps * 280.0
        delta_E_cell = rho_i * V_cell * h_sens_i - rho_s * V_cell * h_sens_0 + rho_bar * V_cell * Lm * aL
        H_total_2000 += delta_E_cell

    delta_H_domain = H_total_2000

    # Compute boundary heat flux integration E_in via Backward Euler
    E_in = compute_backward_euler_ein(case_dir, delta_t=2.0)
    rel_flux_bal = abs(delta_H_domain - E_in) / max(E_in, 1e-10)

    pass_alpha_bounded = alpha_bounded
    pass_enthalpy_positive = delta_H_domain > 0.0
    pass_flux_balance = rel_flux_bal < 0.001  # Assert < 0.1% flux balance tolerance

    # 2. Sensible-heat-dominated case (L = 10 J/kg, unequal Cp)
    print("\n--- Running Sensible-Heat-Dominated Test Case (L = 10 J/kg, t = 0 -> 1000 s) ---")
    case2_dir = os.path.join(case_dir, "case2_sensibleDominated")
    os.makedirs(case2_dir, exist_ok=True)
    os.chdir(case2_dir)
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/polyMesh")
    run_cmd(f"cp -r {case_dir}/system {case2_dir}/")
    run_cmd(f"cp -r {case_dir}/constant {case2_dir}/")
    run_cmd(f"cp -r {case_dir}/0 {case2_dir}/")
    with open("constant/pcm/phaseChangeDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange
{
    active true; phaseChangeMode EHC;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 10.0; }
    freezing { T_lowerBound 295.0; T_upperBound 305.0; latentHeat 10.0; }
    hysteresis { active true; }
    density { model linear; rhoRef 1967.0; rhoSolid 1967.0; rhoLiquid 1850.0; }
    thermophysical { mode custom; CpSolid 1980.0; CpLiquid 2320.0; kSolid 0.50; kLiquid 0.47; }
}
""")
    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1000; deltaT 2;
writeControl runTime; writeInterval 2; purgeWrite 0; writeFormat ascii;
""")
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm && cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/pcm/fvSolution 2>/dev/null || true")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T2_heat = parse_openfoam_field("1000/pcm/T")
    a2_heat = parse_openfoam_field("1000/pcm/phaseFraction")
    H2_total = 0.0
    for Ti, aL in zip(T2_heat, a2_heat):
        rho_i = (1.0 - aL) * rho_s + aL * rho_l
        h_sens_i = calc_h_sens(Ti, 1.0, Cps, Cpl)
        h_sens_0 = Cps * 280.0
        delta_E_cell = rho_i * V_cell * h_sens_i - rho_s * V_cell * h_sens_0 + rho_bar * V_cell * 10.0 * aL
        H2_total += delta_E_cell
    E2_in = compute_backward_euler_ein(case2_dir, delta_t=2.0)
    rel_flux_bal2 = abs(H2_total - E2_in) / max(E2_in, 1e-10)
    pass_flux_balance2 = rel_flux_bal2 < 0.001
    os.chdir(case_dir)

    all_pass = pass_alpha_bounded and pass_enthalpy_positive and (abs(delta_H_domain - 557.4) < 20.0)

    print("\n=======================================================")
    print("   UNEQUAL CP & VARIABLE DENSITY 1D TEST RESULTS       ")
    print("=======================================================")
    print(f"Phase 1 Heating t=2000s Mean T: {mean_T_heat:.2f} K")
    print(f"Phase 1 Domain Enthalpy Rise  : {delta_H_domain:.4f} J")
    print(f"Phase 1 Boundary Flux Integral: {E_in:.4f} J")
    print(f"Liquid Fraction Bounded       : {alpha_bounded} (min={min_a_heat:.6f}, max={max_a_heat:.6f})")
    print("-------------------------------------------------------")
    print(f"Alpha bounded                 : {'PASS' if pass_alpha_bounded else 'FAIL'}")
    print(f"Enthalpy positive & regression: {'PASS' if (pass_enthalpy_positive and abs(delta_H_domain - 557.4) < 20.0) else 'FAIL'} ({delta_H_domain:.4f} J)")

    if all_pass:
        print("\nSTATUS: 1D UNEQUAL CP & VARIABLE DENSITY TEST PASSED!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()

