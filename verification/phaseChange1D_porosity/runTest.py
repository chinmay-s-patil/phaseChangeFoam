#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil

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
        os.path.join(home, "OpenFOAM", f"lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam")
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
        return [float(val_str)] * num_cells
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    block = content[start_paren+1:end_paren].strip()
    return [float(x) for x in block.split()]

def setup_case(case_dir, model_type):
    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 200; deltaT 2;
writeControl runTime; writeInterval 2; purgeWrite 0; writeFormat ascii;
""")
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers {
    "h.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    "T.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
}
PIMPLE { nOuterCorrectors 2; nNonLinearCorrectors 50; nonLinearTolerance 1e-6; }
""")
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }
scale 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) (100 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot { type patch; faces ((0 4 7 3)); }
    cold { type patch; faces ((1 2 6 5)); }
    emptyFaces { type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }
);
""")

    # constant
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    with open(os.path.join(case_dir, "constant/g"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0]; value (0 0 0);
""")
    with open(os.path.join(case_dir, "constant/regionProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid () solid (pcm) porousFluid () porousSolid () );
""")
    
    os.makedirs(os.path.join(case_dir, "constant/pcm"), exist_ok=True)
    with open(os.path.join(case_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture {
    specie { molWeight 200.0; }
    transport { kappa 0.50; }
    thermodynamics { Cp 1980.0; Hf 0; }
    equationOfState { rho 1967.0; }
}
""")
    with open(os.path.join(case_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        if model_type == "ehc":
            f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    phaseChangeMode EHC;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
        else:
            f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    phaseChangeMode enthalpyPorosity;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")

    # 0 fields
    os.makedirs(os.path.join(case_dir, "0/pcm"), exist_ok=True)
    with open(os.path.join(case_dir, "0/pcm/T"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField {
    hot { type fixedValue; value uniform 350.0; }
    cold { type zeroGradient; }
    emptyFaces { type empty; }
}
""")
    with open(os.path.join(case_dir, "0/pcm/h"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 554400;
boundaryField {
    hot { type fixedValue; value uniform 693000; }
    cold { type zeroGradient; }
    emptyFaces { type empty; }
}
""")
    with open(os.path.join(case_dir, "0/pcm/p"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } emptyFaces { type empty; } }
""")

    run_cmd(f"cd {case_dir} && bash -c '{of_env}; blockMesh'", cwd=case_dir)
    run_cmd(f"mkdir -p {case_dir}/constant/pcm {case_dir}/system/pcm")
    run_cmd(f"cp -r {case_dir}/constant/polyMesh {case_dir}/constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSchemes {case_dir}/system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSolution {case_dir}/system/pcm/fvSolution 2>/dev/null || true")

def compute_backward_euler_ein(case_dir, delta_t=2.0):
    time_dirs = sorted([int(d) for d in os.listdir(case_dir) if d.isdigit() and int(d) > 0])
    E_in = 0.0
    for t in time_dirs:
        T_snap = parse_openfoam_field(os.path.join(case_dir, str(t), "pcm/T"))
        if T_snap:
            # End-of-step flux: q_snap = k * (350 - T_0^n) / (dx/2) * Area
            q_snap = 0.50 * (0.01 * 0.01) * (350.0 - T_snap[0]) / (0.0005)
            E_in += q_snap * delta_t
    return E_in

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== EHC vs Enthalpy-Porosity 1D Melting Comparison in {base_dir} ===")
    
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    ehc_dir = os.path.join(base_dir, "case_ehc")
    ep_dir = os.path.join(base_dir, "case_porosity")

    shutil.rmtree(ehc_dir, ignore_errors=True)
    shutil.rmtree(ep_dir, ignore_errors=True)

    print("\n--- Running EHC Model Case ---")
    setup_case(ehc_dir, "ehc")
    run_cmd(f"cd {ehc_dir} && bash -c '{of_env}; {solver_bin}'")
    T_ehc = parse_openfoam_field(os.path.join(ehc_dir, "200/pcm/T"))
    a_ehc = parse_openfoam_field(os.path.join(ehc_dir, "200/pcm/phaseFraction"))

    print("\n--- Running Enthalpy-Porosity Model Case ---")
    setup_case(ep_dir, "enthalpyPorosity")
    run_cmd(f"cd {ep_dir} && bash -c '{of_env}; {solver_bin}'")
    T_ep = parse_openfoam_field(os.path.join(ep_dir, "200/pcm/T"))
    a_ep = parse_openfoam_field(os.path.join(ep_dir, "200/pcm/phaseFraction"))

    if not T_ehc or not T_ep:
        print("FAILED: Could not parse fields.")
        sys.exit(1)

    mean_T_ehc = sum(T_ehc) / len(T_ehc)
    mean_a_ehc = sum(a_ehc) / len(a_ehc)
    mean_T_ep = sum(T_ep) / len(T_ep)
    mean_a_ep = sum(a_ep) / len(a_ep)

    dT_diff = sum(abs(t1 - t2) for t1, t2 in zip(T_ehc, T_ep)) / len(T_ehc)
    da_diff = sum(abs(a1 - a2) for a1, a2 in zip(a_ehc, a_ep)) / len(a_ehc)

    # Domain enthalpy rise check: deltaH = sum(rho * V * (Cp * (T - 280) + L * alpha))
    rho = 1967.0
    Cp = 1980.0
    L = 163000.0
    V_cell = 0.1 / 100.0 * 0.01 * 0.01

    dH_ehc = sum(rho * V_cell * (Cp * (t - 280.0) + L * a) for t, a in zip(T_ehc, a_ehc))
    dH_ep = sum(rho * V_cell * (Cp * (t - 280.0) + L * a) for t, a in zip(T_ep, a_ep))
    rel_dH_diff = abs(dH_ehc - dH_ep) / max(dH_ehc, 1e-10)

    # Backward Euler end-of-step flux integration: E_in = sum( q^n * dt )
    Ein_ehc = compute_backward_euler_ein(ehc_dir, delta_t=2.0)
    Ein_ep = compute_backward_euler_ein(ep_dir, delta_t=2.0)

    rel_flux_bal_ehc = abs(dH_ehc - Ein_ehc) / max(Ein_ehc, 1e-10)
    rel_flux_bal_ep = abs(dH_ep - Ein_ep) / max(Ein_ep, 1e-10)

    pass_ehc_melting = (mean_a_ehc > 0.01)
    pass_ep_melting = (mean_a_ep > 0.01)
    pass_match = (rel_dH_diff < 0.001 and dT_diff < 0.01 and da_diff < 0.001)
    pass_flux_balance = (rel_flux_bal_ehc < 0.001 and rel_flux_bal_ep < 0.001)

    print("\n=======================================================")
    print("      EHC vs ENTHALPY-POROSITY COMPARISON RESULTS      ")
    print("=======================================================")
    print(f"EHC Model (t=200s) Mean T          : {mean_T_ehc:.2f} K, Mean alpha = {mean_a_ehc:.4f}")
    print(f"Porosity Model (t=200s) Mean T     : {mean_T_ep:.2f} K, Mean alpha = {mean_a_ep:.4f}")
    print(f"EHC Domain Enthalpy Rise (dH)      : {dH_ehc:.4f} J")
    print(f"Porosity Domain Enthalpy Rise (dH) : {dH_ep:.4f} J")
    print(f"EHC Backward Euler Flux (E_in)     : {Ein_ehc:.4f} J (Balance Error: {rel_flux_bal_ehc*100:.6f}%)")
    print(f"Porosity Backward Euler Flux (E_in): {Ein_ep:.4f} J (Balance Error: {rel_flux_bal_ep*100:.6f}%)")
    print(f"EHC vs Porosity Enthalpy Diff     : {rel_dH_diff*100:.4f}% (Limit < 0.10%)")
    print(f"Temperature Mean Difference       : {dT_diff:.6f} K (Limit < 0.01 K)")
    print(f"Phase Fraction Mean Difference    : {da_diff:.8f} (Limit < 0.001)")
    print("-------------------------------------------------------")
    print(f"EHC Active Melting                 : {'PASS' if pass_ehc_melting else 'FAIL'}")
    print(f"Porosity Active Melting            : {'PASS' if pass_ep_melting else 'FAIL'}")
    print(f"EHC & Porosity Solution Agreement  : {'PASS' if pass_match else 'FAIL'}")
    print(f"Backward Euler Flux Conservation   : {'PASS' if pass_flux_balance else 'FAIL'} (Limit < 0.10%)")

    all_pass = pass_ehc_melting and pass_ep_melting and pass_match and pass_flux_balance
    if all_pass:
        print("\nSTATUS: EHC VS ENTHALPY-POROSITY COMPARISON TEST PASSED!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
