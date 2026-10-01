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

def parse_openfoam_field(file_path):
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
        return [float(val_str)]
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    block = content[start_paren+1:end_paren].strip()
    return [float(x) for x in block.split()]

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== 1D PCM Fixed Heat Flux BC & Conservation Test in {case_dir} ===")
    
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

    t_end = 1000.0
    dt = 1.0
    q_flux = 500.0 # W/m^2
    rho = 1000.0 # kg/m^3
    Cp = 2000.0 # J/(kg K)
    Lf = 100000.0 # J/kg
    T0 = 280.0 # K
    L = 0.1 # m

    with open("system/controlDict", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {int(t_end)}; deltaT {dt};
writeControl runTime; writeInterval {int(t_end)}; purgeWrite 0; writeFormat ascii;
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
    freezing { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    hysteresis { active false; }
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
    hot { type fixedGradient; gradient uniform 500.0; }
    cold { type zeroGradient; }
    "(top|bottom|front|back)" { type empty; }
}
""")
    with open("0/pcm/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField {
    hot { type fixedGradient; gradient uniform 1000000.0; }
    cold { type zeroGradient; }
    "(top|bottom|front|back)" { type empty; }
}
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } "(top|bottom|front|back)" { type empty; } }
""")

    # Mesh and setup
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/pcm/fvSolution 2>/dev/null || true")

    print(f"\n--- Running Solver with Neumann Fixed Heat Flux q = {q_flux} W/m^2 (t = 0 -> {int(t_end)} s) ---")
    run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    out_dir = f"{int(t_end)}/pcm"
    T_vals = parse_openfoam_field(f"{out_dir}/T")
    alpha_vals = parse_openfoam_field(f"{out_dir}/liquidFraction")
    N = len(T_vals)
    dx = L / N

    # Calculate domain energy change per unit cross-sectional area:
    # E_domain = sum_i rho * [ Cp * (T_i - T0) + alpha_i * Lf ] * dx
    if len(alpha_vals) == 1:
        alpha_vals = alpha_vals * len(T_vals)
    print("len(T_vals):", len(T_vals), "len(alpha_vals):", len(alpha_vals))
    print("T_vals min/max:", min(T_vals), max(T_vals))
    print("alpha_vals min/max:", min(alpha_vals), max(alpha_vals))
    delta_H_domain = sum(rho * (Cp * (T - T0) + a * Lf) * dx for T, a in zip(T_vals, alpha_vals))
    
    # Total input energy per unit area: E_in = q * t_end
    E_in = q_flux * t_end

    rel_err = abs(delta_H_domain - E_in) / E_in

    print(f"\nEnergy Conservation Analysis at t = {t_end} s:")
    print(f"Total Input Heat Energy (E_in)     : {E_in:.2f} J/m^2")
    print(f"Domain Enthalpy Rise (delta_H)     : {delta_H_domain:.2f} J/m^2")
    print(f"Absolute Energy Conservation Error  : {abs(delta_H_domain - E_in):.4f} J/m^2")
    print(f"Relative Energy Error               : {rel_err*100:.6f}%")

    print("\n=======================================================")
    print("      FIXED HEAT FLUX BC & CONSERVATION RESULTS        ")
    print("=======================================================")
    print(f"E_in  : {E_in:.2f} J/m^2")
    print(f"delta_H: {delta_H_domain:.2f} J/m^2")
    print(f"Error  : {rel_err*100:.4f}%")
    print("-------------------------------------------------------")

    # Pass criterion: relative error < 0.1% (1e-3)
    if rel_err < 1e-3:
        print("\nSTATUS: FIXED HEAT FLUX CONSERVATION TEST PASSED!")
        print("Global domain energy balance perfectly matches integrated Neumann boundary heat flux.")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
