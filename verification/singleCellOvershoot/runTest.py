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

def setup_single_cell_case(case_dir, T0, Q_source, L_heat, T_lm=300.0, T_um=310.0, T_lf=295.0, T_uf=305.0):
    os.makedirs(case_dir, exist_ok=True)
    os.chdir(case_dir)
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/polyMesh")
    
    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices ((0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 1) (1 0 1) (1 1 1) (0 1 1));
blocks ( hex (0 1 2 3 4 5 6 7) pcm (1 1 1) simpleGrading (1 1 1) );
edges ();
boundary ( walls { type empty; faces ((0 1 5 4) (3 2 6 7) (0 3 7 4) (1 2 6 5) (0 1 2 3) (4 5 6 7)); } );
""")

    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 100; deltaT 100;
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
PIMPLE { nNonOrthogonalCorrectors 0; nNonLinearCorrectors 50; nonLinearTolerance 1e-6; }
solvers { "h.*" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; } }
""")

    with open("system/fvOptions", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvOptions; }}
heatSource
{{
    type scalarSemiImplicitSource;
    active true;
    selectionMode all;
    volumeMode absolute;
    injectionRateSuSp
    {{
        h ({Q_source} 0);
    }}
}}
""")

    os.makedirs("system/pcm", exist_ok=True)
    with open("system/pcm/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/pcm/fvSolution", "w") as f: f.write(open("system/fvSolution").read())
    with open("system/pcm/fvOptions", "w") as f: f.write(open("system/fvOptions").read())

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
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }}
active true;
phaseChange
{{
    active true; phaseChangeMode EHC;
    melting {{ T_lowerBound {T_lm}; T_upperBound {T_um}; latentHeat {L_heat}; }}
    freezing {{ T_lowerBound {T_lf}; T_upperBound {T_uf}; latentHeat {L_heat}; }}
    hysteresis {{ active true; historyLength 5; directionDetection {{ method deltaT; tolerance 1e-6; }} }}
    density {{ model linear; rhoRef 1000.0; rhoSolid 1000.0; rhoLiquid 1000.0; }}
    thermophysical {{ mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 1.0; kLiquid 1.0; }}
}}
""")

    h0 = 2000.0 * T0
    alpha0 = 0.0 if T0 <= T_lm else (1.0 if T0 >= T_um else (T0 - T_lm) / (T_um - T_lm))
    traj0 = 1.0 if Q_source >= 0 else 0.0

    os.makedirs("0/pcm", exist_ok=True)
    with open("0/pcm/T", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }}
dimensions [0 0 0 1 0 0 0]; internalField uniform {T0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/pcm/liquidFraction", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object liquidFraction; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {alpha0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/pcm/heatingTrajectory", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object heatingTrajectory; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {traj0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/pcm/h", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }}
dimensions [0 2 -2 0 0 0 0]; internalField uniform {h0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { walls { type empty; } }
""")

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()

    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/pcm/fvSolution 2>/dev/null || true")
    run_cmd("cp system/fvOptions system/pcm/fvOptions 2>/dev/null || true")

    output = run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_sim = parse_openfoam_field("100/pcm/T")[0]
    alpha_sim = parse_openfoam_field("100/pcm/liquidFraction")[0]

    return T_sim, alpha_sim, output

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== Comprehensive Single-Cell Phase Change Verification Suite in {base_dir} ===")

    all_passed = True

    # --- Case 1: Single-Step Heating Overshoot ---
    # T0 = 280 K, Q = 2.4 MW, dt = 100 s => E_in = 2.4e8 J => delta_H_mass = 240 kJ/kg
    # Sensible heat 280->300 = 20*2000 = 40 kJ/kg
    # Latent heat = 100 kJ/kg (300->310 K)
    # Remaining sensible heat = 240 - 40 - 100 = 100 kJ/kg => Delta T = 100000/2000 = 50 K => T_exact = 350.0 K
    print("\n--- Case 1: Heating Jump Overshoot (T0 = 280 K -> 350 K in 1 step) ---")
    case1_dir = os.path.join(base_dir, "case1_heating")
    T1, a1, log1 = setup_single_cell_case(case1_dir, T0=280.0, Q_source=2400000.0, L_heat=100000.0)
    err1 = abs(T1 - 350.0)
    print(f"Simulated T = {T1:.6f} K, alphaL = {a1:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err1:.6f} K")
    if err1 < 0.001 and abs(a1 - 1.0) < 0.001:
        print("STATUS: CASE 1 PASSED!")
    else:
        print("STATUS: CASE 1 FAILED!")
        all_passed = False

    # --- Case 2: Single-Step Cooling Overshoot ---
    # T0 = 350 K, Q = -2.2 MW, dt = 100 s => E_out = -2.2e8 J/m^3 => delta_H_mass = -220 kJ/kg
    # Sensible heat 350->305 = 45*2000 = 90 kJ/kg
    # Latent heat = 100 kJ/kg (305->295 K)
    # Remaining sensible heat = 220 - 90 - 100 = 30 kJ/kg => Delta T = -30000/2000 = -15 K => T_exact = 295 - 15 = 280.0 K for exact piecewise, but 290.0 K for Backward-Euler secant formulation over single step.
    # Enthalpy balance: 1000 * 2000 * (T_sim - 350) + 1000 * 100000 * (0 - 1) = -2.2e8 J/m^3 => T_exact = 290.000000 K
    print("\n--- Case 2: Cooling Jump Overshoot (T0 = 350 K -> 290 K in 1 step) ---")
    case2_dir = os.path.join(base_dir, "case2_cooling")
    T2, a2, log2 = setup_single_cell_case(case2_dir, T0=350.0, Q_source=-2200000.0, L_heat=100000.0)
    err2 = abs(T2 - 290.0)
    print(f"Simulated T = {T2:.6f} K, alphaL = {a2:.6f}")
    print(f"Exact T     = 290.000000 K, alphaL = 0.000000")
    print(f"Temperature Error = {err2:.6f} K")
    if err2 < 0.001 and abs(a2 - 0.0) < 0.001:
        print("STATUS: CASE 2 PASSED!")
    else:
        print("STATUS: CASE 2 FAILED!")
        all_passed = False

    # --- Case 3: Partial Window Crossing ---
    # T0 = 295 K, Q = 0.7 MW, dt = 100 s => E_in = 7.0e7 J/m^3 => delta_H_mass = 70 kJ/kg
    # Sensible heat 295->300 = 5*2000 = 10 kJ/kg
    # Latent heat 0->0.5 = 50 kJ/kg (300->305 K)
    # Sensible heat 300->305 = 5*2000 = 10 kJ/kg
    # Total enthalpy = 10 + 50 + 10 = 70 kJ/kg => T_exact = 305.000000 K, alphaL_exact = 0.500000
    print("\n--- Case 3: Partial Window Crossing (T0 = 295 K -> 305 K, alphaL = 0.50) ---")
    case3_dir = os.path.join(base_dir, "case3_partial")
    T3, a3, log3 = setup_single_cell_case(case3_dir, T0=295.0, Q_source=700000.0, L_heat=100000.0)
    err3 = abs(T3 - 305.0)
    err_a3 = abs(a3 - 0.50)
    print(f"Simulated T = {T3:.6f} K, alphaL = {a3:.6f}")
    print(f"Exact T     = 305.000000 K, alphaL = 0.500000")
    print(f"Temperature Error = {err3:.6f} K, Alpha Error = {err_a3:.6f}")
    if err3 < 0.001 and err_a3 < 0.001:
        print("STATUS: CASE 3 PASSED!")
    else:
        print("STATUS: CASE 3 FAILED!")
        all_passed = False

    # --- Case 4: Realistic Latent Heat (L = 163 kJ/kg) ---
    # T0 = 280 K, Q = 3.03 MW, dt = 100 s => E_in = 3.03e8 J => delta_H_mass = 303 kJ/kg
    # Sensible heat 280->300 = 20*2000 = 40 kJ/kg
    # Latent heat = 163 kJ/kg (300->310 K)
    # Remaining sensible heat = 303 - 40 - 163 = 100 kJ/kg => Delta T = 100000/2000 = 50 K => T_exact = 350.0 K
    print("\n--- Case 4: Realistic Latent Heat (L = 163 kJ/kg, T0 = 280 K -> 350 K) ---")
    case4_dir = os.path.join(base_dir, "case4_realisticL")
    T4, a4, log4 = setup_single_cell_case(case4_dir, T0=280.0, Q_source=3030000.0, L_heat=163000.0)
    err4 = abs(T4 - 350.0)
    print(f"Simulated T = {T4:.6f} K, alphaL = {a4:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err4:.6f} K")
    if err4 < 0.001 and abs(a4 - 1.0) < 0.001:
        print("STATUS: CASE 4 PASSED!")
    else:
        print("STATUS: CASE 4 FAILED!")
        all_passed = False

    print("\n=======================================================")
    print("      SINGLE-CELL VERIFICATION SUITE SUMMARY           ")
    print("=======================================================")
    print(f"Case 1 (Heating Jump)   : {'PASSED' if err1 < 0.001 else 'FAILED'} (err = {err1:.6f} K)")
    print(f"Case 2 (Cooling Jump)   : {'PASSED' if err2 < 0.001 else 'FAILED'} (err = {err2:.6f} K)")
    print(f"Case 3 (Partial Melt)   : {'PASSED' if err3 < 0.001 else 'FAILED'} (err = {err3:.6f} K)")
    print(f"Case 4 (Realistic L)    : {'PASSED' if err4 < 0.001 else 'FAILED'} (err = {err4:.6f} K)")
    print("-------------------------------------------------------")

    if all_passed:
        print("\nALL SINGLE-CELL VERIFICATION TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSOME TESTS FAILED!")
        sys.exit(1)

if __name__ == "__main__":
    main()
