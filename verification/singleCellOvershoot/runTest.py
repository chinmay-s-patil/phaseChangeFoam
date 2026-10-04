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

def setup_single_cell_case(case_dir, T0, Q_source, L_heat, T_lm=300.0, T_um=310.0, T_lf=295.0, T_uf=305.0, alpha0=None, traj0=None, Cps=2000.0, Cpl=2000.0, rhoS=1000.0, rhoL=1000.0):
    os.makedirs(case_dir, exist_ok=True)
    os.chdir(case_dir)
    run_cmd("rm -rf [1-9]* 0.* constant/phaseChange/polyMesh constant/polyMesh")
    
    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices ((0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 1) (1 0 1) (1 1 1) (0 1 1));
blocks ( hex (0 1 2 3 4 5 6 7) phaseChange (1 1 1) simpleGrading (1 1 1) );
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

    os.makedirs("system/phaseChange", exist_ok=True)
    with open("system/phaseChange/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/phaseChange/fvSolution", "w") as f: f.write(open("system/fvSolution").read())
    with open("system/phaseChange/fvOptions", "w") as f: f.write(open("system/fvOptions").read())

    os.makedirs("constant/phaseChange", exist_ok=True)
    with open("constant/regionProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid () solid (phaseChange) porousFluid () porousSolid () );
""")
    with open("constant/g", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0]; value (0 0 0);
""")
    with open("constant/phaseChange/radiationProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/phaseChange"; object radiationProperties; }
radiationModel none;
""")
    with open("constant/phaseChange/thermophysicalProperties", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/phaseChange"; object thermophysicalProperties; }}
thermoType {{ type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }}
mixture {{ specie {{ molWeight 100; }} transport {{ kappa 1.0; }} thermodynamics {{ Cp {Cps}; Hf 0; }} equationOfState {{ rho {rhoS}; }} }}
""")
    with open("constant/phaseChange/phaseChangeDict", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/phaseChange"; object phaseChangeDict; }}
active true;
phaseChange
{{
    active true; phaseChangeMode EHC;
    melting {{ T_lowerBound {T_lm}; T_upperBound {T_um}; latentHeat {L_heat}; }}
    freezing {{ T_lowerBound {T_lf}; T_upperBound {T_uf}; latentHeat {L_heat}; }}
    hysteresis {{ active true; }}
    density {{ model linear; rhoRef {rhoS}; rhoSolid {rhoS}; rhoLiquid {rhoL}; allowNonConservativeDensity true; }}
    thermophysical {{ mode custom; CpSolid {Cps}; CpLiquid {Cpl}; kSolid 1.0; kLiquid 1.0; }}
}}
""")

    if alpha0 is None:
        alpha0 = 0.0 if T0 <= T_lm else (1.0 if T0 >= T_um else (T0 - T_lm) / (T_um - T_lm))
    if traj0 is None:
        traj0 = 1.0 if Q_source >= 0 else 0.0

    cp_mush = 0.5 * (Cps + Cpl)
    h_lm = Cps * T_lm
    h_um = h_lm + cp_mush * (T_um - T_lm)
    if T0 <= T_lm:
        h0 = Cps * T0
    elif T0 >= T_um:
        h0 = h_um + Cpl * (T0 - T_um)
    else:
        h0 = h_lm + cp_mush * (T0 - T_lm)

    os.makedirs("0/phaseChange", exist_ok=True)
    with open("0/phaseChange/T", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/phaseChange"; object T; }}
dimensions [0 0 0 1 0 0 0]; internalField uniform {T0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/phaseChange/phaseFraction", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/phaseChange"; object phaseFraction; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {alpha0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/phaseChange/heatingTrajectory", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/phaseChange"; object heatingTrajectory; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {traj0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/phaseChange/h", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/phaseChange"; object h; }}
dimensions [0 2 -2 0 0 0 0]; internalField uniform {h0};
boundaryField {{ walls {{ type empty; }} }}
""")
    with open("0/phaseChange/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/phaseChange"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { walls { type empty; } }
""")

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()

    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/phaseChange")
    run_cmd("cp -r constant/polyMesh constant/phaseChange/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/phaseChange/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/phaseChange/fvSolution 2>/dev/null || true")
    run_cmd("cp system/fvOptions system/phaseChange/fvOptions 2>/dev/null || true")

    output = run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_sim = parse_openfoam_field("100/phaseChange/T")[0]
    alpha_sim = parse_openfoam_field("100/phaseChange/phaseFraction")[0]

    return T_sim, alpha_sim, output

def calc_h_sens(T, traj, Cps, Cpl, Tlm=300.0, Tum=310.0, Tlf=295.0, Tuf=305.0):
    Tl = Tlm if traj > 0.5 else Tlf
    Tu = Tum if traj > 0.5 else Tuf
    if T <= Tl:
        return Cps * T
    elif T >= Tu:
        h_l = Cps * Tl
        h_mush = 0.5 * (Cps + Cpl) * (Tu - Tl)
        return h_l + h_mush + Cpl * (T - Tu)
    else:
        h_l = Cps * Tl
        dT = T - Tl
        dTum = Tu - Tl
        if dTum > 1e-8:
            return h_l + Cps * dT + 0.5 * (Cpl - Cps) * (dT**2) / dTum
        else:
            return h_l

def calc_delta_h_sens(T0, T1, alpha0, traj0, traj1, Cps, Cpl, Tlm=300.0, Tum=310.0, Tlf=295.0, Tuf=305.0):
    dT = T1 - T0
    if abs(dT) < 1e-12:
        return 0.0
    
    if True:
        if traj1 > 0.5: # Heating branch
            Texit = Tlm + alpha0 * (Tum - Tlm)
            if Texit < Tlm: Texit = Tlm
            if Texit > Tum: Texit = Tum
            if T0 < Texit:
                if T1 <= Texit:
                    Cp_plateau = (1.0 - alpha0) * Cps + alpha0 * Cpl
                    return Cp_plateau * dT
                else:
                    Cp_plateau = (1.0 - alpha0) * Cps + alpha0 * Cpl
                    dH1 = Cp_plateau * (Texit - T0)
                    dH2 = calc_h_sens(T1, 1.0, Cps, Cpl, Tlm, Tum, Tlf, Tuf) - calc_h_sens(Texit, 1.0, Cps, Cpl, Tlm, Tum, Tlf, Tuf)
                    return dH1 + dH2
        else: # Cooling branch
            Texit = Tlf + alpha0 * (Tuf - Tlf)
            if Texit < Tlf: Texit = Tlf
            if Texit > Tuf: Texit = Tuf
            if T0 > Texit:
                if T1 >= Texit:
                    Cp_plateau = (1.0 - alpha0) * Cps + alpha0 * Cpl
                    return Cp_plateau * dT
                else:
                    Cp_plateau = (1.0 - alpha0) * Cps + alpha0 * Cpl
                    dH1 = Cp_plateau * (Texit - T0)
                    dH2 = calc_h_sens(T1, 0.0, Cps, Cpl, Tlm, Tum, Tlf, Tuf) - calc_h_sens(Texit, 0.0, Cps, Cpl, Tlm, Tum, Tlf, Tuf)
                    return dH1 + dH2

    return calc_h_sens(T1, traj1, Cps, Cpl, Tlm, Tum, Tlf, Tuf) - calc_h_sens(T0, traj1, Cps, Cpl, Tlm, Tum, Tlf, Tuf)

def calc_exact_q_source(T0, T1, alpha0, alpha1, traj0, traj1, Cps, Cpl, rhoS, rhoL, L, dt=100.0, Tlm=300.0, Tum=310.0, Tlf=295.0, Tuf=305.0):
    delta_h_sens = calc_delta_h_sens(T0, T1, alpha0, traj0, traj1, Cps, Cpl, Tlm, Tum, Tlf, Tuf)
    rho_eff1 = (1.0 - alpha1) * rhoS + alpha1 * rhoL
    rho_bar = 0.5 * (rhoS + rhoL)
    Q = (rho_eff1 * delta_h_sens + rho_bar * L * (alpha1 - alpha0)) / dt
    return Q

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== Comprehensive Single-Cell Phase Change Verification Suite in {base_dir} ===")

    all_passed = True

    # --- Case 1: Single-Step Heating Overshoot ---
    print("\n--- Case 1: Heating Jump Overshoot (T0 = 280 K -> 350 K in 1 step) ---")
    case1_dir = os.path.join(base_dir, "case1_heating")
    T1, a1, log1 = setup_single_cell_case(case1_dir, T0=280.0, Q_source=2400000.0, L_heat=100000.0)
    err1 = abs(T1 - 350.0)
    print(f"Simulated T = {T1:.6f} K, alphaL = {a1:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err1:.6f} K")
    if err1 < 0.01 and abs(a1 - 1.0) < 0.001:
        print("STATUS: CASE 1 PASSED!")
    else:
        print("STATUS: CASE 1 FAILED!")
        all_passed = False

    # --- Case 2: Single-Step Cooling Overshoot ---
    print("\n--- Case 2: Cooling Jump Overshoot (T0 = 350 K -> 290 K in 1 step) ---")
    case2_dir = os.path.join(base_dir, "case2_cooling")
    T2, a2, log2 = setup_single_cell_case(case2_dir, T0=350.0, Q_source=-2200000.0, L_heat=100000.0)
    err2 = abs(T2 - 290.0)
    print(f"Simulated T = {T2:.6f} K, alphaL = {a2:.6f}")
    print(f"Exact T     = 290.000000 K, alphaL = 0.000000")
    print(f"Temperature Error = {err2:.6f} K")
    if err2 < 0.01 and abs(a2 - 0.0) < 0.001:
        print("STATUS: CASE 2 PASSED!")
    else:
        print("STATUS: CASE 2 FAILED!")
        all_passed = False

    # --- Case 3: Partial Window Crossing ---
    print("\n--- Case 3: Partial Window Crossing (T0 = 295 K -> 305 K, alphaL = 0.50) ---")
    case3_dir = os.path.join(base_dir, "case3_partial")
    T3, a3, log3 = setup_single_cell_case(case3_dir, T0=295.0, Q_source=700000.0, L_heat=100000.0)
    err3 = abs(T3 - 305.0)
    err_a3 = abs(a3 - 0.50)
    print(f"Simulated T = {T3:.6f} K, alphaL = {a3:.6f}")
    print(f"Exact T     = 305.000000 K, alphaL = 0.500000")
    print(f"Temperature Error = {err3:.6f} K, Alpha Error = {err_a3:.6f}")
    if err3 < 0.01 and err_a3 < 0.001:
        print("STATUS: CASE 3 PASSED!")
    else:
        print("STATUS: CASE 3 FAILED!")
        all_passed = False

    # --- Case 4: Realistic Latent Heat (L = 163 kJ/kg) ---
    print("\n--- Case 4: Realistic Latent Heat (L = 163 kJ/kg, T0 = 280 K -> 350 K) ---")
    case4_dir = os.path.join(base_dir, "case4_realisticL")
    T4, a4, log4 = setup_single_cell_case(case4_dir, T0=280.0, Q_source=3030000.0, L_heat=163000.0)
    err4 = abs(T4 - 350.0)
    print(f"Simulated T = {T4:.6f} K, alphaL = {a4:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err4:.6f} K")
    if err4 < 0.30 and abs(a4 - 1.0) < 0.001:
        print("STATUS: CASE 4 PASSED!")
    else:
        print("STATUS: CASE 4 FAILED!")
        all_passed = False

    # --- Case 5: Mushy Start State & Partial Melt Heating (T0 = 305 K, alpha0 = 0.5 -> T = 308 K, alphaL = 0.80) ---
    print("\n--- Case 5: Mushy Start State & Partial Melt Heating (T0 = 305 K, alpha0 = 0.5 -> T = 308 K, alphaL = 0.80) ---")
    case5_dir = os.path.join(base_dir, "case5_mushyStart")
    T5, a5, log5 = setup_single_cell_case(case5_dir, T0=305.0, Q_source=360000.0, L_heat=100000.0, alpha0=0.5, traj0=1.0)
    err5 = abs(T5 - 308.0)
    err_a5 = abs(a5 - 0.80)
    print(f"Simulated T = {T5:.6f} K, alphaL = {a5:.6f}")
    print(f"Exact T     = 308.000000 K, alphaL = 0.800000")
    print(f"Temperature Error = {err5:.6f} K, Alpha Error = {err_a5:.6f}")
    if err5 < 0.001 and err_a5 < 0.001:
        print("STATUS: CASE 5 PASSED!")
    else:
        print("STATUS: CASE 5 FAILED!")
        all_passed = False

    # --- Case 6: Cooling Reversal Plateau (T0 = 308 K, alpha0 = 0.8 -> T = 306 K, alphaL = 0.80) ---
    print("\n--- Case 6: Cooling Reversal Plateau (T0 = 308 K, alpha0 = 0.8 -> T = 306 K, alphaL = 0.80) ---")
    case6_dir = os.path.join(base_dir, "case6_coolingReversal")
    Q6 = calc_exact_q_source(T0=308.0, T1=306.0, alpha0=0.8, alpha1=0.8, traj0=1.0, traj1=0.0, Cps=2000.0, Cpl=2000.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T6, a6, log6 = setup_single_cell_case(case6_dir, T0=308.0, Q_source=Q6, L_heat=100000.0, alpha0=0.8, traj0=0.0)
    err6 = abs(T6 - 306.0)
    err_a6 = abs(a6 - 0.80)
    print(f"Simulated T = {T6:.6f} K, alphaL = {a6:.6f}")
    print(f"Exact T     = 306.000000 K, alphaL = 0.800000")
    print(f"Temperature Error = {err6:.6f} K, Alpha Error = {err_a6:.6f}")
    if err6 < 0.001 and err_a6 < 0.001:
        print("STATUS: CASE 6 PASSED!")
    else:
        print("STATUS: CASE 6 FAILED!")
        all_passed = False

    # --- Case 7: Heating Reversal Plateau (T0 = 297 K, alpha0 = 0.2 -> T = 298 K, alphaL = 0.20) ---
    print("\n--- Case 7: Heating Reversal Plateau (T0 = 297 K, alpha0 = 0.2 -> T = 298 K, alphaL = 0.20) ---")
    case7_dir = os.path.join(base_dir, "case7_heatingReversal")
    Q7 = calc_exact_q_source(T0=297.0, T1=298.0, alpha0=0.2, alpha1=0.2, traj0=0.0, traj1=1.0, Cps=2000.0, Cpl=2000.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T7, a7, log7 = setup_single_cell_case(case7_dir, T0=297.0, Q_source=Q7, L_heat=100000.0, alpha0=0.2, traj0=1.0)
    err7 = abs(T7 - 298.0)
    err_a7 = abs(a7 - 0.20)
    print(f"Simulated T = {T7:.6f} K, alphaL = {a7:.6f}")
    print(f"Exact T     = 298.000000 K, alphaL = 0.200000")
    print(f"Temperature Error = {err7:.6f} K, Alpha Error = {err_a7:.6f}")
    if err7 < 0.001 and err_a7 < 0.001:
        print("STATUS: CASE 7 PASSED!")
    else:
        print("STATUS: CASE 7 FAILED!")
        all_passed = False

    # --- Case 8: Exact Base Cp (Cps = 1980, Cpl = 2320 J/(kg.K), T0 = 280 K -> 350 K) ---
    print("\n--- Case 8: Exact Piecewise Path Integrated Base Cp (Cps = 1980, Cpl = 2320 J/(kg.K), T0 = 280 K -> 350 K) ---")
    case8_dir = os.path.join(base_dir, "case8_unequalCp")
    Q8 = calc_exact_q_source(T0=280.0, T1=350.0, alpha0=0.0, alpha1=1.0, traj0=1.0, traj1=1.0, Cps=1980.0, Cpl=2320.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T8, a8, log8 = setup_single_cell_case(case8_dir, T0=280.0, Q_source=Q8, L_heat=100000.0, Cps=1980.0, Cpl=2320.0)
    err8 = abs(T8 - 350.0)
    print(f"Simulated T = {T8:.6f} K, alphaL = {a8:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err8:.6f} K")
    if err8 < 0.01 and abs(a8 - 1.0) < 0.001:
        print("STATUS: CASE 8 PASSED!")
    else:
        print("STATUS: CASE 8 FAILED!")
        all_passed = False

    # --- Case 9: Variable Density (rhoS = 1967, rhoL = 1850 kg/m^3, T0 = 280 K -> 350 K) ---
    print("\n--- Case 9: Variable Density (rhoS = 1967, rhoL = 1850 kg/m^3, T0 = 280 K -> 350 K) ---")
    case9_dir = os.path.join(base_dir, "case9_variableDensity")
    Q9 = calc_exact_q_source(T0=280.0, T1=350.0, alpha0=0.0, alpha1=1.0, traj0=1.0, traj1=1.0, Cps=1980.0, Cpl=2320.0, rhoS=1967.0, rhoL=1850.0, L=100000.0)
    T9, a9, log9 = setup_single_cell_case(case9_dir, T0=280.0, Q_source=Q9, L_heat=100000.0, Cps=1980.0, Cpl=2320.0, rhoS=1967.0, rhoL=1850.0)
    err9 = abs(T9 - 350.0)
    print(f"Simulated T = {T9:.6f} K, alphaL = {a9:.6f}")
    print(f"Exact T     = 350.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err9:.6f} K")
    if err9 < 0.01 and abs(a9 - 1.0) < 0.001:
        print("STATUS: CASE 9 PASSED!")
    else:
        print("STATUS: CASE 9 FAILED!")
        all_passed = False

    # --- Case 10: Cooling Jump with Unequal Cp (Cps = 1980, Cpl = 2320 J/(kg.K), T0 = 350 K -> 290 K) ---
    print("\n--- Case 10: Cooling Jump with Unequal Cp (Cps = 1980, Cpl = 2320 J/(kg.K), T0 = 350 K -> 290 K) ---")
    case10_dir = os.path.join(base_dir, "case10_unequalCpCooling")
    Q10 = calc_exact_q_source(T0=350.0, T1=290.0, alpha0=1.0, alpha1=0.0, traj0=0.0, traj1=0.0, Cps=1980.0, Cpl=2320.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T10, a10, log10 = setup_single_cell_case(case10_dir, T0=350.0, Q_source=Q10, L_heat=100000.0, Cps=1980.0, Cpl=2320.0, alpha0=1.0, traj0=0.0)
    err10 = abs(T10 - 290.0)
    print(f"Simulated T = {T10:.6f} K, alphaL = {a10:.6f}")
    print(f"Exact T     = 290.000000 K, alphaL = 0.000000")
    print(f"Temperature Error = {err10:.6f} K")
    if err10 < 0.01 and abs(a10 - 0.0) < 0.001:
        print("STATUS: CASE 10 PASSED!")
    else:
        print("STATUS: CASE 10 FAILED!")
        all_passed = False

    # --- Case 11: Cooling Reversal Plateau with Unequal Cp ---
    print("\n--- Case 11: Cooling Reversal Plateau with Unequal Cp (T0 = 308 K, alpha0 = 0.8 -> T = 306 K, alphaL = 0.80) ---")
    case11_dir = os.path.join(base_dir, "case11_unequalCpReversal")
    Q11 = calc_exact_q_source(T0=308.0, T1=306.0, alpha0=0.8, alpha1=0.8, traj0=1.0, traj1=0.0, Cps=1980.0, Cpl=2320.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T11, a11, log11 = setup_single_cell_case(case11_dir, T0=308.0, Q_source=Q11, L_heat=100000.0, Cps=1980.0, Cpl=2320.0, alpha0=0.8, traj0=1.0)
    err11 = abs(T11 - 306.0)
    err_a11 = abs(a11 - 0.80)
    print(f"Simulated T = {T11:.6f} K, alphaL = {a11:.6f}")
    print(f"Exact T     = 306.000000 K, alphaL = 0.800000")
    print(f"Temperature Error = {err11:.6f} K, Alpha Error = {err_a11:.6f}")
    if err11 < 0.001 and err_a11 < 0.001:
        print("STATUS: CASE 11 PASSED!")
    else:
        print("STATUS: CASE 11 FAILED!")
        all_passed = False

    # --- Case 12: Closed Thermal Cycle Net Enthalpy Conservation ---
    print("\n--- Case 12: Closed Thermal Cycle (280 K -> 350 K -> 280 K) Enthalpy Conservation ---")
    case12_dir = os.path.join(base_dir, "case12_closedCycle")
    T12_step1, a12_step1, log12_step1 = setup_single_cell_case(case12_dir, T0=280.0, Q_source=2400000.0, L_heat=100000.0, Cps=2000.0, Cpl=2000.0)
    with open(os.path.join(case12_dir, "system/phaseChange/fvOptions"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system/phaseChange"; object fvOptions; }
heatSource { type scalarSemiImplicitSource; active true; selectionMode all; volumeMode specific;
    sources { h ( -2400000 0 ); } }
""")
    with open(os.path.join(case12_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam; startFrom latestTime; startTime 100; stopAt endTime; endTime 200; deltaT 100; writeControl timeStep; writeInterval 1;
""")
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()
    run_cmd(f"mkdir -p {case12_dir}/100/phaseChange/polyMesh && cp -r {case12_dir}/constant/phaseChange/polyMesh/* {case12_dir}/100/phaseChange/polyMesh/ 2>/dev/null || true")
    run_cmd(f"cd {case12_dir} && bash -c '{of_env}; {solver_bin}'")
    T12 = parse_openfoam_field(os.path.join(case12_dir, "200/phaseChange/T"))[0]
    a12 = parse_openfoam_field(os.path.join(case12_dir, "200/phaseChange/phaseFraction"))[0]
    err12 = abs(T12 - 280.0)
    err_a12 = abs(a12 - 0.0)
    print(f"Simulated T = {T12:.6f} K, alphaL = {a12:.6f}")
    print(f"Exact T     = 280.000000 K, alphaL = 0.000000")
    print(f"Temperature Error = {err12:.6f} K, Alpha Error = {err_a12:.6f}")
    # --- Case 13: Continuous Multi-Step Cooling Reversal with Unequal Cp (308 -> 306 -> 304 -> 302 K) ---
    print("\n--- Case 13: Continuous Multi-Step Cooling Reversal with Unequal Cp (308 -> 306 -> 304 -> 302 K) ---")
    case13_dir = os.path.join(base_dir, "case13_multiStepReversal")
    Q13a = calc_exact_q_source(308.0, 306.0, 0.8, 0.8, 0.0, 0.0, 1980.0, 2320.0, 1000.0, 1000.0, 100000.0)
    Q13b = calc_exact_q_source(306.0, 304.0, 0.8, 0.8, 0.0, 0.0, 1980.0, 2320.0, 1000.0, 1000.0, 100000.0)
    Q13c = calc_exact_q_source(304.0, 302.0, 0.8, 0.7, 0.0, 0.0, 1980.0, 2320.0, 1000.0, 1000.0, 100000.0)

    setup_single_cell_case(case13_dir, T0=308.0, Q_source=Q13a, L_heat=100000.0, Cps=1980.0, Cpl=2320.0, alpha0=0.8, traj0=0.0)

    with open(f"{case13_dir}/system/controlDict", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 300; deltaT 100;
writeControl runTime; writeInterval 100; purgeWrite 0; writeFormat ascii;
""")
    fvOpt13 = f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvOptions; }}
heatSource
{{
    type scalarSemiImplicitSource;
    active true;
    selectionMode all;
    volumeMode absolute;
    injectionRateSuSp
    {{
        h
        {{
            Su table
            (
                (0 {Q13a})
                (100 {Q13a})
                (100.001 {Q13b})
                (200 {Q13b})
                (200.001 {Q13c})
                (300 {Q13c})
            );
            Sp 0;
        }}
    }}
}}
"""
    with open(f"{case13_dir}/system/fvOptions", "w") as f: f.write(fvOpt13)
    with open(f"{case13_dir}/system/phaseChange/fvOptions", "w") as f: f.write(fvOpt13)

    run_cmd(f"cd {case13_dir} && bash -c '{of_env}; {solver_bin}'")

    T13a = parse_openfoam_field(os.path.join(case13_dir, "100/phaseChange/T"))[0]
    a13a = parse_openfoam_field(os.path.join(case13_dir, "100/phaseChange/phaseFraction"))[0]
    err13a = abs(T13a - 306.0)

    T13b = parse_openfoam_field(os.path.join(case13_dir, "200/phaseChange/T"))[0]
    a13b = parse_openfoam_field(os.path.join(case13_dir, "200/phaseChange/phaseFraction"))[0]
    err13b = abs(T13b - 304.0)

    T13c = parse_openfoam_field(os.path.join(case13_dir, "300/phaseChange/T"))[0]
    a13c = parse_openfoam_field(os.path.join(case13_dir, "300/phaseChange/phaseFraction"))[0]
    err13c = abs(T13c - 302.0)

    err13 = max(err13a, err13b, err13c)

    print(f"Step 1 (308->306 K): Simulated T = {T13a:.6f} K, alpha = {a13a:.6f}")
    print(f"Step 2 (306->304 K): Simulated T = {T13b:.6f} K, alpha = {a13b:.6f}")
    print(f"Step 3 (304->302 K): Simulated T = {T13c:.6f} K, alpha = {a13c:.6f}")
    if err13 < 0.001 and abs(a13a - 0.8) < 0.001 and abs(a13b - 0.8) < 0.001 and abs(a13c - 0.7) < 0.001:
        print("STATUS: CASE 13 PASSED!")
    else:
        print("STATUS: CASE 13 FAILED!")
        all_passed = False

    # --- Case 14: Closed Thermal Cycle with Unequal Cp & Hysteresis (280 -> 350 -> 280 K) ---
    print("\n--- Case 14: Closed Thermal Cycle with Unequal Cp (280 -> 350 -> 280 K) ---")
    case14_dir = os.path.join(base_dir, "case14_unequalCpClosedCycle")
    Q14_heat = calc_exact_q_source(280.0, 350.0, 0.0, 1.0, 1.0, 1.0, 1980.0, 2320.0, 1000.0, 1000.0, 100000.0)
    T14a, a14a, _ = setup_single_cell_case(case14_dir, T0=280.0, Q_source=Q14_heat, L_heat=100000.0, Cps=1980.0, Cpl=2320.0)
    
    Q14_cool = calc_exact_q_source(350.0, 280.0, 1.0, 0.0, 0.0, 0.0, 1980.0, 2320.0, 1000.0, 1000.0, 100000.0)
    with open(f"{case14_dir}/system/controlDict", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 100; stopAt endTime; endTime 200; deltaT 100;
writeControl runTime; writeInterval 100; purgeWrite 0; writeFormat ascii;
""")
    fvOpt14_cool = f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvOptions; }}
heatSource {{ type scalarSemiImplicitSource; active true; selectionMode all; volumeMode absolute; injectionRateSuSp {{ h ({Q14_cool} 0); }} }}
"""
    with open(f"{case14_dir}/system/fvOptions", "w") as f: f.write(fvOpt14_cool)
    with open(f"{case14_dir}/system/phaseChange/fvOptions", "w") as f: f.write(fvOpt14_cool)
    run_cmd(f"mkdir -p {case14_dir}/100/phaseChange/polyMesh && cp -r {case14_dir}/constant/phaseChange/polyMesh/* {case14_dir}/100/phaseChange/polyMesh/ 2>/dev/null || true")
    run_cmd(f"cd {case14_dir} && bash -c '{of_env}; {solver_bin}'")
    T14 = parse_openfoam_field(os.path.join(case14_dir, "200/phaseChange/T"))[0]
    a14 = parse_openfoam_field(os.path.join(case14_dir, "200/phaseChange/phaseFraction"))[0]
    err14 = abs(T14 - 280.0)
    err_a14 = abs(a14 - 0.0)
    print(f"Simulated T = {T14:.6f} K, alphaL = {a14:.6f}")
    print(f"Exact T     = 280.000000 K, alphaL = 0.000000")
    print(f"Temperature Error = {err14:.6f} K, Alpha Error = {err_a14:.6f}")
    if err14 < 0.01 and err_a14 < 0.001:
        print("STATUS: CASE 14 PASSED!")
    else:
        print("STATUS: CASE 14 FAILED!")
        all_passed = False

    # --- Case 15: Liquid Reversal Heating Plateau (T0 = 309 K, alpha0 = 1.0 -> 312 K) ---
    print("\n--- Case 15: Liquid Reversal Heating Plateau (T0 = 309 K, alpha0 = 1.0 -> 312 K) ---")
    case15_dir = os.path.join(base_dir, "case15_liquidReversal")
    Q15 = calc_exact_q_source(T0=309.0, T1=312.0, alpha0=1.0, alpha1=1.0, traj0=0.0, traj1=1.0, Cps=1980.0, Cpl=2320.0, rhoS=1000.0, rhoL=1000.0, L=100000.0)
    T15, a15, log15 = setup_single_cell_case(case15_dir, T0=309.0, Q_source=Q15, L_heat=100000.0, Cps=1980.0, Cpl=2320.0, alpha0=1.0, traj0=0.0)
    err15 = abs(T15 - 312.0)
    err_a15 = abs(a15 - 1.0)
    print(f"Simulated T = {T15:.6f} K, alphaL = {a15:.6f}")
    print(f"Exact T     = 312.000000 K, alphaL = 1.000000")
    print(f"Temperature Error = {err15:.6f} K, Alpha Error = {err_a15:.6f}")
    if err15 < 0.01 and err_a15 < 0.001:
        print("STATUS: CASE 15 PASSED!")
    else:
        print("STATUS: CASE 15 FAILED!")
        all_passed = False

    print("\n=======================================================")
    print("      SINGLE-CELL VERIFICATION SUITE SUMMARY           ")
    print("=======================================================")
    print(f"Case 1 (Heating Jump)          : {'PASSED' if err1 < 0.01 else 'FAILED'} (err = {err1:.6f} K)")
    print(f"Case 2 (Cooling Jump)          : {'PASSED' if err2 < 0.01 else 'FAILED'} (err = {err2:.6f} K)")
    print(f"Case 3 (Partial Melt)          : {'PASSED' if err3 < 0.01 else 'FAILED'} (err = {err3:.6f} K)")
    print(f"Case 4 (Realistic L)           : {'PASSED' if err4 < 0.30 else 'FAILED'} (err = {err4:.6f} K)")
    print(f"Case 5 (Mushy Start & Melt)    : {'PASSED' if err5 < 0.01 else 'FAILED'} (err = {err5:.6f} K)")
    print(f"Case 6 (Cooling Reversal)      : {'PASSED' if err6 < 0.01 else 'FAILED'} (err = {err6:.6f} K)")
    print(f"Case 7 (Heating Reversal)      : {'PASSED' if err7 < 0.01 else 'FAILED'} (err = {err7:.6f} K)")
    print(f"Case 8 (Exact Path Integrated) : {'PASSED' if err8 < 0.01 else 'FAILED'} (err = {err8:.6f} K)")
    print(f"Case 9 (Variable Density)      : {'PASSED' if err9 < 0.01 else 'FAILED'} (err = {err9:.6f} K)")
    print(f"Case 10 (Unequal Cp Cooling)   : {'PASSED' if err10 < 0.01 else 'FAILED'} (err = {err10:.6f} K)")
    print(f"Case 11 (Unequal Cp Reversal)  : {'PASSED' if err11 < 0.01 else 'FAILED'} (err = {err11:.6f} K)")
    print(f"Case 12 (Closed Cycle Net H=0) : {'PASSED' if err12 < 0.01 else 'FAILED'} (err = {err12:.6f} K)")
    print(f"Case 13 (Multi-Step Reversal)  : {'PASSED' if err13 < 0.001 else 'FAILED'} (err = {err13:.6f} K)")
    print(f"Case 14 (Unequal Cp Closed)    : {'PASSED' if err14 < 0.01 else 'FAILED'} (err = {err14:.6f} K)")
    print(f"Case 15 (Liquid Reversal)      : {'PASSED' if err15 < 0.01 else 'FAILED'} (err = {err15:.6f} K)")
    print("-------------------------------------------------------")

    if all_passed:
        print("\nALL SINGLE-CELL VERIFICATION TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSOME TESTS FAILED!")
        sys.exit(1)

if __name__ == "__main__":
    main()
