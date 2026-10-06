#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import math

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

def run_cmd(cmd, cwd=None, allow_failure=False):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0 and not allow_failure:
        print(f"Error executing command: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res

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

def setup_case(case_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=False, p_init=101325.0, c_evap=0.1, c_cond=0.1, dt=1.0):
    os.makedirs(case_dir, exist_ok=True)
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {dt}; deltaT {dt};
writeControl timeStep; writeInterval 1; purgeWrite 0; writeFormat ascii; writePrecision 12;
""")

    # system/fvSchemes
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss linear; div(phi,Yi_h) Gauss upwind; div(phi,he) Gauss upwind; div(phiv,p) Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")

    # system/fvSolution
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers {
    "h.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "T.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "Yi.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "air.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "H2O_l.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "H2O_v.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "p.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "p_rgh.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 1;
    nNonLinearCorrectors 1;
    pRefCell 0;
    pRefValue 101325;
}
""")

    # system/blockMeshDict
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }
scale 1;
vertices ( (0 0 0) (0.01 0 0) (0.01 0.01 0) (0 0.01 0) (0 0 0.01) (0.01 0 0.01) (0.01 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) (1 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    walls { type wall; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }
);
""")

    # constant/g and regionProperties
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    with open(os.path.join(case_dir, "constant/g"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0]; value (0 0 0);
""")
    with open(os.path.join(case_dir, "constant/regionProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid (air) solid () porousFluid () porousSolid () );
""")

    # constant/air/
    air_const = os.path.join(case_dir, "constant/air")
    os.makedirs(air_const, exist_ok=True)

    with open(os.path.join(air_const, "thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/air"; object thermophysicalProperties; }
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

species
(
    air
    H2O_l
    H2O_v
);

inertSpecie air;

air
{
    specie { molWeight 28.96; }
    thermodynamics { Cp 1000; Hf 0; }
    transport { mu 1.8e-5; Pr 0.7; }
}

H2O_l
{
    specie { molWeight 18.015; }
    thermodynamics { Cp 1000; Hf 0; }
    transport { mu 1.8e-5; Pr 0.7; }
}

H2O_v
{
    specie { molWeight 18.015; }
    thermodynamics { Cp 1000; Hf 0; }
    transport { mu 1.8e-5; Pr 0.7; }
}
""")

    enable_tsat_p_str = "true" if enable_tsat_p else "false"
    with open(os.path.join(air_const, "phaseChangeDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/air"; object phaseChangeDict; }}
phaseChange
{{
    active true;
    type Lee;
    liquid H2O_l;
    vapor H2O_v;
    C_evap {c_evap};
    C_cond {c_cond};
    latentHeat 2.26e6;
    Tsat 373.15;
    pRef 101325;
    enableTsatP {enable_tsat_p_str};
}}
""")

    with open(os.path.join(air_const, "turbulenceProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/air"; object turbulenceProperties; }
simulationType laminar;
""")

    with open(os.path.join(air_const, "radiationProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/air"; object radiationProperties; }
radiationModel none;
""")

    # 0/air fields
    air_zero = os.path.join(case_dir, "0/air")
    os.makedirs(air_zero, exist_ok=True)

    h_init = T_init * 1000.0

    with open(os.path.join(air_zero, "T"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object T; }}
dimensions [0 0 0 1 0 0 0]; internalField uniform {T_init};
boundaryField {{ walls {{ type zeroGradient; }} }}
""")

    with open(os.path.join(air_zero, "h"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object h; }}
dimensions [0 2 -2 0 0 0 0]; internalField uniform {h_init};
boundaryField {{ walls {{ type zeroGradient; }} }}
""")

    with open(os.path.join(air_zero, "p"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object p; }}
dimensions [1 -1 -2 0 0 0 0]; internalField uniform {p_init};
boundaryField {{ walls {{ type calculated; value uniform {p_init}; }} }}
""")

    with open(os.path.join(air_zero, "p_rgh"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object p_rgh; }}
dimensions [1 -1 -2 0 0 0 0]; internalField uniform {p_init};
boundaryField {{ walls {{ type fixedFluxPressure; value uniform {p_init}; }} }}
""")

    with open(os.path.join(air_zero, "U"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volVectorField; location "0/air"; object U; }
dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
boundaryField { walls { type fixedValue; value uniform (0 0 0); } }
""")

    Y_air_init = max(0.0, 1.0 - Y_l_init - Y_v_init)
    with open(os.path.join(air_zero, "Ydefault"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object Ydefault; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.0;
boundaryField { walls { type zeroGradient; } }
""")

    with open(os.path.join(air_zero, "air"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object air; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {Y_air_init};
boundaryField {{ walls {{ type zeroGradient; }} }}
""")

    with open(os.path.join(air_zero, "H2O_l"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object H2O_l; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {Y_l_init};
boundaryField {{ walls {{ type zeroGradient; }} }}
""")

    with open(os.path.join(air_zero, "H2O_v"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/air"; object H2O_v; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform {Y_v_init};
boundaryField {{ walls {{ type zeroGradient; }} }}
""")

    # Mesh generation and region directory copies
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    run_cmd(f"cd {case_dir} && bash -c '{of_env}; blockMesh'", cwd=case_dir)
    run_cmd(f"mkdir -p {case_dir}/constant/air {case_dir}/system/air")
    run_cmd(f"cp -r {case_dir}/constant/polyMesh {case_dir}/constant/air/polyMesh 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSchemes {case_dir}/system/air/fvSchemes 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSolution {case_dir}/system/air/fvSolution 2>/dev/null || true")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    print("=== Lee Fluid Phase Change Model Quantitative Verification Suite ===")

    # Test 1: Quantitative Evaporation Check (Small dt=0.1s, C=0.01/s)
    print("\n--- Test 1: Quantitative Small dt Evaporation Verification ---")
    c1_dir = os.path.join(base_dir, "case_evaporation_quant")
    if os.path.exists(c1_dir):
        shutil.rmtree(c1_dir)
    
    C_evap1 = 0.01
    dt1 = 0.1
    T_init1 = 380.0
    Tsat1 = 373.15
    Yl_init1 = 0.20
    L1 = 2.26e6
    Cp1 = 1000.0

    setup_case(c1_dir, T_init=T_init1, Y_l_init=Yl_init1, Y_v_init=0.0, c_evap=C_evap1, dt=dt1)
    run_cmd(f"cd {c1_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c1_dir)

    T_final1 = parse_openfoam_field(os.path.join(c1_dir, f"{dt1}/air/T"))[0]
    Yl_final1 = parse_openfoam_field(os.path.join(c1_dir, f"{dt1}/air/H2O_l"))[0]
    Yv_final1 = parse_openfoam_field(os.path.join(c1_dir, f"{dt1}/air/H2O_v"))[0]

    # Analytical formulas (mixture density rho cancels out!):
    # dYl = - C_evap * Yl * (T - Tsat) / Tsat * dt
    deltaYl_anal1 = - C_evap1 * Yl_init1 * (T_init1 - Tsat1) / Tsat1 * dt1
    # dT = - mDot * L * dt / (rho * Cp) = - C_evap * Yl * (T - Tsat) / Tsat * L * dt / Cp
    deltaT_anal1 = - C_evap1 * Yl_init1 * (T_init1 - Tsat1) / Tsat1 * L1 * dt1 / Cp1

    deltaYl_num1 = Yl_final1 - Yl_init1
    deltaT_num1 = T_final1 - T_init1

    err_Y1 = abs(deltaYl_num1 - deltaYl_anal1) / abs(deltaYl_anal1)
    err_T1 = abs(deltaT_num1 - deltaT_anal1) / abs(deltaT_anal1)

    print(f"Numerical deltaYl:  {deltaYl_num1:.8e}")
    print(f"Analytic deltaYl:   {deltaYl_anal1:.8e} (rel error: {err_Y1:.6e})")
    print(f"Numerical deltaT:   {deltaT_num1:.6f} K")
    print(f"Analytic deltaT:    {deltaT_anal1:.6f} K (rel error: {err_T1:.6e})")

    assert err_Y1 < 1e-3, f"FAIL: deltaYl relative error ({err_Y1}) > 1e-3!"
    assert err_T1 < 1e-3, f"FAIL: deltaT relative error ({err_T1}) > 1e-3!"
    print("PASS: Quantitative Evaporation Verification Passed!")

    # Test 2: Quantitative Condensation Check (Small dt=0.1s, C=0.01/s)
    print("\n--- Test 2: Quantitative Small dt Condensation Verification ---")
    c2_dir = os.path.join(base_dir, "case_condensation_quant")
    if os.path.exists(c2_dir):
        shutil.rmtree(c2_dir)
    
    C_cond2 = 0.01
    dt2 = 0.1
    T_init2 = 360.0
    Yv_init2 = 0.20

    setup_case(c2_dir, T_init=T_init2, Y_l_init=0.0, Y_v_init=Yv_init2, c_cond=C_cond2, dt=dt2)
    run_cmd(f"cd {c2_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c2_dir)

    T_final2 = parse_openfoam_field(os.path.join(c2_dir, f"{dt2}/air/T"))[0]
    Yv_final2 = parse_openfoam_field(os.path.join(c2_dir, f"{dt2}/air/H2O_v"))[0]

    # Analytical formulas:
    # dYv = - C_cond * Yv * (Tsat - T) / Tsat * dt
    deltaYv_anal2 = - C_cond2 * Yv_init2 * (Tsat1 - T_init2) / Tsat1 * dt2
    # dT = + C_cond * Yv * (Tsat - T) / Tsat * L * dt / Cp
    deltaT_anal2 = + C_cond2 * Yv_init2 * (Tsat1 - T_init2) / Tsat1 * L1 * dt2 / Cp1

    deltaYv_num2 = Yv_final2 - Yv_init2
    deltaT_num2 = T_final2 - T_init2

    err_Y2 = abs(deltaYv_num2 - deltaYv_anal2) / abs(deltaYv_anal2)
    err_T2 = abs(deltaT_num2 - deltaT_anal2) / abs(deltaT_anal2)

    print(f"Numerical deltaYv:  {deltaYv_num2:.8e}")
    print(f"Analytic deltaYv:   {deltaYv_anal2:.8e} (rel error: {err_Y2:.6e})")
    print(f"Numerical deltaT:   {deltaT_num2:.6f} K")
    print(f"Analytic deltaT:    {deltaT_anal2:.6f} K (rel error: {err_T2:.6e})")

    assert err_Y2 < 1e-3, f"FAIL: deltaYv relative error ({err_Y2}) > 1e-3!"
    assert err_T2 < 1e-3, f"FAIL: deltaT relative error ({err_T2}) > 1e-3!"
    print("PASS: Quantitative Condensation Verification Passed!")

    # Test 3: Quantitative Tsat(p) Clausius-Clapeyron Verification (2 bar)
    print("\n--- Test 3: Quantitative Tsat(p=200 kPa) Clausius-Clapeyron Verification ---")
    p_test = 200000.0 # Pa
    p_ref = 101325.0  # Pa
    R_vap = 8314.463 / 18.015 # J/(kg K) = 461.53
    L_vap = 2.26e6

    # Clausius-Clapeyron formula: 1/Tsat = 1/Tsat0 - (R/L)*ln(p/pRef)
    Tsat_anal_2bar = 1.0 / ((1.0 / Tsat1) - (R_vap / L_vap) * math.log(p_test / p_ref))
    print(f"Analytical Tsat(200 kPa): {Tsat_anal_2bar:.4f} K")

    # Case 3A: T0 = 390 K < Tsat(200 kPa) => Must NOT evaporate (Tsat > 390 K)
    c3a_dir = os.path.join(base_dir, "case_tsat_sub")
    if os.path.exists(c3a_dir):
        shutil.rmtree(c3a_dir)
    setup_case(c3a_dir, T_init=390.0, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=True, p_init=p_test, dt=0.1)
    run_cmd(f"cd {c3a_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c3a_dir)

    Yl_final3a = parse_openfoam_field(os.path.join(c3a_dir, "0.1/air/H2O_l"))[0]
    assert abs(Yl_final3a - 0.20) < 1e-6, "FAIL: Evaporation occurred when T0=390K < Tsat(200kPa)!"
    print("PASS: Sub-saturated state (T0=390K < Tsat=393.9K) correctly prevented evaporation.")

    # Case 3B: T0 = 400 K > Tsat(200 kPa) => Evaporates matching Tsat(200 kPa) formula
    c3b_dir = os.path.join(base_dir, "case_tsat_super")
    if os.path.exists(c3b_dir):
        shutil.rmtree(c3b_dir)
    C_evap3 = 0.01
    dt3 = 0.1
    T0_3b = 400.0
    setup_case(c3b_dir, T_init=T0_3b, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=True, p_init=p_test, c_evap=C_evap3, dt=dt3)
    run_cmd(f"cd {c3b_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c3b_dir)

    Yl_final3b = parse_openfoam_field(os.path.join(c3b_dir, "0.1/air/H2O_l"))[0]
    deltaYl_num3b = Yl_final3b - 0.20
    deltaYl_anal3b = - C_evap3 * 0.20 * (T0_3b - Tsat_anal_2bar) / Tsat_anal_2bar * dt3
    err_tsat = abs(deltaYl_num3b - deltaYl_anal3b) / abs(deltaYl_anal3b)

    print(f"Super-saturated (T0=400K) Numerical deltaYl: {deltaYl_num3b:.8e}")
    print(f"Super-saturated (T0=400K) Analytic deltaYl:  {deltaYl_anal3b:.8e} (rel error: {err_tsat:.6e})")

    assert err_tsat < 1e-3, f"FAIL: Tsat(p) evaporation rate relative error ({err_tsat}) > 1e-3!"
    print("PASS: Quantitative Tsat(p) Clausius-Clapeyron Verification Passed!")

    # Test 4: Large C*dt Limiter Boundary Assertions (C=100/s, dt=1s)
    print("\n--- Test 4: Large C*dt Limiter Boundary Assertions (C=100/s, dt=1s) ---")
    c4_dir = os.path.join(base_dir, "case_large_Cdt")
    if os.path.exists(c4_dir):
        shutil.rmtree(c4_dir)
    setup_case(c4_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, c_evap=100.0, dt=1.0)
    run_cmd(f"cd {c4_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c4_dir)

    T_final4 = parse_openfoam_field(os.path.join(c4_dir, "1/air/T"))[0]
    Yl_final4 = parse_openfoam_field(os.path.join(c4_dir, "1/air/H2O_l"))[0]

    print(f"Large C*dt Final T:    {T_final4:.4f} K (Tsat = 373.1500 K)")
    print(f"Large C*dt Final Yl:   {Yl_final4:.6f} (Yl_init = 0.200000)")

    # Assertion 1: Yl must stay non-negative (mass cap)
    assert Yl_final4 >= 0.0, f"FAIL: Mass limiter breached! Yl ({Yl_final4}) went negative!"
    # Assertion 2: T must NOT drop below Tsat (thermal cap)
    assert T_final4 >= 373.15 - 1e-4, f"FAIL: Thermal limiter breached! T ({T_final4:.4f} K) dropped below Tsat (373.15 K)!"

    # Test 5: Missing Species FatalIOError Verification
    print("\n--- Test 5: Missing Species FatalIOError Verification ---")
    c5_dir = os.path.join(base_dir, "case_missing_species")
    if os.path.exists(c5_dir):
        shutil.rmtree(c5_dir)
    setup_case(c5_dir, T_init=380.0)
    # Overwrite phaseChangeDict with non-existent liquid species name
    with open(os.path.join(c5_dir, "constant/air/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/air"; object phaseChangeDict; }
phaseChange
{
    active true;
    type Lee;
    liquid nonExistentLiquid;
    vapor H2O_v;
    C_evap 0.1;
    C_cond 0.1;
    latentHeat 2.26e6;
    Tsat 373.15;
}
""")
    res5 = run_cmd(f"cd {c5_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c5_dir, allow_failure=True)
    out5 = res5.stdout + res5.stderr
    pass_t5 = (res5.returncode != 0 and "nonExistentLiquid" in out5 and "not found in thermo composition" in out5)
    print(f"Test 5 FatalIOError Triggered: {pass_t5}")
    assert pass_t5, "FAIL: Missing species in active Lee model did not trigger FatalIOError!"
    print("PASS: Missing Species FatalIOError Check Passed!")

    print("\n=========================================================================")
    print("      LEE FLUID PHASE CHANGE MODEL QUANTITATIVE VERIFICATION PASS        ")
    print("=========================================================================")

if __name__ == "__main__":
    main()
