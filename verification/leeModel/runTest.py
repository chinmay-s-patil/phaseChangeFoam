#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import math

def find_solver():
    home = os.environ.get("HOME", "/home/lavender")
    matches = [
        os.path.join(home, "OpenFOAM", "lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam")
    ]
    for m in matches:
        if os.path.exists(m):
            return m
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        path = os.path.join(user_appbin, "phaseChangeMultiRegionFoam")
        if os.path.exists(path):
            return path
    path = shutil.which("phaseChangeMultiRegionFoam")
    if path:
        return path
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

def setup_case(case_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=False, p_init=101325.0, c_evap=0.1, c_cond=0.1, dt=1.0, n_outer_corr=1, q_source=0.0, end_time=None, cp_l=1000.0, cp_v=1000.0, hf_l=0.0, hf_v=0.0):
    os.makedirs(case_dir, exist_ok=True)
    end_t = end_time if end_time is not None else dt
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end_t}; deltaT {dt};
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
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }}
solvers {{
    "h.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "T.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "Yi.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "air.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "H2O_l.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "H2O_v.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "p.*" {{ solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }}
    "p_rgh.*" {{ solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }}
    "U.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }}
    "rho.*" {{ solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }}
}}
PIMPLE
{{
    nOuterCorrectors {n_outer_corr};
    nCorrectors 1;
    nNonLinearCorrectors 1;
    pRefCell 0;
    pRefValue 101325;
}}
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
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/air"; object thermophysicalProperties; }}
thermoType
{{
    type            heRhoThermo;
    mixture         multiComponentMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}}
dpdt            no;

species
(
    air
    H2O_l
    H2O_v
);

inertSpecie air;

air
{{
    specie {{ molWeight 28.96; }}
    thermodynamics {{ Cp 1000; Hf 0; }}
    transport {{ mu 1.8e-5; Pr 0.7; }}
}}

H2O_l
{{
    specie {{ molWeight 18.015; }}
    thermodynamics {{ Cp {cp_l}; Hf {hf_l}; }}
    transport {{ mu 1.8e-5; Pr 0.7; }}
}}

H2O_v
{{
    specie {{ molWeight 18.015; }}
    thermodynamics {{ Cp {cp_v}; Hf {hf_v}; }}
    transport {{ mu 1.8e-5; Pr 0.7; }}
}}
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

    if q_source > 0.0:
        with open(os.path.join(air_const, "fvOptions"), "w") as f:
            f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/air"; object fvOptions; }}
heatSource
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   all;
    scalarSemiImplicitSourceCoeffs
    {{
        selectionMode   all;
        volumeMode      specific;
        injectionRateSuSp
        {{
            h           ({q_source} 0);
        }}
    }}
}}
""")

    # Mesh generation and region directory copies
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    run_cmd(f"cd {case_dir} && bash -c '{of_env}; blockMesh'", cwd=case_dir)
    run_cmd(f"mkdir -p {case_dir}/constant/air {case_dir}/system/air")
    run_cmd(f"cp -r {case_dir}/constant/polyMesh {case_dir}/constant/air/polyMesh 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSchemes {case_dir}/system/air/fvSchemes 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSolution {case_dir}/system/air/fvSolution 2>/dev/null || true")
    if q_source > 0.0:
        run_cmd(f"cp {air_const}/fvOptions {case_dir}/system/air/fvOptions 2>/dev/null || true")

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

    # Test 6: Evaporation Multi-Outer-Corrector Invariance Verification (nOuterCorrectors in [2, 3, 4] vs 1)
    print("\n--- Test 6: Evaporation Multi-Outer-Corrector Invariance Verification (nOuterCorrectors in [2, 3, 4]) ---")
    for n_corr in [2, 3, 4]:
        c6_dir = os.path.join(base_dir, f"case_evap_multi_outer_{n_corr}")
        if os.path.exists(c6_dir):
            shutil.rmtree(c6_dir)
        setup_case(c6_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, c_evap=100.0, dt=1.0, n_outer_corr=n_corr)
        run_cmd(f"cd {c6_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c6_dir)

        T_final6 = parse_openfoam_field(os.path.join(c6_dir, "1/air/T"))[0]
        Yl_final6 = parse_openfoam_field(os.path.join(c6_dir, "1/air/H2O_l"))[0]

        print(f"Evap nOuter={n_corr} Final T:  {T_final6:.6f} K (nOuter=1 Final T: {T_final4:.6f} K)")
        print(f"Evap nOuter={n_corr} Final Yl: {Yl_final6:.6f} (nOuter=1 Final Yl: {Yl_final4:.6f})")

        assert abs(T_final6 - T_final4) < 1e-4, f"FAIL: Evap nOuter={n_corr} final T ({T_final6:.4f} K) differs from nOuter=1 ({T_final4:.4f} K)!"
        assert abs(Yl_final6 - Yl_final4) < 1e-6, f"FAIL: Evap nOuter={n_corr} final Yl ({Yl_final6:.6f}) differs from nOuter=1 ({Yl_final4:.6f})!"
        assert abs(T_final6 - 373.15) < 1e-4, f"FAIL: Evap nOuter={n_corr} final T ({T_final6:.4f} K) did not land on Tsat (373.15 K)!"
    print("PASS: Evaporation Multi-Outer-Corrector Invariance Check Passed!")

    # Test 7: Condensation Multi-Outer-Corrector Invariance Verification (C_cond=100/s, dt=1s, nOuterCorrectors in [2, 3, 4] vs 1)
    print("\n--- Test 7: Condensation Multi-Outer-Corrector Invariance Verification (nOuterCorrectors in [2, 3, 4]) ---")
    c7_n1_dir = os.path.join(base_dir, "case_cond_large_Cdt_n1")
    if os.path.exists(c7_n1_dir):
        shutil.rmtree(c7_n1_dir)
    setup_case(c7_n1_dir, T_init=360.0, Y_l_init=0.0, Y_v_init=0.2, c_cond=100.0, dt=1.0, n_outer_corr=1)
    run_cmd(f"cd {c7_n1_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c7_n1_dir)

    T_final7_n1 = parse_openfoam_field(os.path.join(c7_n1_dir, "1/air/T"))[0]
    Yv_final7_n1 = parse_openfoam_field(os.path.join(c7_n1_dir, "1/air/H2O_v"))[0]

    for n_corr in [2, 3, 4]:
        c7_n_dir = os.path.join(base_dir, f"case_cond_multi_outer_{n_corr}")
        if os.path.exists(c7_n_dir):
            shutil.rmtree(c7_n_dir)
        setup_case(c7_n_dir, T_init=360.0, Y_l_init=0.0, Y_v_init=0.2, c_cond=100.0, dt=1.0, n_outer_corr=n_corr)
        run_cmd(f"cd {c7_n_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c7_n_dir)

        T_final7_n = parse_openfoam_field(os.path.join(c7_n_dir, "1/air/T"))[0]
        Yv_final7_n = parse_openfoam_field(os.path.join(c7_n_dir, "1/air/H2O_v"))[0]

        print(f"Cond nOuter={n_corr} Final T:  {T_final7_n:.6f} K (nOuter=1 Final T: {T_final7_n1:.6f} K)")
        print(f"Cond nOuter={n_corr} Final Yv: {Yv_final7_n:.6f} (nOuter=1 Final Yv: {Yv_final7_n1:.6f})")

        assert abs(T_final7_n - T_final7_n1) < 1e-4, f"FAIL: Cond nOuter={n_corr} final T ({T_final7_n:.4f} K) differs from nOuter=1 ({T_final7_n1:.4f} K)!"
        assert abs(Yv_final7_n - Yv_final7_n1) < 1e-6, f"FAIL: Cond nOuter={n_corr} final Yv ({Yv_final7_n:.6f}) differs from nOuter=1 ({Yv_final7_n1:.6f})!"
        assert abs(T_final7_n - 373.15) < 1e-4, f"FAIL: Cond nOuter={n_corr} final T ({T_final7_n:.4f} K) did not land on Tsat (373.15 K)!"
    print("PASS: Condensation Multi-Outer-Corrector Invariance Check Passed!")

    # Test 8: Heated Pot Physical Energy Balance & Thermal Cap Plateau Verification
    print("\n--- Test 8: Heated Pot Physical Energy Balance Verification ---")
    c8_dir = os.path.join(base_dir, "case_heated_pot")
    if os.path.exists(c8_dir):
        shutil.rmtree(c8_dir)

    Q_in = 100000.0  # 1e5 W/m^3
    dt8 = 0.01
    t_end8 = 3.0
    T_init8 = 350.0
    Tsat8 = 373.15
    Yl_init8 = 0.20
    Yv_init8 = 0.0
    L8 = 2.26e6
    Cp8 = 1000.0

    setup_case(c8_dir, T_init=T_init8, Y_l_init=Yl_init8, Y_v_init=Yv_init8, c_evap=1000.0, dt=dt8, end_time=t_end8, q_source=Q_in)
    run_cmd(f"cd {c8_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c8_dir)

    # 1. Intermediate Plateau Check at t = 0.25s and t = 0.50s (cap solidly binding with C_evap = 1000/s)
    def find_time_dir(target_t, case_path):
        dirs = [d for d in os.listdir(case_path) if d.replace('.', '', 1).isdigit() and abs(float(d) - target_t) < 1e-4]
        return dirs[0] if dirs else None

    t25_dir = find_time_dir(0.25, c8_dir)
    t50_dir = find_time_dir(0.50, c8_dir)
    assert t25_dir, "FAIL: Time directory for t=0.25s not found in case_heated_pot!"
    assert t50_dir, "FAIL: Time directory for t=0.50s not found in case_heated_pot!"
    latest_time_dir = max([d for d in os.listdir(c8_dir) if d.replace('.', '', 1).isdigit() and float(d) > 0], key=lambda x: float(x))
    assert latest_time_dir, "FAIL: Final time directory not found in case_heated_pot!"

    T_25 = parse_openfoam_field(os.path.join(c8_dir, f"{t25_dir}/air/T"))[0]
    Yl_25 = parse_openfoam_field(os.path.join(c8_dir, f"{t25_dir}/air/H2O_l"))[0]
    rho_25 = parse_openfoam_field(os.path.join(c8_dir, f"{t25_dir}/air/rho"))[0]

    T_50 = parse_openfoam_field(os.path.join(c8_dir, f"{t50_dir}/air/T"))[0]
    Yl_50 = parse_openfoam_field(os.path.join(c8_dir, f"{t50_dir}/air/H2O_l"))[0]

    # Exact analytical boiling plateau temperature under explicit heat input Q:
    # T_plateau = Tsat + Q * dt / (rho * Cp)
    T_plateau_anal = Tsat8 + (Q_in * dt8) / (rho_25 * Cp8)

    print(f"Intermediate t = 0.25s: T = {T_25:.6f} K, Yl = {Yl_25:.6f}")
    print(f"Intermediate t = 0.50s: T = {T_50:.6f} K, Yl = {Yl_50:.6f}")
    print(f"Analytic Plateau T:     {T_plateau_anal:.6f} K (Tsat = {Tsat8} K)")

    # Assert flat boiling plateau: T(t=0.25s) == T(t=0.50s) to within 0.01 K
    assert abs(T_25 - T_50) < 0.01, f"FAIL: Boiling plateau temperature drifted between t=0.25s ({T_25:.4f} K) and t=0.50s ({T_50:.4f} K)!"
    # Assert simulation plateau temperature matches exact analytical T_plateau to within 0.02 K
    assert abs(T_25 - T_plateau_anal) < 0.02, f"FAIL: Plateau T ({T_25:.4f} K) differs from analytical T_plateau ({T_plateau_anal:.4f} K) by > 0.02 K!"

    # 2. Final Energy Balance Check at t = 3.0s
    T_final8 = parse_openfoam_field(os.path.join(c8_dir, f"{latest_time_dir}/air/T"))[0]
    Yl_final8 = parse_openfoam_field(os.path.join(c8_dir, f"{latest_time_dir}/air/H2O_l"))[0]
    Yv_final8 = parse_openfoam_field(os.path.join(c8_dir, f"{latest_time_dir}/air/H2O_v"))[0]
    rho8 = parse_openfoam_field(os.path.join(c8_dir, f"{latest_time_dir}/air/rho"))[0]

    E_input = Q_in * float(latest_time_dir)
    E_sensible = rho8 * Cp8 * (T_final8 - T_init8)
    E_latent = rho8 * L8 * (Yv_final8 - Yv_init8)
    E_stored = E_sensible + E_latent
    rel_err8 = abs(E_input - E_stored) / E_input

    print(f"Total Heat Injected (Q*t): {E_input:.2f} J/m^3")
    print(f"Sensible Energy Stored:    {E_sensible:.2f} J/m^3")
    print(f"Latent Energy Stored:      {E_latent:.2f} J/m^3")
    print(f"Total Energy Stored:       {E_stored:.2f} J/m^3 (rel error: {rel_err8:.4%})")
    print(f"Final T (t=3.0s):          {T_final8:.4f} K")
    print(f"Final Yl (t=3.0s):         {Yl_final8:.6f} (Yl_init = {Yl_init8})")

    assert rel_err8 < 0.01, f"FAIL: Energy balance relative error ({rel_err8:.4%}) > 1%!"
    print("PASS: Heated Pot Physical Energy Balance Verification Passed!")

    # Test 8c: dt-Halving Convergence Check (dt=0.01s -> 0.005s) to Physical Tsat
    print("\n--- Test 8c: dt-Halving Plateau Offset Convergence Verification (dt = 0.01s -> 0.005s) ---")
    c8c_dir = os.path.join(base_dir, "case_heated_pot_half_dt")
    if os.path.exists(c8c_dir):
        shutil.rmtree(c8c_dir)

    dt8c = 0.005
    setup_case(c8c_dir, T_init=T_init8, Y_l_init=Yl_init8, Y_v_init=Yv_init8, c_evap=1000.0, dt=dt8c, end_time=0.5, q_source=Q_in)
    run_cmd(f"cd {c8c_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c8c_dir)

    t50_dir_8c = find_time_dir(0.50, c8c_dir)
    assert t50_dir_8c, "FAIL: Time directory for t=0.50s not found in case_heated_pot_half_dt!"
    T_50_8c = parse_openfoam_field(os.path.join(c8c_dir, f"{t50_dir_8c}/air/T"))[0]

    offset_dt1 = T_25 - Tsat8          # Offset with dt = 0.01s at t=0.25s
    offset_dt2 = T_50_8c - Tsat8       # Offset with dt = 0.005s at t=0.50s

    ratio = offset_dt2 / offset_dt1
    print(f"Plateau offset (dt=0.01s):   delta1 = {offset_dt1:.6f} K")
    print(f"Plateau offset (dt=0.005s):  delta2 = {offset_dt2:.6f} K")
    print(f"Offset Ratio (delta2/delta1): {ratio:.4f} (Expected: 0.5000)")

    assert abs(ratio - 0.50) < 0.02, f"FAIL: Offset ratio ({ratio:.4f}) differs from expected first-order dt scaling (0.5000) by > 2%!"
    print("PASS: dt-Halving Plateau Offset Convergence to Physical Tsat Passed!")

    # Test 8b: Heated Pot Liquid Exhaustion & Superheating Verification (Yl_init = 0.05)
    print("\n--- Test 8b: Heated Pot Liquid Exhaustion & Superheating Verification (Yl_init = 0.05) ---")
    c8b_dir = os.path.join(base_dir, "case_heated_pot_exhaustion")
    if os.path.exists(c8b_dir):
        shutil.rmtree(c8b_dir)

    Yl_init8b = 0.05
    setup_case(c8b_dir, T_init=T_init8, Y_l_init=Yl_init8b, Y_v_init=Yv_init8, c_evap=100.0, dt=dt8, end_time=t_end8, q_source=Q_in)
    run_cmd(f"cd {c8b_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c8b_dir)

    latest_time_dir_8b = max([d for d in os.listdir(c8b_dir) if d.replace('.', '', 1).isdigit() and float(d) > 0], key=lambda x: float(x))
    assert latest_time_dir_8b, "FAIL: Final time directory not found in case_heated_pot_exhaustion!"

    t05_dir_8b = find_time_dir(0.50, c8b_dir)
    t10_dir_8b = find_time_dir(1.00, c8b_dir)
    assert t05_dir_8b, "FAIL: Time directory for t=0.50s not found in case_heated_pot_exhaustion!"
    assert t10_dir_8b, "FAIL: Time directory for t=1.00s not found in case_heated_pot_exhaustion!"

    # 1. Intermediate Rate Model Verification (during active boiling before exhaustion)
    T_05_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t05_dir_8b}/air/T"))[0]
    Yl_05_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t05_dir_8b}/air/H2O_l"))[0]
    rho_05_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t05_dir_8b}/air/rho"))[0]

    T_10_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t10_dir_8b}/air/T"))[0]
    Yl_10_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t10_dir_8b}/air/H2O_l"))[0]
    rho_10_8b = parse_openfoam_field(os.path.join(c8b_dir, f"{t10_dir_8b}/air/rho"))[0]

    # Analytical intermediate liquid mass fractions matching energy split at t=0.5s and t=1.0s:
    E_sens_05 = rho_05_8b * Cp8 * (T_05_8b - T_init8)
    E_sens_10 = rho_10_8b * Cp8 * (T_10_8b - T_init8)

    deltaYl_05_anal = - (Q_in * 0.50 - E_sens_05) / (rho_05_8b * L8)
    deltaYl_10_anal = - (Q_in * 1.00 - E_sens_10) / (rho_10_8b * L8)
    Yl_05_anal = Yl_init8b + deltaYl_05_anal
    Yl_10_anal = Yl_init8b + deltaYl_10_anal

    err_Yl_05 = abs(Yl_05_8b - Yl_05_anal) / Yl_05_anal
    err_Yl_10 = abs(Yl_10_8b - Yl_10_anal) / Yl_10_anal

    print(f"Intermediate t = 0.50s: Yl = {Yl_05_8b:.6f} (Analytic: {Yl_05_anal:.6f}, rel err: {err_Yl_05:.6e})")
    print(f"Intermediate t = 1.00s: Yl = {Yl_10_8b:.6f} (Analytic: {Yl_10_anal:.6f}, rel err: {err_Yl_10:.6e})")

    assert err_Yl_05 < 1e-4, f"FAIL: Intermediate Yl at t=0.5s ({Yl_05_8b:.6f}) differs from analytical ({Yl_05_anal:.6f})!"
    assert err_Yl_10 < 1e-4, f"FAIL: Intermediate Yl at t=1.0s ({Yl_10_8b:.6f}) differs from analytical ({Yl_10_anal:.6f})!"

    # 2. Analytical final temperature calculation after liquid depletion:
    # All liquid Yl_init (0.05) is evaporated (consuming rho * L * Yl_init energy).
    # The remaining energy Q*t - rho*L*Yl_init goes into sensible heating: rho * Cp * (T_final - T_init).
    # T_final_anal = T_init + (Q*t - rho*L*Yl_init) / (rho * Cp)
    T_final8b = parse_openfoam_field(os.path.join(c8b_dir, f"{latest_time_dir_8b}/air/T"))[0]
    Yl_final8b = parse_openfoam_field(os.path.join(c8b_dir, f"{latest_time_dir_8b}/air/H2O_l"))[0]
    Yv_final8b = parse_openfoam_field(os.path.join(c8b_dir, f"{latest_time_dir_8b}/air/H2O_v"))[0]
    rho8b = parse_openfoam_field(os.path.join(c8b_dir, f"{latest_time_dir_8b}/air/rho"))[0]

    E_input_8b = Q_in * float(latest_time_dir_8b)
    E_latent_8b = rho8b * L8 * Yl_init8b
    E_sensible_needed_8b = E_input_8b - E_latent_8b
    T_final_anal_8b = T_init8 + E_sensible_needed_8b / (rho8b * Cp8)

    E_sensible_8b = rho8b * Cp8 * (T_final8b - T_init8)
    E_stored_8b = E_sensible_8b + rho8b * L8 * (Yv_final8b - Yv_init8)
    rel_err_8b = abs(E_input_8b - E_stored_8b) / E_input_8b
    rel_T_err_8b = abs(T_final8b - T_final_anal_8b) / T_final_anal_8b

    print(f"Liquid Exhaustion Final Yl: {Yl_final8b:.6f} (Expected: 0.000000)")
    print(f"Liquid Exhaustion Final T:  {T_final8b:.4f} K (Analytic T_final: {T_final_anal_8b:.4f} K, rel error: {rel_T_err_8b:.6%})")
    print(f"Sensible Energy fraction:   {E_sensible_8b / E_input_8b:.2%}")
    print(f"Energy Balance rel error:   {rel_err_8b:.6%}")

    # Assertions:
    assert Yl_final8b >= 0.0, f"FAIL: Liquid mass fraction ({Yl_final8b}) fell below 0!"
    assert Yl_final8b < 1e-5, f"FAIL: Liquid mass fraction ({Yl_final8b}) was not completely exhausted!"
    assert rel_T_err_8b < 1e-4, f"FAIL: Final temperature ({T_final8b:.4f} K) differs from analytical superheated T ({T_final_anal_8b:.4f} K) by > 1e-4!"
    assert rel_err_8b < 1e-4, f"FAIL: Energy balance relative error ({rel_err_8b:.6%}) > 1e-4!"
    print("PASS: Heated Pot Liquid Exhaustion & Superheating Verification Passed!")

    # Test 9: Closed Adiabatic Unequal-Cp Absolute Enthalpy Conservation Verification
    print("\n--- Test 9: Closed Adiabatic Unequal-Cp Absolute Enthalpy Conservation Verification ---")
    c9_dir = os.path.join(base_dir, "case_unequal_cp_enthalpy_conservation")
    if os.path.exists(c9_dir):
        shutil.rmtree(c9_dir)

    C_evap9 = 0.01
    dt9 = 0.1
    T_init9 = 380.0
    Yl_init9 = 0.20
    Yv_init9 = 0.0
    Yair_init9 = 0.80
    L9 = 2.26e6

    Cp_air = 1000.0
    Cp_l = 4184.0   # Liquid Cp (water)
    Cp_v = 2030.0   # Vapor Cp (steam) - distinctly unequal from Cp_l!
    Hf_l = 0.0
    Hf_v = L9       # Vapor formation enthalpy encodes latent heat L

    T_std = 298.15

    setup_case(
        c9_dir,
        T_init=T_init9,
        Y_l_init=Yl_init9,
        Y_v_init=Yv_init9,
        c_evap=C_evap9,
        dt=dt9,
        cp_l=Cp_l,
        cp_v=Cp_v,
        hf_l=Hf_l,
        hf_v=Hf_v
    )

    run_cmd(f"cd {c9_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c9_dir)

    T_final9 = parse_openfoam_field(os.path.join(c9_dir, f"{dt9}/air/T"))[0]
    Yl_final9 = parse_openfoam_field(os.path.join(c9_dir, f"{dt9}/air/H2O_l"))[0]
    Yv_final9 = parse_openfoam_field(os.path.join(c9_dir, f"{dt9}/air/H2O_v"))[0]
    Yair_final9 = parse_openfoam_field(os.path.join(c9_dir, f"{dt9}/air/air"))[0]

    W_air = 28.96
    W_H2O = 18.015
    R_univ = 8314.463
    W_mix0 = 1.0 / (Yair_init9 / W_air + Yl_init9 / W_H2O + Yv_init9 / W_H2O)
    R_mix0 = R_univ / W_mix0
    rho9_0 = 101325.0 / (R_mix0 * T_init9)

    rho9_1 = parse_openfoam_field(os.path.join(c9_dir, f"{dt9}/air/rho"))[0]

    # Calculate absolute enthalpy per unit mass at initial state t=0:
    # h_abs = sum_i Y_i * [ Cp_i * (T - T_std) + Hf_i ]
    h_abs_init = (
        Yair_init9 * (Cp_air * (T_init9 - T_std)) +
        Yl_init9   * (Cp_l   * (T_init9 - T_std) + Hf_l) +
        Yv_init9   * (Cp_v   * (T_init9 - T_std) + Hf_v)
    )
    E_abs_init = rho9_0 * h_abs_init

    # Calculate absolute enthalpy per unit mass at final state t=dt:
    h_abs_final = (
        Yair_final9 * (Cp_air * (T_final9 - T_std)) +
        Yl_final9   * (Cp_l   * (T_final9 - T_std) + Hf_l) +
        Yv_final9   * (Cp_v   * (T_final9 - T_std) + Hf_v)
    )
    E_abs_final = rho9_1 * h_abs_final

    err_E9 = abs(E_abs_final - E_abs_init) / E_abs_init

    print(f"Initial T:              {T_init9:.4f} K")
    print(f"Final T:                {T_final9:.4f} K")
    print(f"Evaporated Yv:          {Yv_final9:.8e}")
    print(f"Initial E_abs:          {E_abs_init:.6f} J/m^3")
    print(f"Final E_abs:            {E_abs_final:.6f} J/m^3")
    print(f"Enthalpy Rel Error:     {err_E9:.6e}")

    assert err_E9 < 1e-4, f"FAIL: Closed adiabatic unequal-Cp absolute enthalpy conservation error ({err_E9:.6e}) > 1e-4!"
    print("PASS: Closed Adiabatic Unequal-Cp Absolute Enthalpy Conservation Passed!")

    print("\n=========================================================================")
    print("      LEE FLUID PHASE CHANGE MODEL QUANTITATIVE VERIFICATION PASS        ")
    print("=========================================================================")

if __name__ == "__main__":
    main()
