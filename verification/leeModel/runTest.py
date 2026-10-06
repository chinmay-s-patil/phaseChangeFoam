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

def setup_case(case_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=False, p_init=101325.0):
    os.makedirs(case_dir, exist_ok=True)
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
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
    C_evap 0.1;
    C_cond 0.1;
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

    print("=== Lee Fluid Phase Change Model Verification Suite ===")

    # Test 1: In-region Evaporation Check (T0 = 380 K > Tsat = 373.15 K)
    print("\n--- Test 1: In-Region Evaporation Verification (T0 > Tsat) ---")
    c1_dir = os.path.join(base_dir, "case_evaporation")
    if os.path.exists(c1_dir):
        shutil.rmtree(c1_dir)
    setup_case(c1_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0)

    run_cmd(f"cd {c1_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c1_dir)

    T_init1 = parse_openfoam_field(os.path.join(c1_dir, "0/air/T"))[0]
    T_final1 = parse_openfoam_field(os.path.join(c1_dir, "1/air/T"))[0]
    Yl_init1 = parse_openfoam_field(os.path.join(c1_dir, "0/air/H2O_l"))[0]
    Yl_final1 = parse_openfoam_field(os.path.join(c1_dir, "1/air/H2O_l"))[0]
    Yv_init1 = parse_openfoam_field(os.path.join(c1_dir, "0/air/H2O_v"))[0]
    Yv_final1 = parse_openfoam_field(os.path.join(c1_dir, "1/air/H2O_v"))[0]

    print(f"Evaporation T:    {T_init1:.2f} K -> {T_final1:.2f} K (deltaT = {T_final1 - T_init1:.4f} K)")
    print(f"Evaporation H2O_l: {Yl_init1:.4f} -> {Yl_final1:.4f} (deltaYl = {Yl_final1 - Yl_init1:.6e})")
    print(f"Evaporation H2O_v: {Yv_init1:.4f} -> {Yv_final1:.4f} (deltaYv = {Yv_final1 - Yv_init1:.6e})")

    assert T_final1 < T_init1, "FAIL: Temperature did not drop during evaporation (latent heat absorption missing)!"
    assert Yl_final1 < Yl_init1, "FAIL: Liquid species fraction did not deplete during evaporation!"
    assert Yv_final1 > Yv_init1, "FAIL: Vapor species fraction did not increase during evaporation!"
    print("PASS: In-Region Evaporation Verification Passed!")

    # Test 2: In-region Condensation Check (T0 = 360 K < Tsat = 373.15 K)
    print("\n--- Test 2: In-Region Condensation Verification (T0 < Tsat) ---")
    c2_dir = os.path.join(base_dir, "case_condensation")
    if os.path.exists(c2_dir):
        shutil.rmtree(c2_dir)
    setup_case(c2_dir, T_init=360.0, Y_l_init=0.0, Y_v_init=0.2)

    run_cmd(f"cd {c2_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c2_dir)

    T_init2 = parse_openfoam_field(os.path.join(c2_dir, "0/air/T"))[0]
    T_final2 = parse_openfoam_field(os.path.join(c2_dir, "1/air/T"))[0]
    Yl_init2 = parse_openfoam_field(os.path.join(c2_dir, "0/air/H2O_l"))[0]
    Yl_final2 = parse_openfoam_field(os.path.join(c2_dir, "1/air/H2O_l"))[0]
    Yv_init2 = parse_openfoam_field(os.path.join(c2_dir, "0/air/H2O_v"))[0]
    Yv_final2 = parse_openfoam_field(os.path.join(c2_dir, "1/air/H2O_v"))[0]

    print(f"Condensation T:    {T_init2:.2f} K -> {T_final2:.2f} K (deltaT = {T_final2 - T_init2:.4f} K)")
    print(f"Condensation H2O_l: {Yl_init2:.4f} -> {Yl_final2:.4f} (deltaYl = {Yl_final2 - Yl_init2:.6e})")
    print(f"Condensation H2O_v: {Yv_init2:.4f} -> {Yv_final2:.4f} (deltaYv = {Yv_final2 - Yv_init2:.6e})")

    assert T_final2 > T_init2, "FAIL: Temperature did not rise during condensation (latent heat release missing)!"
    assert Yl_final2 > Yl_init2, "FAIL: Liquid species fraction did not increase during condensation!"
    assert Yv_final2 < Yv_init2, "FAIL: Vapor species fraction did not deplete during condensation!"
    print("PASS: In-Region Condensation Verification Passed!")

    # Test 3: Pressure-Dependent Tsat(p) Check
    print("\n--- Test 3: Pressure-Dependent Tsat(p) Check ---")
    c3_dir = os.path.join(base_dir, "case_tsat_p")
    if os.path.exists(c3_dir):
        shutil.rmtree(c3_dir)
    # At p = 200 kPa, Tsat(p) > 373.15 K (~393 K). T0 = 380 K is below Tsat(p=200kPa), so it should NOT evaporate!
    setup_case(c3_dir, T_init=380.0, Y_l_init=0.2, Y_v_init=0.0, enable_tsat_p=True, p_init=200000.0)

    run_cmd(f"cd {c3_dir} && bash -c '{of_env}; {solver_bin}'", cwd=c3_dir)

    Yl_init3 = parse_openfoam_field(os.path.join(c3_dir, "0/air/H2O_l"))[0]
    Yl_final3 = parse_openfoam_field(os.path.join(c3_dir, "1/air/H2O_l"))[0]
    print(f"Elevated Pressure (2 bar) H2O_l: {Yl_init3:.4f} -> {Yl_final3:.4f}")

    assert abs(Yl_final3 - Yl_init3) < 1e-6, "FAIL: Evaporation occurred at T=380K when p=200kPa elevated Tsat above 380K!"
    print("PASS: Tsat(p) Pressure Dependence Verification Passed!")

    print("\n=======================================================")
    print("      LEE FLUID PHASE CHANGE MODEL VERIFICATION PASS   ")
    print("=======================================================")

if __name__ == "__main__":
    main()
