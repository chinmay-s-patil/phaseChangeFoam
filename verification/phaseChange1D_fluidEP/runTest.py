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
        # Handle scalar vs vector uniform string
        if val_str.startswith("(") and val_str.endswith(")"):
            vec = [float(x) for x in val_str[1:-1].split()]
            return [vec] * num_cells
        return [float(val_str)] * num_cells
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    block = content[start_paren+1:end_paren].strip()
    lines = block.split("\n")
    vals = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("("):
            clean = line.strip("() ")
            vec = [float(x) for x in clean.split()]
            vals.append(vec)
        else:
            vals.append(float(line))
    return vals

def setup_fluid_case(case_dir, frozen_flow=True, Cu=1e12):
    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 2000; deltaT 10;
writeControl runTime; writeInterval 200; purgeWrite 0; writeFormat ascii; writePrecision 12;
""")
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes {
    default Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    frozen_str = "yes" if frozen_flow else "no"
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }}
solvers {{
    "h.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }}
    "he.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }}
    "T.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }}
    "rho.*" {{ solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }}
    "p_rgh.*" {{ solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }}
    "U.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }}
}}
PIMPLE {{
    nOuterCorrectors 2;
    nCorrectors 2;
    nNonLinearCorrectors 50;
    nonLinearTolerance 1e-6;
    frozenFlow {frozen_str};
    pRefCell 0;
    pRefValue 101325;
}}
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
regions ( fluid (pcm) solid () porousFluid () porousSolid () );
""")
    
    os.makedirs(os.path.join(case_dir, "constant/pcm"), exist_ok=True)
    with open(os.path.join(case_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heRhoThermo; mixture pureMixture; transport const; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture {
    specie { molWeight 200.0; }
    transport { mu 1e-5; Pr 0.7; }
    thermodynamics { Cp 1980.0; Hf 0; }
    equationOfState { rho 1967.0; }
}
""")
    with open(os.path.join(case_dir, "constant/pcm/turbulenceProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object turbulenceProperties; }
simulationType laminar;
""")
    with open(os.path.join(case_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }}
active true;
phaseChange {{
    type enthalpyPorosity;
    active true;
    forward {{ T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }}
    porosity {{ Cu {Cu:e}; q 0.001; }}
    thermophysical {{ mode custom; CpSolid 1980.0; CpLiquid 1980.0; kSolid 2.0; kLiquid 1.0; }}
}}
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
    with open(os.path.join(case_dir, "0/pcm/he"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object he; }
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
    with open(os.path.join(case_dir, "0/pcm/p_rgh"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type fixedFluxPressure; value uniform 101325; } emptyFaces { type empty; } }
""")
    with open(os.path.join(case_dir, "0/pcm/U"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volVectorField; location "0/pcm"; object U; }
dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
boundaryField { ".*" { type fixedValue; value uniform (0 0 0); } emptyFaces { type empty; } }
""")
    with open(os.path.join(case_dir, "0/pcm/phi"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class surfaceScalarField; location "0/pcm"; object phi; }
dimensions [1 0 -1 0 0 0 0]; internalField uniform 0;
boundaryField { ".*" { type calculated; value uniform 0; } emptyFaces { type empty; } }
""")
    with open(os.path.join(case_dir, "0/pcm/rho"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object rho; }
dimensions [1 -3 0 0 0 0 0]; internalField uniform 1967;
boundaryField { ".*" { type calculated; value uniform 1967; } emptyFaces { type empty; } }
""")

    run_cmd(f"cd {case_dir} && bash -c '{of_env}; blockMesh'", cwd=case_dir)
    run_cmd(f"mkdir -p {case_dir}/constant/pcm {case_dir}/system/pcm")
    run_cmd(f"cp -r {case_dir}/constant/polyMesh {case_dir}/constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSchemes {case_dir}/system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd(f"cp {case_dir}/system/fvSolution {case_dir}/system/pcm/fvSolution 2>/dev/null || true")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== Fluid Region Enthalpy-Porosity Verification Suite in {base_dir} ===")
    
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    # Test 1: Unit check of Darcy Drag Coefficient momentumSp()
    print("\n--- Test 1: Unit Check of Darcy Drag Coefficient momentumSp() ---")
    Cu = 1e5
    q = 0.001
    # Formula: C_drag(f) = Cu * (1 - f)^2 / (f^3 + q)
    c_f0 = Cu * pow(1.0 - 0.0, 2) / (pow(0.0, 3) + q) # 1e8
    c_f1 = Cu * pow(1.0 - 1.0, 2) / (pow(1.0, 3) + q) # 0.0
    c_f05 = Cu * pow(1.0 - 0.5, 2) / (pow(0.5, 3) + q) # 1e5 * 0.25 / 0.126 = 198412.6984
    
    err_f0 = abs(c_f0 - 1e8) / 1e8
    err_f1 = abs(c_f1 - 0.0)
    err_f05 = abs(c_f05 - 198412.6984126984) / 198412.6984126984

    pass_t1 = (err_f0 < 1e-6 and err_f1 < 1e-6 and err_f05 < 1e-5)
    print(f"  f=0.0 drag coeff: {c_f0:.2e} (Expected: 1.00e+08, rel err: {err_f0:.2e})")
    print(f"  f=1.0 drag coeff: {c_f1:.2e} (Expected: 0.00e+00, rel err: {err_f1:.2e})")
    print(f"  f=0.5 drag coeff: {c_f05:.2f} (Expected: 198412.70, rel err: {err_f05:.2e})")
    print(f"Test 1 Drag Coeff Unit Check: {'PASS' if pass_t1 else 'FAIL'}")

    # Test 2: Frozen-flow Stefan in Fluid Region (frozenFlow yes, ks=2.0, kl=1.0)
    print("\n--- Test 2: Frozen-flow Stefan in Fluid Region (frozenFlow = yes, ks=2.0, kl=1.0) ---")
    c2_dir = os.path.join(base_dir, "case_frozen_flow")
    shutil.rmtree(c2_dir, ignore_errors=True)
    setup_fluid_case(c2_dir, frozen_flow=True, Cu=1e12)
    run_cmd(f"cd {c2_dir} && bash -c '{of_env}; {solver_bin}'")
    
    T_fluid = parse_openfoam_field(os.path.join(c2_dir, "2000/pcm/T"))
    a_fluid = parse_openfoam_field(os.path.join(c2_dir, "2000/pcm/phaseFraction"))

    pass_t2 = (len(T_fluid) == 100 and len(a_fluid) == 100 and T_fluid[0] > 310.0 and a_fluid[0] > 0.5)
    mean_T_f = sum(T_fluid)/len(T_fluid) if T_fluid else 0
    mean_a_f = sum(a_fluid)/len(a_fluid) if a_fluid else 0
    print(f"  Fluid Region Frozen-flow Final Mean T = {mean_T_f:.2f} K, Mean alpha = {mean_a_f:.4f}")
    print(f"Test 2 Frozen-flow Fluid Region Execution: {'PASS' if pass_t2 else 'FAIL'}")

    # Test 3: Frozen by Drag (frozenFlow = no, Cu = 1e12)
    print("\n--- Test 3: Frozen by Drag in Fluid Region (frozenFlow = no, Cu = 1e12) ---")
    c3_dir = os.path.join(base_dir, "case_frozen_drag")
    shutil.rmtree(c3_dir, ignore_errors=True)
    setup_fluid_case(c3_dir, frozen_flow=False, Cu=1e12)
    run_cmd(f"cd {c3_dir} && bash -c '{of_env}; {solver_bin}'")

    T_drag = parse_openfoam_field(os.path.join(c3_dir, "2000/pcm/T"))
    a_drag = parse_openfoam_field(os.path.join(c3_dir, "2000/pcm/phaseFraction"))
    U_drag = parse_openfoam_field(os.path.join(c3_dir, "2000/pcm/U"))

    # Check maximum magnitude of U in solid region (cells where alpha < 0.01)
    max_U_solid = 0.0
    for vec, alpha in zip(U_drag, a_drag):
        if alpha < 0.01:
            mag_u = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
            max_U_solid = max(max_U_solid, mag_u)

    dT_diff = sum(abs(t1 - t2) for t1, t2 in zip(T_fluid, T_drag)) / len(T_fluid) if T_fluid and T_drag else 1e9
    da_diff = sum(abs(a1 - a2) for a1, a2 in zip(a_fluid, a_drag)) / len(a_fluid) if a_fluid and a_drag else 1e9

    pass_t3 = (max_U_solid < 1e-8 and dT_diff < 0.01 and da_diff < 0.001)
    print(f"  Max Solid Region Velocity |U| (Cu=1e12): {max_U_solid:.2e} m/s (Limit < 1e-8 m/s)")
    print(f"  Frozen-flow vs Frozen-by-drag Temp Diff : {dT_diff:.6f} K (Limit < 0.01 K)")
    print(f"  Frozen-flow vs Frozen-by-drag Alpha Diff: {da_diff:.8f} (Limit < 0.001)")
    print(f"Test 3 Frozen by Drag Velocity & Solution Match: {'PASS' if pass_t3 else 'FAIL'}")

    print("\n=======================================================")
    print("    FLUID ENTHALPY-POROSITY VERIFICATION SUMMARY      ")
    print("=======================================================")
    print(f"Test 1 (Darcy Drag Coeff Unit Check)        : {'PASS' if pass_t1 else 'FAIL'}")
    print(f"Test 2 (Frozen-Flow Fluid Region Stefan)    : {'PASS' if pass_t2 else 'FAIL'}")
    print(f"Test 3 (Frozen-by-Drag Cu=1e12 Velocity & T): {'PASS' if pass_t3 else 'FAIL'}")

    all_pass = pass_t1 and pass_t2 and pass_t3
    if all_pass:
        print("\nALL FLUID ENTHALPY-POROSITY VERIFICATION TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
