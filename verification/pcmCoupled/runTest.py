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

def run_cmd(cmd, cwd=None, allow_failure=False):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0 and not allow_failure:
        print(f"Error executing command: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res

def parse_openfoam_field(file_path, num_cells=50):
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

def setup_coupled_case(base_dir, nOuterCorr=2):
    os.chdir(base_dir)
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/solid2/polyMesh constant/polyMesh system/pcm system/solid2")

    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices (
    (0 0 0) (0.05 0 0) (0.05 0.01 0) (0 0.01 0) (0 0 0.01) (0.05 0 0.01) (0.05 0.01 0.01) (0 0.01 0.01)
    (0.05 0 0) (0.10 0 0) (0.10 0.01 0) (0.05 0.01 0) (0.05 0 0.01) (0.10 0 0.01) (0.10 0.01 0.01) (0.05 0.01 0.01)
);
blocks (
    hex (0 1 2 3 4 5 6 7) pcm (50 1 1) simpleGrading (1 1 1)
    hex (8 9 10 11 12 13 14 15) solid2 (50 1 1) simpleGrading (1 1 1)
);
edges ();
boundary (
    pcm_cold { type patch; faces ((0 3 7 4)); }
    pcm_to_solid2 { type mappedWall; sampleMode nearestPatchFace; sampleRegion solid2; samplePatch solid2_to_pcm; faces ((1 5 6 2)); }
    solid2_to_pcm { type mappedWall; sampleMode nearestPatchFace; sampleRegion pcm; samplePatch pcm_to_solid2; faces ((8 11 15 12)); }
    solid2_hot { type patch; faces ((9 13 14 10)); }
    emptyFaces { type empty; faces ((0 1 5 4) (3 2 6 7) (0 1 2 3) (4 5 6 7) (8 9 13 12) (11 10 14 15) (8 9 10 11) (12 13 14 15)); }
);
""")

    with open("system/controlDict", "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 10; deltaT 1;
writeControl runTime; writeInterval 10; purgeWrite 0; writeFormat ascii;
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
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }}
PIMPLE {{ nOuterCorrectors {nOuterCorr}; nNonOrthogonalCorrectors 0; nNonLinearCorrectors 20; nonLinearTolerance 1e-5; }}
solvers {{ "h.*" {{ solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }} }}
""")

    os.makedirs("constant", exist_ok=True)
    with open("constant/regionProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid () solid (pcm solid2) porousFluid () porousSolid () );
""")
    with open("constant/g", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0]; value (0 0 0);
""")

    # Region pcm
    os.makedirs("constant/pcm", exist_ok=True)
    with open("constant/pcm/radiationProperties", "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; location \"constant/pcm\"; object radiationProperties; } radiationModel none;\n")
    with open("constant/pcm/thermophysicalProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture { specie { molWeight 100; } transport { kappa 0.50; } thermodynamics { Cp 2000; Hf 0; } equationOfState { rho 1000; } }
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
    hysteresis { active true; }
    density { model constant; rhoRef 1000.0; }
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 0.50; kLiquid 0.50; }
}
""")

    # Region solid2
    os.makedirs("constant/solid2", exist_ok=True)
    with open("constant/solid2/radiationProperties", "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; location \"constant/solid2\"; object radiationProperties; } radiationModel none;\n")
    with open("constant/solid2/thermophysicalProperties", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/solid2"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture { specie { molWeight 100; } transport { kappa 0.50; } thermodynamics { Cp 2000; Hf 0; } equationOfState { rho 1000; } }
""")
    with open("constant/solid2/phaseChangeDict", "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; location \"constant/solid2\"; object phaseChangeDict; } active false;\n")

    os.makedirs("system/pcm", exist_ok=True)
    os.makedirs("system/solid2", exist_ok=True)
    with open("system/pcm/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/pcm/fvSolution", "w") as f: f.write(open("system/fvSolution").read())
    with open("system/solid2/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/solid2/fvSolution", "w") as f: f.write(open("system/fvSolution").read())

    # Create meshes
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd(f"bash -c '{of_env}; splitMeshRegions -cellZones -overwrite'")

    # Set 0 fields for pcm
    os.makedirs("0/pcm", exist_ok=True)
    with open("0/pcm/T", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField {
    pcm_cold { type fixedValue; value uniform 280.0; }
    pcm_to_solid2 { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 280.0; Tnbr T; kappaMethod lookup; kappa kPCM; }
    emptyFaces { type empty; }
}
""")
    with open("0/pcm/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField {
    pcm_cold { type fixedValue; value uniform 560000; }
    pcm_to_solid2 { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 560000; Tnbr T; kappaMethod solidThermo; }
    emptyFaces { type empty; }
}
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } emptyFaces { type empty; } }
""")

    # Set 0 fields for solid2
    os.makedirs("0/solid2", exist_ok=True)
    with open("0/solid2/T", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/solid2"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 350.0;
boundaryField {
    solid2_hot { type fixedValue; value uniform 350.0; }
    solid2_to_pcm { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 350.0; Tnbr T; kappaMethod solidThermo; }
    emptyFaces { type empty; }
}
""")
    with open("0/solid2/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/solid2"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 700000;
boundaryField {
    solid2_hot { type fixedValue; value uniform 700000; }
    solid2_to_pcm { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 700000; Tnbr T; kappaMethod solidThermo; }
    emptyFaces { type empty; }
}
""")
    with open("0/solid2/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/solid2"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } emptyFaces { type empty; } }
""")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== Coupled Multi-Region PCM Phase Change Verification Test in {base_dir} ===")

    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    # Test 1: nOuterCorr = 2 (Valid coupled run)
    print("\n--- Test 1: Coupled Solve with nOuterCorr = 2 ---")
    case1_dir = os.path.join(base_dir, "case_coupled_nOuter2")
    os.makedirs(case1_dir, exist_ok=True)
    setup_coupled_case(case1_dir, nOuterCorr=2)

    res1 = run_cmd(f"cd {case1_dir} && bash -c '{of_env}; {solver_bin}'")
    T_pcm = parse_openfoam_field(os.path.join(case1_dir, "10/pcm/T"))
    a_pcm = parse_openfoam_field(os.path.join(case1_dir, "10/pcm/liquidFraction"))
    mean_T = sum(T_pcm) / len(T_pcm) if T_pcm else 0.0
    mean_a = sum(a_pcm) / len(a_pcm) if a_pcm else 0.0

    print(f"nOuterCorr = 2 Run Completed: Mean pcm T = {mean_T:.2f} K, Mean alphaL = {mean_a:.4f}")
    pass_t1 = (mean_T > 280.0 and mean_T < 350.0 and mean_a > 0.0)
    print(f"Test 1 Status: {'PASSED' if pass_t1 else 'FAILED'}")

    # Test 2: nOuterCorr = 1 (Should trigger FatalError)
    print("\n--- Test 2: Coupled Solve with nOuterCorr = 1 (Expecting FatalError) ---")
    case2_dir = os.path.join(base_dir, "case_coupled_nOuter1")
    os.makedirs(case2_dir, exist_ok=True)
    setup_coupled_case(case2_dir, nOuterCorr=1)

    res2 = run_cmd(f"cd {case2_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    output_text = res2.stdout + res2.stderr
    pass_t2 = ("Phase change is active in a multi-region coupled domain, but nOuterCorrectors = 1 < 2" in output_text)

    print(f"nOuterCorr = 1 FatalError Triggered: {pass_t2}")
    print(f"Test 2 Status: {'PASSED' if pass_t2 else 'FAILED'}")

    print("\n=======================================================")
    print("      COUPLED MULTI-REGION VERIFICATION RESULTS        ")
    print("=======================================================")
    print(f"Test 1 (nOuterCorr = 2 Convergence) : {'PASSED' if pass_t1 else 'FAILED'}")
    print(f"Test 2 (nOuterCorr = 1 FatalError)  : {'PASSED' if pass_t2 else 'FAILED'}")
    print("-------------------------------------------------------")

    if pass_t1 and pass_t2:
        print("\nALL COUPLED MULTI-REGION VERIFICATION TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nCOUPLED TEST FAILED!")
        sys.exit(1)

if __name__ == "__main__":
    main()
