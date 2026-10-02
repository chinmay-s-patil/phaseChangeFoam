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

def setup_single_region_case(base_dir):
    os.chdir(base_dir)
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/polyMesh system/pcm system/solid2")

    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices (
    (0 0 0) (0.10 0 0) (0.10 0.01 0) (0 0.01 0) (0 0 0.01) (0.10 0 0.01) (0.10 0.01 0.01) (0 0.01 0.01)
);
blocks (
    hex (0 1 2 3 4 5 6 7) pcm (100 1 1) simpleGrading (1 1 1)
);
edges ();
boundary (
    pcm_cold { type patch; faces ((0 3 7 4)); }
    pcm_hot { type patch; faces ((1 5 6 2)); }
    emptyFaces { type empty; faces ((0 1 5 4) (3 2 6 7) (0 1 2 3) (4 5 6 7)); }
);
""")

    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 2000; deltaT 2;
writeControl runTime; writeInterval 2000; purgeWrite 0; writeFormat ascii;
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
PIMPLE { nOuterCorrectors 2; nNonOrthogonalCorrectors 0; nNonLinearCorrectors 20; nonLinearTolerance 1e-5; }
solvers { "h.*" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; } }
""")

    os.makedirs("constant", exist_ok=True)
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
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 2.0; kLiquid 0.50; }
}
""")

    os.makedirs("system/pcm", exist_ok=True)
    with open("system/pcm/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
    with open("system/pcm/fvSolution", "w") as f: f.write(open("system/fvSolution").read())

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")

    os.makedirs("0/pcm", exist_ok=True)
    with open("0/pcm/T", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField {
    pcm_cold { type fixedValue; value uniform 280.0; }
    pcm_hot { type fixedValue; value uniform 350.0; }
    emptyFaces { type empty; }
}
""")
    with open("0/pcm/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField {
    pcm_cold { type fixedValue; value uniform 560000; }
    pcm_hot { type fixedValue; value uniform 700000; }
    emptyFaces { type empty; }
}
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { ".*" { type calculated; value uniform 101325; } emptyFaces { type empty; } }
""")


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
startFrom startTime; startTime 0; stopAt endTime; endTime 2000; deltaT 2;
writeControl runTime; writeInterval 2000; purgeWrite 0; writeFormat ascii;
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

    phaseChangeDictContent = """
FoamFile { version 2.0; format ascii; class dictionary; location "constant/REGION"; object phaseChangeDict; }
active true;
phaseChange
{
    active true; phaseChangeMode EHC;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    freezing { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    hysteresis { active true; }
    density { model constant; rhoRef 1000.0; }
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 2.0; kLiquid 0.50; }
}
"""

    for region in ["pcm", "solid2"]:
        os.makedirs(f"constant/{region}", exist_ok=True)
        with open(f"constant/{region}/radiationProperties", "w") as f:
            f.write(f"FoamFile {{ version 2.0; format ascii; class dictionary; location \"constant/{region}\"; object radiationProperties; }} radiationModel none;\n")
        with open(f"constant/{region}/thermophysicalProperties", "w") as f:
            f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/{region}"; object thermophysicalProperties; }}
thermoType {{ type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }}
mixture {{ specie {{ molWeight 100; }} transport {{ kappa 0.50; }} thermodynamics {{ Cp 2000; Hf 0; }} equationOfState {{ rho 1000; }} }}
""")
        with open(f"constant/{region}/phaseChangeDict", "w") as f:
            f.write(phaseChangeDictContent.replace("REGION", region))

        os.makedirs(f"system/{region}", exist_ok=True)
        with open(f"system/{region}/fvSchemes", "w") as f: f.write(open("system/fvSchemes").read())
        with open(f"system/{region}/fvSolution", "w") as f: f.write(open("system/fvSolution").read())

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
    pcm_to_solid2 { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 280.0; Tnbr T; kappaMethod lookup; kappa kEff; }
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
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField {
    solid2_hot { type fixedValue; value uniform 350.0; }
    solid2_to_pcm { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 280.0; Tnbr T; kappaMethod lookup; kappa kEff; }
    emptyFaces { type empty; }
}
""")
    with open("0/solid2/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/solid2"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField {
    solid2_hot { type fixedValue; value uniform 700000; }
    solid2_to_pcm { type compressible::turbulentTemperatureRadCoupledMixed; value uniform 560000; Tnbr T; kappaMethod solidThermo; }
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

    # Reference: Single-region
    print("\n--- Running Single-Region Reference Case ---")
    single_dir = os.path.join(base_dir, "case_single")
    os.makedirs(single_dir, exist_ok=True)
    setup_single_region_case(single_dir)
    run_cmd(f"cd {single_dir} && bash -c '{of_env}; {solver_bin}'")

    T_single = parse_openfoam_field(os.path.join(single_dir, "2000/pcm/T"), num_cells=100)
    a_single = parse_openfoam_field(os.path.join(single_dir, "2000/pcm/phaseFraction"), num_cells=100)

    # Test 1: nOuterCorr = 2 (Valid coupled run)
    print("\n--- Test 1: Coupled Solve with nOuterCorr = 2 ---")
    case1_dir = os.path.join(base_dir, "case_coupled_nOuter2")
    os.makedirs(case1_dir, exist_ok=True)
    setup_coupled_case(case1_dir, nOuterCorr=2)
    res1 = run_cmd(f"cd {case1_dir} && bash -c '{of_env}; {solver_bin}'")
    
    T_pcm_left = parse_openfoam_field(os.path.join(case1_dir, "2000/pcm/T"), num_cells=50)
    a_pcm_left = parse_openfoam_field(os.path.join(case1_dir, "2000/pcm/phaseFraction"), num_cells=50)
    T_pcm_right = parse_openfoam_field(os.path.join(case1_dir, "2000/solid2/T"), num_cells=50)
    a_pcm_right = parse_openfoam_field(os.path.join(case1_dir, "2000/solid2/phaseFraction"), num_cells=50)

    # The coupled case has two regions of 50 cells each
    if T_pcm_left and T_pcm_right and T_single:
        T_coupled = T_pcm_left + T_pcm_right  # concatenate
        a_coupled = a_pcm_left + a_pcm_right
        
        max_dT = max(abs(Tc - Ts) for Tc, Ts in zip(T_coupled, T_single))
        max_da = max(abs(ac - as_) for ac, as_ in zip(a_coupled, a_single))
        
        # Interface temperature continuity
        T_interface_left = T_pcm_left[-1]  # last cell of left region
        T_interface_right = T_pcm_right[0]  # first cell of right region  
        interface_dT = abs(T_interface_left - T_interface_right)
        
        mean_a = sum(a_coupled) / len(a_coupled)
        
        pass_t1 = (max_dT < 0.05 and max_da < 0.005 and interface_dT < 0.5 and mean_a > 0.05)
        
        print(f"Max T diff compared to single-region: {max_dT:.4f} K")
        print(f"Max alpha diff compared to single-region: {max_da:.6f}")
        print(f"Interface T discontinuity: {interface_dT:.4f} K")
        print(f"Mean alpha = {mean_a:.4f}")
    else:
        pass_t1 = False
        print("Failed to read fields.")

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
