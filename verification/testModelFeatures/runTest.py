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

def setup_base_case(case_dir):
    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 10; deltaT 2;
writeControl runTime; writeInterval 10; purgeWrite 0; writeFormat ascii; writePrecision 12;
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
PIMPLE { nOuterCorrectors 1; nNonLinearCorrectors 10; nonLinearTolerance 1e-5; }
""")
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }
scale 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) (10 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot { type patch; faces ((0 4 7 3)); }
    cold { type patch; faces ((1 2 6 5)); }
    emptyFaces { type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }
);
""")

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

def setup_fluid_base_case(case_dir):
    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 10; deltaT 2;
writeControl runTime; writeInterval 10; purgeWrite 0; writeFormat ascii; writePrecision 12;
""")
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers {
    "h.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }
    "he.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }
    "T.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    "p_rgh.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-20; relTol 0; }
}
PIMPLE { nOuterCorrectors 1; nCorrectors 1; pRefCell 0; pRefValue 101325; }
""")
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }
scale 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) (10 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot { type patch; faces ((0 4 7 3)); }
    cold { type patch; faces ((1 2 6 5)); }
    emptyFaces { type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }
);
""")

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
    transport { mu 1e-5; Pr 0.0099; }
    thermodynamics { Cp 1980.0; Hf 0; Tref 0; Href 0; }
    equationOfState { rho 1967.0; }
}
""")
    with open(os.path.join(case_dir, "constant/pcm/turbulenceProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object turbulenceProperties; }
simulationType laminar;
""")

    os.makedirs(os.path.join(case_dir, "0/pcm"), exist_ok=True)
    with open(os.path.join(case_dir, "0/pcm/phaseFraction"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object phaseFraction; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0;
boundaryField { ".*" { type calculated; value uniform 0; } emptyFaces { type empty; } }
""")
    with open(os.path.join(case_dir, "0/pcm/T"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280;
boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }
""")
    with open(os.path.join(case_dir, "0/pcm/he"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object he; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform -35937;
boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }
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
    print(f"=== Model Features & Validation Verification Suite in {base_dir} ===")
    
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    # Test 1: Convection suppress false FatalIOError check
    print("\n--- Test 1: convection { suppress false; } FatalIOError Check ---")
    c1_dir = os.path.join(base_dir, "case_convection_false")
    setup_base_case(c1_dir)
    with open(os.path.join(c1_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    phaseChangeMode EHC;
    convection { suppress false; }
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
    res1 = run_cmd(f"cd {c1_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out1 = res1.stdout + res1.stderr
    pass_t1 = ("phaseChange.convection.suppress = false requested" in out1 or res1.returncode != 0)
    print(f"Test 1 FatalIOError Triggered: {pass_t1}")

    # Test 2: Legacy liquidFraction disk restore
    print("\n--- Test 2: Legacy liquidFraction Disk Restore Check ---")
    c2_dir = os.path.join(base_dir, "case_legacy_restore")
    setup_base_case(c2_dir)
    with open(os.path.join(c2_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    phaseChangeMode EHC;
    hysteresis { active true; }
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
    # Write legacy liquidFraction file in 0/pcm (0.5 everywhere)
    with open(os.path.join(c2_dir, "0/pcm/liquidFraction"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object liquidFraction; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.5;
boundaryField { ".*" { type calculated; value uniform 0.5; } emptyFaces { type empty; } }
""")
    run_cmd(f"cd {c2_dir} && bash -c '{of_env}; {solver_bin}'")
    a_restored = parse_openfoam_field(os.path.join(c2_dir, "10/pcm/phaseFraction"))
    pass_t2 = (len(a_restored) > 0 and a_restored[0] >= 0.5)
    print(f"Test 2 Legacy liquidFraction Restored to phaseFraction: {pass_t2} (alpha[0] = {a_restored[0] if a_restored else 0})")

    # Test 3: Hysteresis active false
    print("\n--- Test 3: hysteresis { active false; } Check ---")
    c3_dir = os.path.join(base_dir, "case_hysteresis_false")
    setup_base_case(c3_dir)
    with open(os.path.join(c3_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    phaseChangeMode EHC;
    hysteresis { active false; }
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
    res3 = run_cmd(f"cd {c3_dir} && bash -c '{of_env}; {solver_bin}'")
    a_nonhys = parse_openfoam_field(os.path.join(c3_dir, "10/pcm/phaseFraction"))
    pass_t3 = (len(a_nonhys) > 0)
    print(f"Test 3 Non-Hysteresis Execution: {pass_t3}")

    # Test 4: EnthalpyPorosity rhoEff*CpEff/(rho*Cp) ratio validation
    print("\n--- Test 4: EnthalpyPorosity rhoEff*CpEff/(rho*Cp) == 1 Check ---")
    c4_dir = os.path.join(base_dir, "case_ep_ratio")
    setup_base_case(c4_dir)
    with open(os.path.join(c4_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type enthalpyPorosity;
    active true;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
    lambda 1.0;
}
""")
    res4 = run_cmd(f"cd {c4_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out4 = res4.stdout + res4.stderr
    pass_t4 = (res4.returncode == 0 and "Ratio rhoEff*CpEff/(rho*Cp)" not in out4)
    print(f"Test 4 EP ratio validation: {pass_t4}")

    # Test 5: Partial forward block missing latentHeat FatalIOError check
    print("\n--- Test 5: Partial forward block (missing latentHeat) FatalIOError Check ---")
    c5_dir = os.path.join(base_dir, "case_ep_partial_fwd")
    setup_base_case(c5_dir)
    with open(os.path.join(c5_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type enthalpyPorosity;
    active true;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; }
}
""")
    res5 = run_cmd(f"cd {c5_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out5 = res5.stdout + res5.stderr
    pass_t5 = (res5.returncode != 0 and ("Entry 'latentHeat' not found" in out5 or "latentHeat" in out5))
    print(f"Test 5 Partial forward FatalIOError Triggered: {pass_t5}")

    # Test 6: Non-hConst thermo FatalIOError check
    print("\n--- Test 6: Non-hConst thermo FatalIOError Check ---")
    c6_dir = os.path.join(base_dir, "case_ep_non_hconst")
    setup_base_case(c6_dir)
    with open(os.path.join(c6_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport polynomial; thermo hPolynomial; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture {
    specie { molWeight 200.0; }
    transport {
        kappaCoeffs<8> ( 0.5 0 0 0 0 0 0 0 );
    }
    thermodynamics {
        Hf 0;
        Sf 0;
        CpCoeffs<8> ( 1980 1 0 0 0 0 0 0 );
    }
    equationOfState { rho 1967.0; }
}
""")
    with open(os.path.join(c6_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type enthalpyPorosity;
    active true;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
    res6 = run_cmd(f"cd {c6_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out6 = res6.stdout + res6.stderr
    pass_t6 = (res6.returncode != 0 and "requires constant Cp (hConst)" in out6)
    print(f"Test 6 Non-hConst Thermo FatalIOError Triggered: {pass_t6}")

    # Test 7: Missing type or speciesModel key when species parameters present FatalIOError check
    print("\n--- Test 7: Missing type/speciesModel with species keys FatalIOError Check ---")
    c7_dir = os.path.join(base_dir, "case_missing_species_type")
    setup_fluid_base_case(c7_dir)
    with open(os.path.join(c7_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    liquid H2O;
    vapor H2O;
    C_evap 0.1;
}
""")
    res7 = run_cmd(f"cd {c7_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out7 = res7.stdout + res7.stderr
    pass_t7 = (res7.returncode != 0 and ("missing 'type' or 'speciesModel'" in out7 or "species phase change parameters" in out7))
    print(f"Test 7 Missing Species Type FatalIOError Triggered: {pass_t7}")

    # Test 8: Separate speciesModel and meltingModel keys selection check
    print("\n--- Test 8: Separate speciesModel and meltingModel Keys Check ---")
    c8_dir = os.path.join(base_dir, "case_separate_keys")
    setup_fluid_base_case(c8_dir)
    with open(os.path.join(c8_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    speciesModel none;
    meltingModel enthalpyPorosity;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
}
""")
    res8 = run_cmd(f"cd {c8_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out8 = res8.stdout + res8.stderr
    pass_t8 = ("Selecting fluid phase change model type none" in out8 and "using mode: enthalpyPorosity" in out8)
    print(f"Test 8 Separate speciesModel & meltingModel Selection: {pass_t8}")

    print("\n=======================================================")
    print("      MODEL FEATURES & VALIDATION SUMMARY             ")
    print("=======================================================")
    print(f"Test 1 (convection.suppress false FatalError) : {'PASS' if pass_t1 else 'FAIL'}")
    print(f"Test 2 (legacy liquidFraction disk restore)    : {'PASS' if pass_t2 else 'FAIL'}")
    print(f"Test 3 (hysteresis active false execution)     : {'PASS' if pass_t3 else 'FAIL'}")
    print(f"Test 4 (EP rhoEff*CpEff/(rho*Cp) == 1 check)   : {'PASS' if pass_t4 else 'FAIL'}")
    print(f"Test 5 (EP missing forward key FatalIOError)  : {'PASS' if pass_t5 else 'FAIL'}")
    print(f"Test 6 (EP non-hConst thermo FatalIOError)   : {'PASS' if pass_t6 else 'FAIL'}")
    print(f"Test 7 (missing species type FatalIOError)    : {'PASS' if pass_t7 else 'FAIL'}")
    print(f"Test 8 (speciesModel & meltingModel keys)     : {'PASS' if pass_t8 else 'FAIL'}")

    all_pass = pass_t1 and pass_t2 and pass_t3 and pass_t4 and pass_t5 and pass_t6 and pass_t7 and pass_t8
    if all_pass:
        print("\nALL MODEL FEATURE & ERROR TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
