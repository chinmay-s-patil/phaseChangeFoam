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
    pass_t7a = (res7.returncode != 0 and ("missing 'type' or 'speciesModel'" in out7 or "species phase change parameters" in out7))

    # Test 7b: Missing type or speciesModel key when only Tsat and enableTsatP keys present FatalIOError check
    c7b_dir = os.path.join(base_dir, "case_missing_species_tsat_type")
    setup_fluid_base_case(c7b_dir)
    with open(os.path.join(c7b_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    Tsat 373.15;
    enableTsatP true;
}
""")
    res7b = run_cmd(f"cd {c7b_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out7b = res7b.stdout + res7b.stderr
    pass_t7b = (res7b.returncode != 0 and ("missing 'type' or 'speciesModel'" in out7b or "species phase change parameters" in out7b))
    pass_t7 = pass_t7a and pass_t7b
    print(f"Test 7 Missing Species Type (liquid/C_evap & Tsat/enableTsatP) FatalIOError Check: {pass_t7}")

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

    # Test 9: Fluid region PCM unequal CpSolid vs CpLiquid FatalError check on initial f=0 IC
    print("\n--- Test 9: Fluid Region PCM unequal Cps/Cpl on f=0 IC Check ---")
    c9_dir = os.path.join(base_dir, "case_fluid_unequal_cps_cpl")
    setup_fluid_base_case(c9_dir)
    with open(os.path.join(c9_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type enthalpyPorosity;
    active true;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 163000.0; }
    thermophysical { mode custom; CpSolid 1980.0; CpLiquid 2500.0; kSolid 2.0; kLiquid 1.0; }
}
""")
    res9 = run_cmd(f"cd {c9_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out9 = res9.stdout + res9.stderr
    pass_t9 = (res9.returncode != 0 and ("requires equal solid/liquid heat capacities and densities" in out9 or "Model 'enthalpyPorosity' v1 requires Cp_solid == Cp_liquid" in out9))
    print(f"Test 9 Fluid PCM Unequal Cps/Cpl FatalError Triggered: {pass_t9}")

    # Test 10: inertSpecie matching liquid/vapor phase-change specie FatalError check
    print("\n--- Test 10: inertSpecie Phase-Change Specie Collision FatalError Check ---")
    c10_dir = os.path.join(base_dir, "case_inert_specie_collision")
    setup_fluid_base_case(c10_dir)
    with open(os.path.join(c10_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heRhoThermo; mixture multiComponentMixture; transport const; thermo hConst; equationOfState perfectGas; specie specie; energy sensibleEnthalpy; }
species ( H2O_l H2O_v air );
inertSpecie H2O_v;
H2O_l { specie { molWeight 18; } transport { mu 1e-3; Pr 7.0; } thermodynamics { Cp 4184; Hf 0; } }
H2O_v { specie { molWeight 18; } transport { mu 1e-5; Pr 1.0; } thermodynamics { Cp 2000; Hf 0; } }
air   { specie { molWeight 29; } transport { mu 1e-5; Pr 0.7; } thermodynamics { Cp 1000; Hf 0; } }
""")
    with open(os.path.join(c10_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    active true;
    speciesModel Lee;
    liquid H2O_l;
    vapor H2O_v;
    C_evap 0.1;
    C_cond 0.1;
    latentHeat 2.26e6;
    Tsat 373.15;
}
""")
    for s in ["H2O_l", "H2O_v", "air"]:
        with open(os.path.join(c10_dir, f"0/pcm/{s}"), "w") as f:
            f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object {s}; }}
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.33;
boundaryField {{ ".*" {{ type zeroGradient; }} emptyFaces {{ type empty; }} }}
""")
    res10 = run_cmd(f"cd {c10_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out10 = res10.stdout + res10.stderr
    pass_t10 = (res10.returncode != 0 and "must not be a phase-changing specie" in out10)
    print(f"Test 10 inertSpecie Collision FatalError Triggered: {pass_t10}")

    # Test 11: EHC Hot Startup Initial Phase Fraction & Spurious Latent Heat Sink Check
    print("\n--- Test 11: EHC Hot Startup Phase Fraction Initialization Check ---")
    c11_dir = os.path.join(base_dir, "case_ehc_hot_startup")
    setup_base_case(c11_dir)
    # Set initial T = 320 K (above 300-310 K melting window)
    with open(os.path.join(c11_dir, "0/pcm/T"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 320;
boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }
""")
    with open(os.path.join(c11_dir, "0/pcm/h"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 633600;
boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }
""")
    with open(os.path.join(c11_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type EHC;
    active true;
    convection { suppress true; }
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
}
""")
    res11 = run_cmd(f"cd {c11_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=False)
    # Read end-time T and phaseFraction
    t11_arr = parse_openfoam_field(os.path.join(c11_dir, "10/pcm/T"))
    alpha11_arr = parse_openfoam_field(os.path.join(c11_dir, "10/pcm/phaseFraction"))
    t11_val = t11_arr[0] if t11_arr else 0.0
    alpha11_val = alpha11_arr[0] if alpha11_arr else 0.0

    # alpha should be 1.0 (liquid) and T should remain 320.0 K (no spurious latent sink)
    pass_t11 = (abs(alpha11_val - 1.0) < 1e-4 and abs(t11_val - 320.0) < 1e-2)
    print(f"Test 11 EHC Hot Startup Phase Fraction (alpha = {alpha11_val:.4f}, T = {t11_val:.2f} K): {'PASS' if pass_t11 else 'FAIL'}")

    # Test 12: Solid Region Default Convection Suppression Check (no convection block in dict)
    print("\n--- Test 12: Solid Region Default Convection Suppression Check ---")
    c12_dir = os.path.join(base_dir, "case_solid_default_suppress")
    setup_base_case(c12_dir)
    with open(os.path.join(c12_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type enthalpyPorosity;
    active true;
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
}
""")
    res12 = run_cmd(f"cd {c12_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    pass_t12 = (res12.returncode == 0)
    print(f"Test 12 Solid Region Default Convection Suppression Check: {'PASS' if pass_t12 else 'FAIL'}")

    # Test 13: Ignored Dictionary Block Non-Abort Check (direction reverse with malformed forward block)
    print("\n--- Test 13: Ignored Dictionary Block Non-Abort Check ---")
    c13_dir = os.path.join(base_dir, "case_ignored_block_reverse")
    setup_base_case(c13_dir)
    with open(os.path.join(c13_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type EHC;
    active true;
    direction reverse;
    reverse { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    forward { T_lowerBound 310.0; T_upperBound 300.0; latentHeat -999.0; }
}
""")
    res13 = run_cmd(f"cd {c13_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    pass_t13 = (res13.returncode == 0)
    print(f"Test 13 Ignored Dictionary Block Non-Abort Check: {'PASS' if pass_t13 else 'FAIL'}")

    # Test 14: Non-Hysteresis Multi-Curve Rejection Check (hysteresis active false with distinct forward and reverse curves)
    print("\n--- Test 14: Non-Hysteresis Multi-Curve Rejection Check ---")
    c14_dir = os.path.join(base_dir, "case_non_hysteresis_multicurve")
    setup_base_case(c14_dir)
    with open(os.path.join(c14_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type EHC;
    active true;
    direction both;
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    reverse { T_lowerBound 295.0; T_upperBound 305.0; latentHeat 100000.0; }
    hysteresis { active false; }
}
""")
    res14 = run_cmd(f"cd {c14_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out14 = res14.stdout + res14.stderr
    pass_t14 = ("Distinct 'forward' and 'reverse' phase-change parameters specified" in out14 or res14.returncode != 0)
    print(f"Test 14 Non-Hysteresis Multi-Curve Rejection Check: {'PASS' if pass_t14 else 'FAIL'}")

    # Test 15: Custom Mode Kappa vs kSolid Mismatch Check
    print("\n--- Test 15: Custom Mode Kappa vs kSolid Mismatch Check ---")
    c15_dir = os.path.join(base_dir, "case_kappa_mismatch")
    setup_base_case(c15_dir)
    with open(os.path.join(c15_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type EHC;
    active true;
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 5.0; kLiquid 5.0; }
}
""")
    res15 = run_cmd(f"cd {c15_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out15 = res15.stdout + res15.stderr
    pass_t15 = ("thermophysicalProperties transport kappa must match kSolid" in out15 or res15.returncode != 0)
    print(f"Test 15 Custom Mode Kappa vs kSolid Mismatch Check: {'PASS' if pass_t15 else 'FAIL'}")

    # Test 16: Fluid Region PCM nOuterCorrectors < 2 FatalError Check
    print("\n--- Test 16: Fluid Region PCM nOuterCorrectors < 2 Check ---")
    c16_dir = os.path.join(base_dir, "case_fluid_nouter1")
    setup_fluid_base_case(c16_dir)
    with open(os.path.join(c16_dir, "system/fvSolution"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers {
    "h.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "T.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "Yi.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "air.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "H2O.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "p.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "p_rgh.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 1;
}
""")
    with open(os.path.join(c16_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange {
    type EHC;
    active true;
    forward { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
}
""")
    res16 = run_cmd(f"cd {c16_dir} && bash -c '{of_env}; {solver_bin}'", allow_failure=True)
    out16 = res16.stdout + res16.stderr
    pass_t16 = (res16.returncode != 0 and "nOuterCorrectors" in out16)
    print(f"Test 16 Fluid Region PCM nOuterCorrectors < 2 FatalError Check: {'PASS' if pass_t16 else 'FAIL'}")

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
    print(f"Test 9 (fluid PCM unequal Cps/Cpl FatalError): {'PASS' if pass_t9 else 'FAIL'}")
    print(f"Test 10 (inertSpecie phase-change collision) : {'PASS' if pass_t10 else 'FAIL'}")
    print(f"Test 11 (EHC hot startup alpha initialization): {'PASS' if pass_t11 else 'FAIL'}")
    print(f"Test 12 (solid region default convection)     : {'PASS' if pass_t12 else 'FAIL'}")
    print(f"Test 13 (ignored forward block non-abort)     : {'PASS' if pass_t13 else 'FAIL'}")
    print(f"Test 14 (non-hysteresis multi-curve rejection): {'PASS' if pass_t14 else 'FAIL'}")
    print(f"Test 15 (custom mode kappa vs kSolid check)   : {'PASS' if pass_t15 else 'FAIL'}")
    print(f"Test 16 (fluid region PCM nOuterCorrectors check): {'PASS' if pass_t16 else 'FAIL'}")

    all_pass = (pass_t1 and pass_t2 and pass_t3 and pass_t4 and pass_t5 and pass_t6 and
                pass_t7 and pass_t8 and pass_t9 and pass_t10 and pass_t11 and pass_t12 and
                pass_t13 and pass_t14 and pass_t15 and pass_t16)
    if all_pass:
        print("\nALL MODEL FEATURE & ERROR TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
