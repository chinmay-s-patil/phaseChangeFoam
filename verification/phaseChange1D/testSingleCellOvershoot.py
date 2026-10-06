#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil

def find_solver():
    # 1. Check environment variable
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        path = os.path.join(user_appbin, "phaseChangeMultiRegionFoam")
        if os.path.exists(path):
            return path
    
    # 2. Check PATH
    path = shutil.which("phaseChangeMultiRegionFoam")
    if path:
        return path

    # 3. Fallback search in home platforms
    home = os.environ.get("HOME", "/home/lavender")
    matches = glob_matches = [
        os.path.join(home, "OpenFOAM", f"lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam")
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

def main():
    test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "singleCellOvershoot_run")
    print(f"=== Single-Cell Overshoot Energy Balance Test in {test_dir} ===")
    
    run_cmd(f"rm -rf {test_dir}")
    os.makedirs(test_dir, exist_ok=True)
    os.chdir(test_dir)
    
    # 1. Create 1-cell blockMeshDict (all zeroGradient / adiabatic boundaries)
    os.makedirs("system", exist_ok=True)
    with open("system/blockMeshDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
vertices ((0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 1) (1 0 1) (1 1 1) (0 1 1));
blocks ( hex (0 1 2 3 4 5 6 7) pcm (1 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    walls { type zeroGradient; faces ((0 3 7 4) (1 5 6 2) (0 1 5 4) (3 2 6 7) (0 1 2 3) (4 5 6 7)); }
);
""")

    with open("system/controlDict", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 100; deltaT 100;
writeControl runTime; writeInterval 100; purgeWrite 0; writeFormat ascii; writePrecision 12;
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
PIMPLE
{
    nNonOrthogonalCorrectors 0;
    nNonLinearCorrectors 25;
    nonLinearTolerance 1e-6;
}
solvers { "h.*" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; } }
""")

    os.makedirs("system/pcm", exist_ok=True)
    with open("system/pcm/fvSchemes", "w") as f:
        f.write(open("system/fvSchemes").read())
    with open("system/pcm/fvSolution", "w") as f:
        f.write(open("system/fvSolution").read())

    # constant & fvOptions
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
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }
active true;
phaseChange
{
    active true; phaseChangeMode EHC;
    melting { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    freezing { T_lowerBound 300.0; T_upperBound 310.0; latentHeat 100000.0; }
    hysteresis { active false; }
    density { model linear; rhoRef 1000.0; rhoSolid 1000.0; rhoLiquid 1000.0; }
    thermophysical { mode custom; CpSolid 2000.0; CpLiquid 2000.0; kSolid 1.0; kLiquid 1.0; }
}
""")

    with open("system/pcm/fvOptions", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "system/pcm"; object fvOptions; }
heatSource
{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   all;
    volumeMode      specific;
    injectionRateSuSp
    {
        h           (2.4e6 0);
    }
}
""")

    os.makedirs("0/pcm", exist_ok=True)
    with open("0/pcm/T", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 280.0;
boundaryField { walls { type zeroGradient; } }
""")
    with open("0/pcm/h", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 560000;
boundaryField { walls { type zeroGradient; } }
""")
    with open("0/pcm/p", "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/pcm"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { walls { type calculated; value uniform 101325; } }
""")

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()

    print(f"\n--- Running Single-Cell Overshoot Test (dt = 100 s in 1 STEP) ---")
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    out = run_cmd(f"bash -c '{of_env}; {solver_bin}'")

    T_sim = parse_openfoam_field("100/pcm/T")[0]
    alpha_sim = parse_openfoam_field("100/pcm/phaseFraction")[0]

    T_analytical = 350.00
    alpha_analytical = 1.0000
    err = abs(T_sim - T_analytical)

    print("\n=======================================================")
    print("      SINGLE-CELL OVERSHOOT TEST RESULTS (1 STEP)      ")
    print("=======================================================")
    print(f"Analytical Exact T_final  : {T_analytical:.6f} K")
    print(f"Simulated T_final (dt=100s): {T_sim:.6f} K")
    print(f"Absolute Temperature Error: {err:.6f} K")
    print(f"Analytical alphaL         : {alpha_analytical:.4f}")
    print(f"Simulated alphaL          : {alpha_sim:.4f}")
    print("-------------------------------------------------------")

    if err < 0.001:
        print("\nSTATUS: SINGLE-CELL OVERSHOOT TEST PASSED!")
        print("Convergence-driven non-linear iterations achieved exact < 0.001 K energy conservation in 1 step!")
    else:
        print(f"\nSTATUS: TEST FAILED (Error {err:.6f} K >= 0.001 K)")
        sys.exit(1)

if __name__ == "__main__":
    main()
