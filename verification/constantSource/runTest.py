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

def setup_case(case_dir):
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
    "H2O.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "p.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "p_rgh.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; }
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 2;
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
    H2O
);

inertSpecie air;

air
{
    specie { molWeight 28.96; }
    thermodynamics { Cp 1000; Hf 0; }
    transport { mu 1.8e-5; Pr 0.7; }
}

H2O
{
    specie { molWeight 18.015; }
    thermodynamics { Cp 1000; Hf 0; }
    transport { mu 1.8e-5; Pr 0.7; }
}
""")

    with open(os.path.join(air_const, "phaseChangeDict"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/air"; object phaseChangeDict; }
phaseChange
{
    active true;
    type constantSource;
    energySource -1e5;
    speciesName H2O;
    speciesSource 1e-3;
    massSource 1e-3;
}
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

    with open(os.path.join(air_zero, "T"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 300;
boundaryField { walls { type zeroGradient; } }
""")

    with open(os.path.join(air_zero, "h"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object h; }
dimensions [0 2 -2 0 0 0 0]; internalField uniform 1850;
boundaryField { walls { type zeroGradient; } }
""")

    with open(os.path.join(air_zero, "p"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { walls { type calculated; value uniform 101325; } }
""")

    with open(os.path.join(air_zero, "p_rgh"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField { walls { type fixedFluxPressure; value uniform 101325; } }
""")

    with open(os.path.join(air_zero, "U"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volVectorField; location "0/air"; object U; }
dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
boundaryField { walls { type fixedValue; value uniform (0 0 0); } }
""")

    with open(os.path.join(air_zero, "rho"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object rho; }
dimensions [1 -3 0 0 0 0 0]; internalField uniform 1.1768;
boundaryField { walls { type calculated; value uniform 1.1768; } }
""")

    with open(os.path.join(air_zero, "Ydefault"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object Ydefault; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.0;
boundaryField { walls { type zeroGradient; } }
""")

    with open(os.path.join(air_zero, "air"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object air; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 1.0;
boundaryField { walls { type zeroGradient; } }
""")

    with open(os.path.join(air_zero, "H2O"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object H2O; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.0;
boundaryField { walls { type zeroGradient; } }
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
    case_dir = os.path.join(base_dir, "case_constantSource")
    print(f"=== constantSource Fluid Phase Change Plumbing Verification in {case_dir} ===")
    
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    
    setup_case(case_dir)

    print("\n--- Running phaseChangeMultiRegionFoam ---")
    res = run_cmd(f"cd {case_dir} && bash -c '{of_env}; {solver_bin}'", cwd=case_dir, allow_failure=False)
    print(res.stdout[-1000:])

    # Extract numerical fields
    T_init = parse_openfoam_field(os.path.join(case_dir, "0/air/T"))[0]
    T_final = parse_openfoam_field(os.path.join(case_dir, "1/air/T"))[0]

    Y_init = parse_openfoam_field(os.path.join(case_dir, "0/air/H2O"))[0]
    Y_final = parse_openfoam_field(os.path.join(case_dir, "1/air/H2O"))[0]

    rho_init = parse_openfoam_field(os.path.join(case_dir, "0/air/rho"))[0]

    # Model parameters
    dt = 1.0       # s
    S_energy = -1e5 # W/m^3
    S_specie = 1e-3 # kg/(m^3 s)
    Cp = 1000.0    # J/(kg K)

    # Analytic predictions
    # rho * Cp * dT/dt = S_energy  => dT = S_energy * dt / (rho * Cp)
    deltaT_analytic = S_energy * dt / (rho_init * Cp)
    T_analytic = T_init + deltaT_analytic

    # rho * dY/dt = S_specie  => dY = S_specie * dt / rho
    deltaY_analytic = S_specie * dt / rho_init
    Y_analytic = Y_init + deltaY_analytic

    deltaT_num = T_final - T_init
    deltaY_num = Y_final - Y_init

    print("\n--- Verification Results ---")
    print(f"Initial T:           {T_init:.4f} K")
    print(f"Final T (Numerical): {T_final:.4f} K (deltaT = {deltaT_num:.4f} K)")
    print(f"Final T (Analytic):  {T_analytic:.4f} K (deltaT = {deltaT_analytic:.4f} K)")
    print(f"T Rel Error:         {abs(T_final - T_analytic) / abs(T_analytic):.6e}")

    print(f"\nInitial YH2O:        {Y_init:.6f}")
    print(f"Final YH2O (Num):    {Y_final:.6f} (deltaY = {deltaY_num:.6e})")
    print(f"Final YH2O (Anal):   {Y_analytic:.6f} (deltaY = {deltaY_analytic:.6e})")
    print(f"YH2O Rel Error:      {abs(Y_final - Y_analytic) / abs(Y_analytic):.6e}")

    # Validation checks
    assert T_final < T_init, f"ERROR: Temperature rose ({T_final:.2f} K) instead of dropping! Sign is flipped!"
    assert Y_final > Y_init, f"ERROR: Vapor fraction dropped ({Y_final:.6f}) instead of rising! Sign is flipped!"
    
    t_err = abs(T_final - T_analytic) / abs(T_analytic)
    y_err = abs(Y_final - Y_analytic) / abs(Y_analytic)

    assert t_err < 1e-2, f"Temperature error too high: {t_err}"
    assert y_err < 1e-2, f"Species fraction error too high: {y_err}"

    print("\nSUCCESS: All plumbing sign checks and analytic verification passed!")

if __name__ == "__main__":
    main()
