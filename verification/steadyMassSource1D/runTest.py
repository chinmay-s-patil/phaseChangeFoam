#!/usr/bin/env python3
"""
Verification test for steady-state mass source sign and magnitude in phaseChangeMultiRegionSimpleFoam.
Solves a 1D steady channel with uniform positive mass source mdot and fixed-pressure outlet.
Verifies that:
1. Outflow velocity is positive: u_out > 0.
2. Outlet mass flow rate matches analytic prediction: M_dot_out = mdot * L * A within 1% tolerance.
"""

import os
import sys
import shutil
import subprocess

def find_solver(name="phaseChangeMultiRegionSimpleFoam"):
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        path = os.path.join(user_appbin, name)
        if os.path.exists(path):
            return path
    path = shutil.which(name)
    if path:
        return path
    home = os.environ.get("HOME", "/home/lavender")
    matches = [
        os.path.join(home, "OpenFOAM", f"lavender-v2412/platforms/linux64GccDPInt32Opt/bin/{name}")
    ]
    for m in matches:
        if os.path.exists(m):
            return m
    return name

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Command failed: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
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
        # vector or scalar
        val_str = val_str.replace("(", "").replace(")", "").strip()
        parts = [float(x) for x in val_str.split()]
        return parts
    
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    
    block = content[start_paren+1:end_paren].strip()
    return [float(x) for x in block.split()]

def setup_case(case_dir, mdot=0.1, L=1.0, W=0.1, H=0.1, n_cells=20, iters=50):
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir, exist_ok=True)

    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionSimpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {iters}; deltaT 1;
writeControl timeStep; writeInterval {iters}; purgeWrite 0; writeFormat ascii; writePrecision 8;
""")

    # system/blockMeshDict
    bm_content = f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }}
scale 1;
vertices
(
    (0 0 0)
    ({L} 0 0)
    ({L} {H} 0)
    (0 {H} 0)
    (0 0 {W})
    ({L} 0 {W})
    ({L} {H} {W})
    (0 {H} {W})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) air ({n_cells} 1 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    inlet
    {{
        type wall;
        faces ((0 4 7 3));
    }}
    outlet
    {{
        type patch;
        faces ((1 2 6 5));
    }}
    topAndBottom
    {{
        type empty;
        faces ((3 7 6 2) (0 1 5 4));
    }}
    frontAndBack
    {{
        type empty;
        faces ((0 3 2 1) (4 5 6 7));
    }}
);
"""
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(bm_content)
    os.makedirs(os.path.join(case_dir, "system/air"), exist_ok=True)
    with open(os.path.join(case_dir, "system/air/blockMeshDict"), "w") as f:
        f.write(bm_content.replace('"system"', '"system/air"'))

    schemes_content = """
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSchemes; }
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss upwind;
    div(phi,h) bounded Gauss upwind;
    div(phi,e) bounded Gauss upwind;
    div(phi,K) bounded Gauss upwind;
    div(phi,Yi_h) Gauss upwind;
    div(phi,H2O) Gauss upwind;
    div(phi,air) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write(schemes_content)

    solution_content = """
FoamFile { version 2.0; format ascii; class dictionary; location "system"; object fvSolution; }
solvers
{
    "h.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "Yi.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "air.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "H2O.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "p.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0.01; }
    "p_rgh.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0.01; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0.01; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 101325;
}
relaxationFactors
{
    fields
    {
        p 0.3;
        p_rgh 0.3;
        rho 0.5;
    }
    equations
    {
        U 0.7;
        h 0.7;
        H2O 0.7;
        air 0.7;
    }
}
"""
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(solution_content)

    # system/air/fvSchemes & fvSolution
    os.makedirs(os.path.join(case_dir, "system/air"), exist_ok=True)
    with open(os.path.join(case_dir, "system/air/fvSchemes"), "w") as f:
        f.write(schemes_content.replace('"system"', '"system/air"'))

    with open(os.path.join(case_dir, "system/air/fvSolution"), "w") as f:
        f.write(solution_content.replace('"system"', '"system/air"'))

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
regions ( fluid (air) solid () porousFluid () porousSolid () );
""")

    # constant/air
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
species ( air H2O );
inertSpecie air;
air { specie { molWeight 28.96; } thermodynamics { Cp 1000; Hf 0; } transport { mu 1.8e-5; Pr 0.7; } }
H2O { specie { molWeight 18.015; } thermodynamics { Cp 1000; Hf 0; } transport { mu 1.8e-5; Pr 0.7; } }
""")

    with open(os.path.join(air_const, "phaseChangeDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/air"; object phaseChangeDict; }}
phaseChange
{{
    active true;
    type constantSource;
    energySource 0;
    speciesName H2O;
    massSource {mdot};
    speciesSource {mdot};
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

    # 0/air
    air_zero = os.path.join(case_dir, "0/air")
    os.makedirs(air_zero, exist_ok=True)

    with open(os.path.join(air_zero, "p"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object p; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type calculated; value uniform 101325; }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(air_zero, "p_rgh"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101325;
boundaryField
{
    inlet { type fixedFluxPressure; value uniform 101325; }
    outlet { type fixedValue; value uniform 101325; }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(air_zero, "U"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volVectorField; location "0/air"; object U; }
dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (0 0 0); }
    outlet { type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(air_zero, "T"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object T; }
dimensions [0 0 0 1 0 0 0]; internalField uniform 300;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type zeroGradient; }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(air_zero, "air"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object air; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 1.0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type zeroGradient; }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(air_zero, "H2O"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class volScalarField; location "0/air"; object H2O; }
dimensions [0 0 0 0 0 0 0]; internalField uniform 0.0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type zeroGradient; }
    topAndBottom { type empty; }
    frontAndBack { type empty; }
}
""")

def run_test():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(test_dir, "case_steady1D")
    
    mdot = 0.5  # kg/(m^3 s)
    L = 1.0     # m
    W = 0.1     # m
    H = 0.1     # m
    A = W * H   # m^2 = 0.01 m^2
    n_cells = 20
    iters = 100

    print("=" * 80)
    print(" 1D Steady Mass Source Verification: phaseChangeMultiRegionSimpleFoam")
    print("=" * 80)
    print(f"Domain: L={L} m, A={A} m^2 (W={W} m, H={H} m), cells={n_cells}")
    print(f"Volumetric mass source: mdot = {mdot} kg/(m^3 s)")
    print(f"Expected mass outflow: M_dot_expected = mdot * L * A = {mdot * L * A:.6e} kg/s\n")

    setup_case(case_dir, mdot=mdot, L=L, W=W, H=H, n_cells=n_cells, iters=iters)

    # Mesh generation directly into constant/air/polyMesh
    print("--- Running blockMesh -region air ---")
    run_cmd("blockMesh -region air", cwd=case_dir)

    # Run solver
    solver = find_solver("phaseChangeMultiRegionSimpleFoam")
    print(f"--- Running {solver} ---")
    res = run_cmd(f"{solver}", cwd=case_dir)

    # Check outlet velocity and flux at final time
    time_dir = os.path.join(case_dir, f"{iters}/air")
    if not os.path.exists(time_dir):
        print(f"ERROR: Final time directory {time_dir} does not exist!")
        sys.exit(1)

    # Parse outlet patch from phi (exact face mass flux)
    with open(os.path.join(time_dir, "phi"), "r") as f:
        phi_content = f.read()

    outlet_idx = phi_content.find("outlet")
    assert outlet_idx != -1, "ERROR: outlet patch not found in phi field!"
    sub_phi = phi_content[outlet_idx:outlet_idx+200]
    val_idx = sub_phi.find("value")
    assert val_idx != -1, "ERROR: value not found in outlet patch of phi!"
    phi_val_str = sub_phi[val_idx:].split(";")[0].replace("value", "").replace("uniform", "").strip()
    phi_out = float(phi_val_str)

    # Parse outlet patch from U
    with open(os.path.join(time_dir, "U"), "r") as f:
        u_content = f.read()

    outlet_u_idx = u_content.find("outlet")
    assert outlet_u_idx != -1, "ERROR: outlet patch not found in U field!"
    sub_u = u_content[outlet_u_idx:outlet_u_idx+200]
    val_idx = sub_u.find("value")
    assert val_idx != -1, "ERROR: value not found in outlet patch of U!"
    u_val_str = sub_u[val_idx:].split(";")[0].replace("value", "").replace("uniform", "").replace("(", "").replace(")", "").strip()
    ux_out = float(u_val_str.split()[0])

    expected_flux = mdot * L * A
    rel_error = abs(phi_out - expected_flux) / expected_flux

    print("\n--- Quantitative Verification Results ---")
    print(f"Outlet velocity Ux(L):       {ux_out:.6f} m/s")
    print(f"Outlet mass outflow (phi):   {phi_out:.6e} kg/s")
    print(f"Expected mass outflow rate:  {expected_flux:.6e} kg/s")
    print(f"Relative error:              {rel_error * 100:.4f}%")

    # Assertions
    assert ux_out > 0, f"FAIL: Outlet velocity {ux_out} <= 0! Mass source drove flow backwards (wrong sign in pEqn)!"
    assert phi_out > 0, f"FAIL: Outflow {phi_out} <= 0!"
    assert rel_error < 0.01, f"FAIL: Outflow discrepancy {rel_error*100:.2f}% exceeds 1% tolerance!"

    # Clean up case dir to keep repository clean
    shutil.rmtree(case_dir, ignore_errors=True)

    print("\nSUCCESS: Uniform positive mass source drives positive outflow matching mdot*L*A within tolerance!")
    return True

if __name__ == "__main__":
    run_test()
