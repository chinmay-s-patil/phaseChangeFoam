#!/usr/bin/env python3
"""
Run Brent et al. (1988) Gallium Melting Benchmark using EHCDarcy model.
Compares EHCDarcy directly against EnthalpyPorosity and Brent numerical benchmark.
"""

import os
import sys
import shutil
import subprocess
import numpy as np

REF_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference_data")
BENCHMARK_LIQUID_FRACTIONS = {
    120: 0.1102,   # 2 min
    360: 0.2346,   # 6 min
    600: 0.3465,   # 10 min
    1020: 0.5334,  # 17 min
}

def find_solver(solver_name="phaseChangeMultiRegionFoam"):
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        p = os.path.join(user_appbin, solver_name)
        if os.path.exists(p):
            return p
    which_p = shutil.which(solver_name)
    if which_p:
        return which_p
    home = os.environ.get("HOME", "/home/lavender")
    cand = os.path.join(home, "OpenFOAM", "lavender-v2412/platforms/linux64GccDPInt32Opt/bin", solver_name)
    if os.path.exists(cand):
        return cand
    return solver_name

def setup_ehc_case(case_dir, nx=42, ny=32, end_time=1020, deltaT_melt=0.5, cu=1.6e6):
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0/fluid"), exist_ok=True)

    W = 0.0889
    H = 0.0635
    D = 0.005

    # 1. blockMeshDict
    block_mesh = f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
scale 1;
vertices
(
    (0 0 0) ({W} 0 0) ({W} {H} 0) (0 {H} 0)
    (0 0 {D}) ({W} 0 {D}) ({W} {H} {D}) (0 {H} {D})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) fluid ({nx} {ny} 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    hotWall
    {{
        type wall;
        faces ((0 4 7 3));
    }}
    coldWall
    {{
        type wall;
        faces ((1 2 6 5));
    }}
    topWall
    {{
        type wall;
        faces ((3 7 6 2));
    }}
    bottomWall
    {{
        type wall;
        faces ((0 1 5 4));
    }}
    frontAndBack
    {{
        type empty;
        faces ((0 3 2 1) (4 5 6 7));
    }}
);
"""
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(block_mesh)

    # 2. controlDict
    control_dict = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT 0.1;
writeControl runTime;
writeInterval 60;
purgeWrite 0;
writeFormat ascii;
writePrecision 6;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable true;
adjustTimeStep false;
"""
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(control_dict)

    # 3. fvSchemes
    fv_schemes = """FoamFile { version 2.0; format ascii; class dictionary; location "system/fluid"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss upwind;
    div(phi,h) Gauss upwind;
    div(phi,K) Gauss upwind;
    div(phi,phaseFraction) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
    with open(os.path.join(case_dir, "system/fluid/fvSchemes"), "w") as f:
        f.write(fv_schemes)
    with open(os.path.join(case_dir, "system/fvSchemes"), "w") as f:
        f.write(fv_schemes)

    # 4. fvSolution
    fv_solution = """FoamFile { version 2.0; format ascii; class dictionary; location "system/fluid"; object fvSolution; }
solvers
{
    "(p|p_rgh).*"
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-8;
        relTol 0.001;
    }
    "(p|p_rgh)Final"
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-8;
        relTol 0;
    }
    "(U|h|rho).*"
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-8;
        relTol 0.001;
    }
    "(U|h|rho)Final"
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-8;
        relTol 0;
    }
    "rho.*"
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-8;
        relTol 0;
    }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 8;
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 101325;
}
relaxationFactors
{
    equations
    {
        U 0.7;
        h 0.7;
    }
}
"""
    with open(os.path.join(case_dir, "system/fluid/fvSolution"), "w") as f:
        f.write(fv_solution)
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(fv_solution)

    # 5. constant dictionaries
    with open(os.path.join(case_dir, "constant/regionProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object regionProperties; }
regions ( fluid (fluid) solid () );
""")

    with open(os.path.join(case_dir, "constant/g"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; location "constant"; object g; }
dimensions [0 1 -2 0 0 0 0];
value (0 -9.81 0);
""")

    with open(os.path.join(case_dir, "constant/fluid/turbulenceProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object turbulenceProperties; }
simulationType laminar;
""")

    with open(os.path.join(case_dir, "constant/fluid/radiationProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object radiationProperties; }
radiationModel none;
""")

    with open(os.path.join(case_dir, "constant/fluid/thermophysicalProperties"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object thermophysicalProperties; }
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie { molWeight 69.723; }
    equationOfState { rho0 6095.0; T0 302.93; beta 0.00012; }
    thermodynamics { Cp 381.5; Hf 0; Tref 0; Href 0; }
    transport { mu 0.00181; Pr 0.02157859375; }
}
""")

    T_m = 302.93
    T_l = T_m - deltaT_melt / 2.0
    T_u = T_m + deltaT_melt / 2.0

    # EHCDarcy phaseChangeDict
    with open(os.path.join(case_dir, "constant/fluid/phaseChangeDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
active true;
phaseChange
{{
    type EHCDarcy;
    active true;
    forward
    {{
        T_lowerBound {T_l:.4f};
        T_upperBound {T_u:.4f};
        latentHeat 80160.0;
    }}
    porosity
    {{
        Cu {cu:e};
        q 0.001;
    }}
    thermophysical
    {{
        mode custom;
        CpSolid 381.5;
        CpLiquid 381.5;
        kSolid 32.0;
        kLiquid 32.0;
    }}
    density
    {{
        rhoSolid 6095.0;
        rhoLiquid 6095.0;
    }}
}}
""")

    # 6. 0/fluid fields
    T_hot = 311.15
    T_cold = 301.45
    Cp = 381.5
    h_hot = Cp * T_hot
    h_cold = Cp * T_cold

    with open(os.path.join(case_dir, "0/fluid/T"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/fluid"; object T; }}
dimensions [0 0 0 1 0 0 0];
internalField uniform {T_cold};
boundaryField
{{
    hotWall {{ type fixedValue; value uniform {T_hot}; }}
    coldWall {{ type fixedValue; value uniform {T_cold}; }}
    topWall {{ type zeroGradient; }}
    bottomWall {{ type zeroGradient; }}
    frontAndBack {{ type empty; }}
}}
""")

    with open(os.path.join(case_dir, "0/fluid/h"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/fluid"; object h; }}
dimensions [0 2 -2 0 0 0 0];
internalField uniform {h_cold};
boundaryField
{{
    hotWall {{ type fixedValue; value uniform {h_hot}; }}
    coldWall {{ type fixedValue; value uniform {h_cold}; }}
    topWall {{ type zeroGradient; }}
    bottomWall {{ type zeroGradient; }}
    frontAndBack {{ type empty; }}
}}
""")

    with open(os.path.join(case_dir, "0/fluid/p"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object p; }
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 101325;
boundaryField
{
    hotWall { type zeroGradient; }
    coldWall { type zeroGradient; }
    topWall { type zeroGradient; }
    bottomWall { type zeroGradient; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/p_rgh"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 101325;
boundaryField
{
    hotWall { type fixedFluxPressure; value uniform 101325; }
    coldWall { type fixedFluxPressure; value uniform 101325; }
    topWall { type fixedFluxPressure; value uniform 101325; }
    bottomWall { type fixedFluxPressure; value uniform 101325; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/U"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volVectorField; location "0/fluid"; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    hotWall { type fixedValue; value uniform (0 0 0); }
    coldWall { type fixedValue; value uniform (0 0 0); }
    topWall { type fixedValue; value uniform (0 0 0); }
    bottomWall { type fixedValue; value uniform (0 0 0); }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/phaseFraction"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object phaseFraction; }
dimensions [0 0 0 0 0 0 0];
internalField uniform 0.0;
boundaryField
{
    hotWall { type calculated; value uniform 0.0; }
    coldWall { type calculated; value uniform 0.0; }
    topWall { type calculated; value uniform 0.0; }
    bottomWall { type calculated; value uniform 0.0; }
    frontAndBack { type empty; }
}
""")

def parse_scalar_field(filepath):
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return None
    sub = content[idx:idx+300]
    if "uniform" in sub and "nonuniform" not in sub:
        val_str = sub.split("uniform")[1].split(";")[0].strip()
        return np.array([float(val_str)])
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def evaluate_melt_front(alpha_2d, nx, ny, W, H):
    dx = W / float(nx)
    dy = H / float(ny)
    y_coords = (np.arange(ny) + 0.5) * dy
    x_front = np.zeros(ny)

    for j in range(ny):
        row = alpha_2d[j, :]
        idx = np.where(row >= 0.5)[0]
        if len(idx) == 0:
            if row[0] > 0.0:
                x_front[j] = row[0] * dx
            else:
                x_front[j] = 0.0
        else:
            i_last = idx[-1]
            if i_last + 1 < nx:
                a1 = row[i_last]
                a2 = row[i_last + 1]
                if abs(a2 - a1) > 1e-9:
                    frac = (0.5 - a1) / (a2 - a1)
                else:
                    frac = 0.0
                x_front[j] = (i_last + 0.5 + frac) * dx
            else:
                x_front[j] = (i_last + 0.5) * dx

    return y_coords, x_front

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(base_dir, "case_ehcDarcy")
    print(f"Setting up and running EHCDarcy case in: {case_dir}")

    setup_ehc_case(case_dir, nx=42, ny=32, end_time=1020, deltaT_melt=0.5, cu=1.6e6)

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver("phaseChangeMultiRegionFoam")

    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && {solver_bin} > log.simulation 2>&1'"
    print("Executing simulation...")
    res = subprocess.run(cmd, shell=True)
    print(f"Simulation completed with exit code: {res.returncode}")
