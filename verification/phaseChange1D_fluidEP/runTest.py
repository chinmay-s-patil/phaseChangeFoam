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
        if val_str.startswith("(") and val_str.endswith(")"):
            vec = [float(x) for x in val_str[1:-1].split()]
            return [vec] * num_cells
        return [float(val_str)] * num_cells

    idx_list = content.find("List<", idx)
    if idx_list == -1:
        return []
    start_paren = content.find("(", idx_list)
    if start_paren == -1:
        return []
    depth = 0
    end_paren = -1
    for i in range(start_paren, len(content)):
        if content[i] == '(':
            depth += 1
        elif content[i] == ')':
            depth -= 1
            if depth == 0:
                end_paren = i
                break
    if end_paren == -1:
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
            for p in line.split():
                vals.append(float(p))
    return vals

def setup_fluid_case(case_dir, frozen_flow=True, Cu=1e12, init_T=None, init_U=None, end_time=2000, delta_t=10, write_interval=200):
    if init_T and isinstance(init_T, list):
        num_cells = len(init_T)
    elif init_U and isinstance(init_U, list):
        num_cells = len(init_U)
    else:
        num_cells = 100

    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    if init_U:
        if isinstance(init_U, list):
            u_field = f"nonuniform List<vector> {num_cells} (\n" + "\n".join(f"({u[0]} {u[1]} {u[2]})" for u in init_U) + "\n);"
            u_bc_type = "zeroGradient"
            u_bc_val = ""
        else:
            u_vec_str = f"({init_U[0]} {init_U[1]} {init_U[2]})"
            u_field = f"uniform {u_vec_str};"
            u_bc_type = "zeroGradient"
            u_bc_val = ""
    else:
        u_field = "uniform (0 0 0);"
        u_bc_type = "fixedValue"
        u_bc_val = "value uniform (0 0 0);"
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end_time}; deltaT {delta_t};
writeControl runTime; writeInterval {write_interval}; purgeWrite 0; writeFormat ascii; writePrecision 12;
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
    "U.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-20; relTol 0; }}
}}
PIMPLE {{
    nOuterCorrectors 50;
    nCorrectors 2;
    nNonLinearCorrectors 50;
    nonLinearTolerance 1e-6;
    momentumPredictor yes;
    frozenFlow {frozen_str};
    pRefCell 0;
    pRefValue 101325;
}}
""")
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }}
scale 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) ({num_cells} 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot {{ type patch; faces ((0 4 7 3)); }}
    cold {{ type patch; faces ((1 2 6 5)); }}
    emptyFaces {{ type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
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
    if init_T:
        t_field = f"nonuniform List<scalar> {num_cells} (\n" + "\n".join(str(t) for t in init_T) + "\n);"
        t_bc = 'boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }'
        he_field = f"nonuniform List<scalar> {num_cells} (\n" + "\n".join(str((t - 298.15) * 1980.0) for t in init_T) + "\n);"
        he_bc = 'boundaryField { ".*" { type zeroGradient; } emptyFaces { type empty; } }'
        f_init = [(1.0 if t >= 310.0 else (0.0 if t <= 300.0 else (t-300.0)/10.0)) for t in init_T]
        pf_field = f"nonuniform List<scalar> {num_cells} (\n" + "\n".join(str(f) for f in f_init) + "\n);"
    else:
        t_field = "uniform 280.0;"
        t_bc = """boundaryField {
    hot { type fixedValue; value uniform 350.0; }
    cold { type zeroGradient; }
    emptyFaces { type empty; }
}"""
        he_field = "uniform -35937;"
        he_bc = """boundaryField {
    hot { type fixedValue; value uniform 102663; }
    cold { type zeroGradient; }
    emptyFaces { type empty; }
}"""
        pf_field = "uniform 0.0;"

    with open(os.path.join(case_dir, "0/pcm/phaseFraction"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object phaseFraction; }}
dimensions [0 0 0 0 0 0 0]; internalField {pf_field}
boundaryField {{ ".*" {{ type calculated; value uniform 0; }} emptyFaces {{ type empty; }} }}
""")

    with open(os.path.join(case_dir, "0/pcm/T"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object T; }}
dimensions [0 0 0 1 0 0 0]; internalField {t_field}
{t_bc}
""")
    with open(os.path.join(case_dir, "0/pcm/he"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volScalarField; location "0/pcm"; object he; }}
dimensions [0 2 -2 0 0 0 0]; internalField {he_field}
{he_bc}
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
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class volVectorField; location "0/pcm"; object U; }}
dimensions [0 1 -1 0 0 0 0]; internalField {u_field}
boundaryField {{ ".*" {{ type {u_bc_type}; {u_bc_val} }} emptyFaces {{ type empty; }} }}
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

def setup_solid_case(case_dir, Cu=1e12, init_T=None, end_time=2000, delta_t=10, write_interval=200):
    os.makedirs(case_dir, exist_ok=True)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    num_cells = len(init_T) if init_T else 100
    
    # system/controlDict
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end_time}; deltaT {delta_t};
writeControl runTime; writeInterval {write_interval}; purgeWrite 0; writeFormat ascii; writePrecision 12;
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
    "h.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    "he.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    "T.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
}
PIMPLE {
    nOuterCorrectors 50;
    nCorrectors 2;
    nNonLinearCorrectors 1;
    nonLinearTolerance 1e-6;
}
""")
    with open(os.path.join(case_dir, "system/blockMeshDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object blockMeshDict; }}
scale 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 0.01 0) (0 0.01 0) (0 0 0.01) (0.1 0 0.01) (0.1 0.01 0.01) (0 0.01 0.01) );
blocks ( hex (0 1 2 3 4 5 6 7) ({num_cells} 1 1) simpleGrading (1 1 1) );
edges ();
boundary (
    hot {{ type patch; faces ((0 4 7 3)); }}
    cold {{ type patch; faces ((1 2 6 5)); }}
    emptyFaces {{ type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
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
regions ( solid (pcm) fluid () porousFluid () porousSolid () );
""")
    
    os.makedirs(os.path.join(case_dir, "constant/pcm"), exist_ok=True)
    with open(os.path.join(case_dir, "constant/pcm/thermophysicalProperties"), "w") as f:
        f.write("""
FoamFile { version 2.0; format ascii; class dictionary; location "constant/pcm"; object thermophysicalProperties; }
thermoType { type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture {
    specie { molWeight 200.0; }
    transport { kappa 2.0; }
    thermodynamics { Cp 1980.0; Hf 0; }
    equationOfState { rho 1967.0; }
}
""")
    with open(os.path.join(case_dir, "constant/pcm/phaseChangeDict"), "w") as f:
        f.write(f"""
FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/pcm"; object phaseChangeDict; }}
active true;
phaseChange {{
    type enthalpyPorosity;
    active true;
    convection {{ suppress true; }}
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

def solve_stefan_neumann(kl=1.0, ks=2.0, rho=1967.0, cp=1980.0, L=163000.0, Th=350.0, Ts=300.0, Tl=310.0, T0=280.0, t=2000.0, num_cells=100, x_max=0.1):
    Tm = 0.5 * (Ts + Tl)
    al = kl / (rho * cp)
    as_ = ks / (rho * cp)

    stl = cp * (Th - Tm) / L
    sts = cp * (Tm - T0) / L

    def f(lam):
        term1 = stl / (math.exp(lam**2) * math.erf(lam))
        term2 = (sts * math.sqrt(as_/al)) / (math.exp(lam**2 * al/as_) * math.erfc(lam * math.sqrt(al/as_)))
        return term1 - term2 - lam * math.sqrt(math.pi)

    low, high = 0.001, 2.0
    for _ in range(100):
        mid = 0.5 * (low + high)
        if f(mid) > 0:
            low = mid
        else:
            high = mid
    lam = 0.5 * (low + high)
    X = 2 * lam * math.sqrt(al * t)

    dx = x_max / num_cells
    xs = [dx * 0.5 + i * dx for i in range(num_cells)]
    T_analytic = []
    for x in xs:
        if x <= X:
            tx = Th - (Th - Tm) * math.erf(x / (2 * math.sqrt(al * t))) / math.erf(lam)
        else:
            tx = T0 + (Tm - T0) * math.erfc(x / (2 * math.sqrt(as_ * t))) / math.erfc(lam * math.sqrt(al / as_))
        T_analytic.append(tx)

    return lam, X, T_analytic

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== Fluid Region Enthalpy-Porosity Verification Suite in {base_dir} ===")
    
    solver_bin = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"

    # Test 1: OpenFOAM momentumSp() Field Verification
    print("\n--- Test 1: OpenFOAM Solver momentumSp() Field Verification ---")
    c1_dir = os.path.join(base_dir, "case_drag_unit")
    shutil.rmtree(c1_dir, ignore_errors=True)
    setup_fluid_case(c1_dir, frozen_flow=False, Cu=1e5, init_T=[280.0, 305.0, 320.0], end_time=2, delta_t=2, write_interval=2)
    run_cmd(f"cd {c1_dir} && bash -c '{of_env}; {solver_bin}'")

    sp_vals = parse_openfoam_field(os.path.join(c1_dir, "2/pcm/momentumSp"), num_cells=3)
    
    # Formula: C_drag(f) = Cu * (1 - f)^2 / (f^3 + q)
    # Cell 0: T=280K (f=0.0) -> Expected: 1.00e+08
    # Cell 1: T=305K (f=0.5) -> Expected: 198412.70
    # Cell 2: T=320K (f=1.0) -> Expected: 0.00e+00
    c_f0 = sp_vals[0] if len(sp_vals) == 3 else -1.0
    c_f05 = sp_vals[1] if len(sp_vals) == 3 else -1.0
    c_f1 = sp_vals[2] if len(sp_vals) == 3 else -1.0

    err_f0 = abs(c_f0 - 1e8) / 1e8 if c_f0 >= 0 else 1.0
    err_f1 = abs(c_f1 - 0.0) if c_f1 >= 0 else 1.0
    err_f05 = abs(c_f05 - 198412.6984126984) / 198412.6984126984 if c_f05 >= 0 else 1.0

    pass_t1 = (len(sp_vals) == 3 and err_f0 < 1e-5 and err_f1 < 1e-5 and err_f05 < 1e-2)
    print(f"  OpenFOAM field 2/pcm/momentumSp at f=0.0: {c_f0:.2e} (Expected: 1.00e+08, rel err: {err_f0:.2e})")
    print(f"  OpenFOAM field 2/pcm/momentumSp at f=0.5: {c_f05:.2f} (Expected: 198412.70, rel err: {err_f05:.2e})")
    print(f"  OpenFOAM field 2/pcm/momentumSp at f=1.0: {c_f1:.2e} (Expected: 0.00e+00, rel err: {err_f1:.2e})")
    print(f"Test 1 OpenFOAM momentumSp() Field Verification: {'PASS' if pass_t1 else 'FAIL'}")

    # Test 2: Frozen-flow Stefan in Fluid Region vs Neumann Analytic Solution & Solid Benchmark
    print("\n--- Test 2: Frozen-flow Stefan vs Neumann Analytic Solution & Solid Benchmark (ks=2.0, kl=1.0) ---")
    c2_fluid_dir = os.path.join(base_dir, "case_frozen_flow")
    shutil.rmtree(c2_fluid_dir, ignore_errors=True)
    setup_fluid_case(c2_fluid_dir, frozen_flow=True, Cu=1e12)
    run_cmd(f"cd {c2_fluid_dir} && bash -c '{of_env}; {solver_bin}'")
    
    T_fluid = parse_openfoam_field(os.path.join(c2_fluid_dir, "2000/pcm/T"))
    a_fluid = parse_openfoam_field(os.path.join(c2_fluid_dir, "2000/pcm/phaseFraction"))

    # Solid-region reference run of identical problem
    c2_solid_dir = os.path.join(base_dir, "case_solid_ep_ref")
    shutil.rmtree(c2_solid_dir, ignore_errors=True)
    setup_solid_case(c2_solid_dir, Cu=1e12)
    run_cmd(f"cd {c2_solid_dir} && bash -c '{of_env}; {solver_bin}'")

    T_solid = parse_openfoam_field(os.path.join(c2_solid_dir, "2000/pcm/T"))
    a_solid = parse_openfoam_field(os.path.join(c2_solid_dir, "2000/pcm/phaseFraction"))

    # Neumann Analytical 2-Phase Stefan Solution
    lam_ana, X_ana, T_ana = solve_stefan_neumann(kl=1.0, ks=2.0, rho=1967.0, cp=1980.0, L=163000.0, Th=350.0, Ts=300.0, Tl=310.0, T0=280.0, t=2000.0)

    dT_fluid_solid = sum(abs(t1 - t2) for t1, t2 in zip(T_fluid, T_solid)) / len(T_fluid) if T_fluid and T_solid else 1e9
    da_fluid_solid = sum(abs(a1 - a2) for a1, a2 in zip(a_fluid, a_solid)) / len(a_fluid) if a_fluid and a_solid else 1e9
    dT_fluid_ana   = sum(abs(t1 - t2) for t1, t2 in zip(T_fluid, T_ana)) / len(T_fluid) if T_fluid and T_ana else 1e9

    # Find melt front position (x where alpha > 0.01)
    x_front_fluid = 0.0
    for idx, (alpha, cell_x) in enumerate(zip(a_fluid, [0.0005 + i*0.001 for i in range(100)])):
        if alpha > 0.01:
            x_front_fluid = cell_x

    x_front_solid = 0.0
    for idx, (alpha, cell_x) in enumerate(zip(a_solid, [0.0005 + i*0.001 for i in range(100)])):
        if alpha > 0.01:
            x_front_solid = cell_x

    dx_front_solid = abs(x_front_fluid - x_front_solid)
    dx_front_ana   = abs(x_front_fluid - X_ana)

    pass_t2 = (
        len(T_fluid) == 100 and len(T_solid) == 100 and len(T_ana) == 100
        and dT_fluid_solid < 0.01
        and da_fluid_solid < 0.001
        and dx_front_solid < 1e-4
        and dT_fluid_ana < 2.0
        and dx_front_ana < 0.007
    )

    print(f"  Neumann Analytical Melt Front X (t=2000s): {X_ana*1000:.2f} mm (Lambda = {lam_ana:.6f})")
    print(f"  Fluid Melt Front Position (x)             : {x_front_fluid*1000:.2f} mm (Analytic diff: {dx_front_ana*1000:.2f} mm)")
    print(f"  Fluid vs Solid Region Temp Mean Diff      : {dT_fluid_solid:.6f} K (Limit < 0.01 K)")
    print(f"  Fluid vs Neumann Analytic Temp Mean Diff  : {dT_fluid_ana:.4f} K (Limit < 2.00 K)")
    print(f"Test 2 Frozen-flow Stefan vs Neumann Analytic & Solid Physics Match: {'PASS' if pass_t2 else 'FAIL'}")

    # Test 3: Non-Zero Velocity Driver & Selective Darcy Damping (Cu = 1e12, U0 = 0.01 m/s)
    print("\n--- Test 3: Velocity Driver & Selective Darcy Drag Damping (Cu = 1e12, U0 = 0.01 m/s) ---")
    
    # 3a: Solid region case (T = 280K, f = 0.0) initialized with U = 0.01 m/s -> expect U to damp to ~0.0 m/s
    c3_solid_dir = os.path.join(base_dir, "case_drag_solid")
    shutil.rmtree(c3_solid_dir, ignore_errors=True)
    setup_fluid_case(c3_solid_dir, frozen_flow=False, Cu=1e12, init_T=[280.0]*100, init_U=(0.01, 0.0, 0.0), end_time=0.001, delta_t=0.001, write_interval=0.001)
    run_cmd(f"cd {c3_solid_dir} && bash -c '{of_env}; {solver_bin}'")
    U_solid = parse_openfoam_field(os.path.join(c3_solid_dir, "0.001/pcm/U"))
    max_U_solid = max(math.sqrt(v[0]**2 + v[1]**2 + v[2]**2) for v in U_solid) if U_solid else 1.0

    # 3b: Liquid region case (T = 320K, f = 1.0) initialized with U = 0.01 m/s -> expect U to remain 0.01 m/s
    c3_liquid_dir = os.path.join(base_dir, "case_drag_liquid")
    shutil.rmtree(c3_liquid_dir, ignore_errors=True)
    setup_fluid_case(c3_liquid_dir, frozen_flow=False, Cu=1e12, init_T=[320.0]*100, init_U=(0.01, 0.0, 0.0), end_time=0.001, delta_t=0.001, write_interval=0.001)
    run_cmd(f"cd {c3_liquid_dir} && bash -c '{of_env}; {solver_bin}'")
    U_liquid = parse_openfoam_field(os.path.join(c3_liquid_dir, "0.001/pcm/U"))
    avg_u_liquid_err = sum(abs(v[0] - 0.01) for v in U_liquid) / len(U_liquid) if U_liquid else 1.0

    pass_t3 = (len(U_solid) == 100 and len(U_liquid) == 100 and max_U_solid < 1e-8 and avg_u_liquid_err < 1e-4)
    print(f"  Initial Velocity U0                      : 0.0100 m/s (1.0 cm/s)")
    print(f"  Solid Region Velocity after 1 step       : {max_U_solid:.2e} m/s (Limit < 1e-8 m/s, DAMPED)")
    print(f"  Liquid Region Velocity Error after 1 step: {avg_u_liquid_err:.2e} m/s (Limit < 1e-4 m/s, UN-DAMPED)")
    print(f"Test 3 Velocity Driver Selective Darcy Damping: {'PASS' if pass_t3 else 'FAIL'}")

    print("\n=======================================================")
    print("    FLUID ENTHALPY-POROSITY VERIFICATION SUMMARY      ")
    print("=======================================================")
    print(f"Test 1 (Darcy Drag Coeff Unit Check)        : {'PASS' if pass_t1 else 'FAIL'}")
    print(f"Test 2 (Stefan vs Neumann & Solid Benchmark) : {'PASS' if pass_t2 else 'FAIL'}")
    print(f"Test 3 (Velocity Driver Selective Damping)  : {'PASS' if pass_t3 else 'FAIL'}")

    all_pass = pass_t1 and pass_t2 and pass_t3
    if all_pass:
        print("\nALL FLUID ENTHALPY-POROSITY VERIFICATION TESTS PASSED SUCCESSFULLY!")
    else:
        print("\nSTATUS: TEST FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
