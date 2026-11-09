#!/usr/bin/env python3
"""
de Vahl Davis (1983) Natural Convection Benchmark (Ra = 10^5, Pr = 0.71)
Rung 1 of the Natural Convection & Phase Change Validation Ladder.

Validates:
  1. Steady-state buoyancy-driven natural convection with Boussinesq approximation
  2. p_rgh and velocity-pressure coupling in a closed square cavity
  3. Sensible enthalpy transport and wall heat transfer rate (Nusselt number)

Reference:
  G. de Vahl Davis (1983), "Natural convection of air in a square cavity:
  a bench mark numerical solution", Int. J. Num. Meth. Fluids, 3, pp. 249-264.
  Benchmark values at Ra = 10^5, Pr = 0.71:
    Mean Nusselt number:     Nu_mean = 4.519 (approx 4.52)
    Maximum Nusselt number:  Nu_max  = 7.717
    Minimum Nusselt number:  Nu_min  = 0.729
"""

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np

BENCHMARK_NU_MEAN = 4.519
BENCHMARK_NU_MAX = 7.717
BENCHMARK_NU_MIN = 0.729
BENCHMARK_U_MAX = 34.73
BENCHMARK_Y_UMAX = 0.855
BENCHMARK_V_MAX = 68.59
BENCHMARK_X_VMAX = 0.066
TOLERANCE_NU_REL = 0.05    # Mean Nu within 5% of benchmark
TOLERANCE_VMAX_REL = 0.05  # v_max* within 5% of benchmark
TOLERANCE_UMAX_REL = 0.05  # u_max* within 5% of benchmark

def find_solver(solver_name="phaseChangeMultiRegionSimpleFoam"):
    import glob
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        p = os.path.join(user_appbin, solver_name)
        if os.path.exists(p):
            return p
    which_p = shutil.which(solver_name)
    if which_p:
        return which_p
    home = os.path.expanduser("~")
    for cand in glob.glob(os.path.join(home, "OpenFOAM", "*-v*", "platforms", "*", "bin", solver_name)):
        if os.path.exists(cand):
            return cand
    return solver_name

def run_cmd(cmd, cwd=None, allow_failure=False):
    res = subprocess.run(
        cmd,
        shell=True,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if res.returncode != 0 and not allow_failure:
        print(f"Command failed: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res

def setup_case(case_dir, n_cells=50, iters=400):
    os.makedirs(case_dir, exist_ok=True)
    
    # 1. system files
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system/fluid"), exist_ok=True)
    
    L = 0.01  # 1 cm cavity width & height
    depth = L / 10.0
    
    block_mesh_dict = f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
scale 1;
vertices
(
    (0 0 0) ({L} 0 0) ({L} {L} 0) (0 {L} 0)
    (0 0 {depth}) ({L} 0 {depth}) ({L} {L} {depth}) (0 {L} {depth})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) fluid ({n_cells} {n_cells} 1) simpleGrading (1 1 1)
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
        f.write(block_mesh_dict)

    control_dict = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionSimpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {iters};
deltaT 1;
writeControl timeStep;
writeInterval {iters};
runTimeModifiable true;
"""
    with open(os.path.join(case_dir, "system/controlDict"), "w") as f:
        f.write(control_dict)

    fv_schemes = """FoamFile { version 2.0; format ascii; class dictionary; location "system/fluid"; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss limitedLinear 1;
    div(phi,h) Gauss limitedLinear 1;
    div(phi,K) Gauss limitedLinear 1;
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

    fv_solution = """FoamFile { version 2.0; format ascii; class dictionary; location "system/fluid"; object fvSolution; }
solvers
{
    "(p|p_rgh).*"
    {
        solver PCG;
        preconditioner DIC;
        tolerance 1e-7;
        relTol 0.01;
    }
    "(U|h).*"
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-7;
        relTol 0.01;
    }
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
        p_rgh 0.3;
        rho 1.0;
    }
    equations
    {
        U 0.5;
        h 0.7;
    }
}
"""
    with open(os.path.join(case_dir, "system/fluid/fvSolution"), "w") as f:
        f.write(fv_solution)
    with open(os.path.join(case_dir, "system/fvSolution"), "w") as f:
        f.write(fv_solution)

    # 2. constant files
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant/fluid"), exist_ok=True)

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

    with open(os.path.join(case_dir, "constant/fluid/phaseChangeDict"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }
active false;
phaseChange { active false; type none; }
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
    specie { molWeight 28.9; }
    equationOfState { rho0 1.0; T0 300.0; beta 0.003333333333; }
    thermodynamics { Cp 1000.0; Hf 0; }
    transport { mu 1.5237e-06; Pr 0.71; }
}
""")

    # 3. 0/fluid fields
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0/fluid"), exist_ok=True)

    with open(os.path.join(case_dir, "0/fluid/T"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object T; }
dimensions [0 0 0 1 0 0 0];
internalField uniform 300.0;
boundaryField
{
    hotWall { type fixedValue; value uniform 305.0; }
    coldWall { type fixedValue; value uniform 295.0; }
    topWall { type zeroGradient; }
    bottomWall { type zeroGradient; }
    frontAndBack { type empty; }
}
""")

    with open(os.path.join(case_dir, "0/fluid/h"), "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class volScalarField; location "0/fluid"; object h; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 300000.0;
boundaryField
{
    hotWall { type fixedValue; value uniform 305000.0; }
    coldWall { type fixedValue; value uniform 295000.0; }
    topWall { type zeroGradient; }
    bottomWall { type zeroGradient; }
    frontAndBack { type empty; }
}
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

def parse_temperature_field(filepath, n_cells):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"T field not found at {filepath}")
    with open(filepath, 'r') as f:
        content = f.read()
    lines = content.splitlines()
    idx = lines.index('(')
    total_cells = n_cells * n_cells
    t_vals = [float(lines[i]) for i in range(idx + 1, idx + 1 + total_cells)]
    return t_vals

def parse_velocity_field(filepath, n_cells):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"U field not found at {filepath}")
    with open(filepath, 'r') as f:
        content = f.read()
    idx = content.find('(')
    end = content.rfind(')')
    block = content[idx+1:end].strip()
    lines = [l.strip().strip('()') for l in block.split('\n') if l.strip()]
    u_list = []
    for l in lines:
        parts = l.split()
        if len(parts) == 3:
            u_list.append([float(parts[0]), float(parts[1]), float(parts[2])])
    return np.array(u_list)

def run_devahldavis(case_dir, n_cells=50, iters=400):
    print("=" * 70)
    print("  de Vahl Davis (1983) Natural Convection Benchmark (Ra = 10^5)")
    print("=" * 70)
    print(f"Setting up case in: {case_dir}")
    print(f"Mesh resolution:    {n_cells} x {n_cells}")
    print(f"SIMPLE iterations:  {iters}")

    setup_case(case_dir, n_cells=n_cells, iters=iters)

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"
    solver_bin = find_solver("phaseChangeMultiRegionSimpleFoam")

    # Clean old time steps and set up region mesh
    run_cmd(f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh'")
    print("Mesh generation and region setup completed.")

    # Run solver
    print(f"Running {solver_bin}...")
    run_cmd(f"bash -c '{of_env}; cd {case_dir} && {solver_bin} > log.simulation 2>&1'")
    print("Simulation completed.")

    # Parse T at endTime
    t_file = os.path.join(case_dir, str(iters), "fluid", "T")
    t_vals = parse_temperature_field(t_file, n_cells)

    # Physical bounds assert: T in [T_cold - 0.1, T_hot + 0.1]
    t_min = float(np.min(t_vals))
    t_max = float(np.max(t_vals))
    assert t_min >= 295.0 - 0.1 and t_max <= 305.0 + 0.1, f"FAIL: Temperature out of bounds: [{t_min}, {t_max}]"

    # In OpenFOAM blockMesh with (Nx Ny 1), cell ordering is:
    # cell = i + n_cells * j
    # For hot wall (x=0), adjacent cells are i=0:
    t_hot_cells = np.array([t_vals[n_cells * j] for j in range(n_cells)])
    
    L = 0.01
    dx = L / float(n_cells)
    xc = 0.5 * dx
    T_hot = 305.0
    T_cold = 295.0
    delta_T = T_hot - T_cold

    # Local Nusselt number: Nu(y) = (q_wall * L) / (k * delta_T) = (dT/dx)_wall * L / delta_T
    nu_local = ((T_hot - t_hot_cells) / xc) * (L / delta_T)
    nu_mean = float(np.mean(nu_local))
    nu_max = float(np.max(nu_local))
    nu_min = float(np.min(nu_local))

    rel_err = abs(nu_mean - BENCHMARK_NU_MEAN) / BENCHMARK_NU_MEAN

    # Parse U at endTime and evaluate non-dimensional midplane extrema
    u_file = os.path.join(case_dir, str(iters), "fluid", "U")
    U_arr = parse_velocity_field(u_file, n_cells)
    u_x = U_arr[:, 0].reshape((n_cells, n_cells))
    u_y = U_arr[:, 1].reshape((n_cells, n_cells))

    alpha_therm = 1.5237e-06 / 0.71
    u_scale = alpha_therm / L
    u_star = u_x / u_scale
    v_star = u_y / u_scale

    mid_col = n_cells // 2
    if n_cells % 2 == 0:
        mid_x_u = 0.5 * (u_star[:, mid_col - 1] + u_star[:, mid_col])
    else:
        mid_x_u = u_star[:, mid_col]
    y_coords = (np.arange(n_cells) + 0.5) / float(n_cells)
    idx_umax = int(np.argmax(np.abs(mid_x_u)))
    u_max_comp = float(np.abs(mid_x_u[idx_umax]))
    y_umax_comp = float(y_coords[idx_umax])

    mid_row = n_cells // 2
    if n_cells % 2 == 0:
        mid_y_v = 0.5 * (v_star[mid_row - 1, :] + v_star[mid_row, :])
    else:
        mid_y_v = v_star[mid_row, :]
    x_coords = (np.arange(n_cells) + 0.5) / float(n_cells)
    idx_vmax = int(np.argmax(np.abs(mid_y_v)))
    v_max_comp = float(np.abs(mid_y_v[idx_vmax]))
    x_vmax_comp = float(x_coords[idx_vmax])

    rel_err_vmax = abs(v_max_comp - BENCHMARK_V_MAX) / BENCHMARK_V_MAX
    rel_err_umax = abs(u_max_comp - BENCHMARK_U_MAX) / BENCHMARK_U_MAX

    print("\n" + "-" * 75)
    print("  RESULTS COMPARISON WITH DE VAHL DAVIS (1983) BENCHMARK (Ra = 10^5)")
    print("-" * 75)
    print(f"  {'Metric':<30} {'Computed':<15} {'Benchmark':<15} {'Rel. Difference':<15}")
    print("-" * 75)
    print(f"  {'Mean Nusselt (Nu_mean)':<30} {nu_mean:<15.4f} {BENCHMARK_NU_MEAN:<15.4f} {rel_err*100:<14.2f}%")
    print(f"  {'Max Nusselt (Nu_max)':<30} {nu_max:<15.4f} {BENCHMARK_NU_MAX:<15.4f} {(nu_max - BENCHMARK_NU_MAX)/BENCHMARK_NU_MAX*100:<14.2f}%")
    print(f"  {'Min Nusselt (Nu_min)':<30} {nu_min:<15.4f} {BENCHMARK_NU_MIN:<15.4f} {(nu_min - BENCHMARK_NU_MIN)/BENCHMARK_NU_MIN*100:<14.2f}%")
    print(f"  {'Max v* on y/L = 0.5':<30} {v_max_comp:<15.2f} {BENCHMARK_V_MAX:<15.2f} {rel_err_vmax*100:<14.2f}%")
    print(f"  {'Location x/L of v_max*':<30} {x_vmax_comp:<15.3f} {BENCHMARK_X_VMAX:<15.3f} -")
    print(f"  {'Max u* on x/L = 0.5':<30} {u_max_comp:<15.2f} {BENCHMARK_U_MAX:<15.2f} {rel_err_umax*100:<14.2f}%")
    print(f"  {'Location y/L of u_max*':<30} {y_umax_comp:<15.3f} {BENCHMARK_Y_UMAX:<15.3f} -")
    print("-" * 75)

    passed = (rel_err <= TOLERANCE_NU_REL) and (rel_err_vmax <= TOLERANCE_VMAX_REL) and (rel_err_umax <= TOLERANCE_UMAX_REL)
    status_str = "PASSED" if passed else "FAILED"
    print(f"Validation Status: {status_str} (Tolerances: Nu_mean <= {TOLERANCE_NU_REL*100:.1f}%, v_max* <= {TOLERANCE_VMAX_REL*100:.1f}%, u_max* <= {TOLERANCE_UMAX_REL*100:.1f}%, T bounds respected)\n")

    # Save results to text file
    res_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "devahldavis_validation_results.txt")
    with open(res_path, "w") as f:
        f.write("=" * 75 + "\n")
        f.write("  de Vahl Davis (1983) Natural Convection Benchmark (Ra = 10^5)\n")
        f.write("=" * 75 + "\n")
        f.write(f"Mesh resolution: {n_cells} x {n_cells}\n")
        f.write(f"SIMPLE iterations: {iters}\n\n")
        f.write(f"{'Metric':<30} {'Computed':<15} {'Benchmark':<15} {'Rel. Diff':<15}\n")
        f.write("-" * 75 + "\n")
        f.write(f"{'Mean Nusselt (Nu_mean)':<30} {nu_mean:<15.4f} {BENCHMARK_NU_MEAN:<15.4f} {rel_err*100:<14.2f}%\n")
        f.write(f"{'Max Nusselt (Nu_max)':<30} {nu_max:<15.4f} {BENCHMARK_NU_MAX:<15.4f} {(nu_max - BENCHMARK_NU_MAX)/BENCHMARK_NU_MAX*100:<14.2f}%\n")
        f.write(f"{'Min Nusselt (Nu_min)':<30} {nu_min:<15.4f} {BENCHMARK_NU_MIN:<15.4f} {(nu_min - BENCHMARK_NU_MIN)/BENCHMARK_NU_MIN*100:<14.2f}%\n")
        f.write(f"{'Max v* on y/L = 0.5':<30} {v_max_comp:<15.2f} {BENCHMARK_V_MAX:<15.2f} {rel_err_vmax*100:<14.2f}%\n")
        f.write(f"{'Location x/L of v_max*':<30} {x_vmax_comp:<15.3f} {BENCHMARK_X_VMAX:<15.3f} -\n")
        f.write(f"{'Max u* on x/L = 0.5':<30} {u_max_comp:<15.2f} {BENCHMARK_U_MAX:<15.2f} {rel_err_umax*100:<14.2f}%\n")
        f.write(f"{'Location y/L of u_max*':<30} {y_umax_comp:<15.3f} {BENCHMARK_Y_UMAX:<15.3f} -\n")
        f.write("-" * 75 + "\n")
        f.write(f"Validation Status: {status_str} (Tolerances: Nu_mean <= {TOLERANCE_NU_REL*100:.1f}%, v_max* <= {TOLERANCE_VMAX_REL*100:.1f}%, u_max* <= {TOLERANCE_UMAX_REL*100:.1f}%)\n")

    return passed

def run_grid_study(case_dir):
    print("=" * 85)
    print("  de Vahl Davis (1983) Mesh Convergence Study: 40x40 -> 50x50 -> 80x80 (Ra = 10^5)")
    print("=" * 85)
    configs = [
        (40, 500),
        (50, 500),
        (80, 1500),
    ]
    all_passed = True
    study_results = []
    for n, iters in configs:
        p = run_devahldavis(case_dir, n_cells=n, iters=iters)
        if not p:
            all_passed = False
    return all_passed

def main():
    parser = argparse.ArgumentParser(description="de Vahl Davis (1983) Natural Convection Benchmark")
    parser.add_argument("--mesh", type=int, default=50, help="Grid points in each direction (default: 50)")
    parser.add_argument("--iters", type=int, default=500, help="Number of SIMPLE iterations (default: 500)")
    parser.add_argument("--grid-study", action="store_true", help="Run full grid convergence study (40x40, 50x50, 80x80)")
    parser.add_argument("--clean", action="store_true", help="Remove simulation directory after running")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(this_dir, "case")

    if args.grid_study:
        passed = run_grid_study(case_dir)
    else:
        passed = run_devahldavis(case_dir, n_cells=args.mesh, iters=args.iters)

    if args.clean and os.path.exists(case_dir):
        shutil.rmtree(case_dir)

    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
