#!/usr/bin/env python3
"""
Brent, Voller, Reid (1988) Gallium Melting Benchmark
Rung 3 of the Natural Convection & Phase Change Validation Ladder.

Validates:
  1. Coupled thermal phase change (enthalpy-porosity formulation)
  2. Darcy drag damping in mushy and solid regions
  3. Natural convection circulation in the liquid melt pool (Ra ~ 6e5, Pr = 0.0216, Ste = 0.039)
  4. Evolution of melt front interface and total liquid fraction over time

Reference:
  A. D. Brent, V. R. Voller, and K. J. Reid (1988),
  "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change:
   Application to the Melting of a Pure Metal",
  Numerical Heat Transfer, 13(3), pp. 297-318.
  Paper copy: validation/papers/brent1988.pdf
  Digitized reference data: validation/gallium2D/reference_data/
"""

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np

# Benchmark reference values
REF_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_data")
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

def setup_case(case_dir, nx=42, ny=32, end_time=120, deltaT_melt=0.5, cu=1.6e6):
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant/fluid"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "0/fluid"), exist_ok=True)

    W = 0.0889  # 8.89 cm width
    H = 0.0635  # 6.35 cm height
    D = 0.005   # 2D thickness

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
    # Write at check times: 120s, 360s, 600s, 1020s (or every 60s)
    write_interval = 60 if end_time > 120 else 60
    control_dict = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT 0.1;
writeControl runTime;
writeInterval {write_interval};
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

    with open(os.path.join(case_dir, "constant/fluid/phaseChangeDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
active true;
phaseChange
{{
    type enthalpyPorosity;
    active true;
    forward
    {{
        T_lowerBound {T_l:.4f};
        T_upperBound {T_u:.4f};
        latentHeat 80160.0;
    }}
    lambda 0.2;
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
    # alpha_2d has shape (ny, nx)
    # For each row j (y-coordinate), find x where alpha = 0.5
    dx = W / float(nx)
    dy = H / float(ny)
    y_coords = (np.arange(ny) + 0.5) * dy
    x_front = np.zeros(ny)

    for j in range(ny):
        row = alpha_2d[j, :]
        # find where row crosses 0.5
        idx = np.where(row >= 0.5)[0]
        if len(idx) == 0:
            # solid everywhere in this row, front near wall
            if row[0] > 0.0:
                x_front[j] = row[0] * dx
            else:
                x_front[j] = 0.0
        else:
            i_last = idx[-1]
            if i_last + 1 < nx:
                # linear interpolation between i_last and i_last + 1
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

def run_gallium(case_dir, nx=42, ny=32, end_time=120, deltaT_melt=0.5, cu=1.6e6, plot=False):
    print("=" * 75)
    print("  Brent, Voller, Reid (1988) Gallium Melting Benchmark")
    print("=" * 75)
    print(f"Setting up case in:   {case_dir}")
    print(f"Mesh resolution:      {nx} x {ny}")
    print(f"End time:             {end_time} s ({end_time/60.0:.1f} min)")
    print(f"Melting window:       {deltaT_melt} K (302.93 +/- {deltaT_melt/2.0} K)")
    print(f"Darcy constant Cu:    {cu:e}")

    setup_case(case_dir, nx=nx, ny=ny, end_time=end_time, deltaT_melt=deltaT_melt, cu=cu)

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver("phaseChangeMultiRegionFoam")

    # Clean and mesh
    print("Generating mesh...")
    run_cmd(f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh'")
    print("Mesh generation completed.")

    # Run simulation
    print(f"Running {solver_bin}...")
    run_cmd(f"bash -c '{of_env}; cd {case_dir} && {solver_bin} > log.simulation 2>&1'")
    print("Simulation completed.")

    # Discover written time directories
    time_dirs = []
    for d in os.listdir(case_dir):
        try:
            t_val = float(d)
            if t_val > 0:
                time_dirs.append((t_val, d))
        except ValueError:
            pass
    time_dirs.sort()

    W = 0.0889
    H = 0.0635

    print("\n" + "-" * 75)
    print("  LIQUID FRACTION EVOLUTION vs BRENT ET AL. (1988)")
    print("-" * 75)
    print(f"  {'Time (s)':<12} {'Time (min)':<12} {'Computed f_L':<15} {'Benchmark f_L':<15} {'Rel. Difference':<15}")
    print("-" * 75)

    extracted_results = []
    overall_pass = True

    for t_val, t_str in time_dirs:
        alpha_file = os.path.join(case_dir, t_str, "fluid", "phaseFraction")
        if not os.path.exists(alpha_file):
            continue
        alpha = parse_scalar_field(alpha_file)
        if alpha is None or len(alpha) != nx * ny:
            continue

        f_L_computed = float(np.mean(alpha))
        nearest_ref_time = None
        for ref_t in sorted(BENCHMARK_LIQUID_FRACTIONS.keys()):
            if abs(t_val - ref_t) <= 5.0:
                nearest_ref_time = ref_t
                break

        if nearest_ref_time is not None:
            f_L_ref = BENCHMARK_LIQUID_FRACTIONS[nearest_ref_time]
            abs_err = abs(f_L_computed - f_L_ref)
            rel_err = abs_err / f_L_ref
            # Benchmark acceptance: <= 25% relative difference or <= 0.03 absolute liquid fraction difference
            passed_time = (rel_err <= 0.25) or (abs_err <= 0.03)
            if not passed_time:
                overall_pass = False
            diff_str = f"{rel_err*100:.2f}% (abs: {abs_err:.4f})"
            ref_str = f"{f_L_ref:.4f}"
        else:
            diff_str = "N/A"
            ref_str = "N/A"
            passed_time = True

        print(f"  {t_val:<12.1f} {t_val/60.0:<12.2f} {f_L_computed:<15.4f} {ref_str:<15} {diff_str:<15}")
        extracted_results.append((t_val, f_L_computed, ref_str, diff_str))

    print("-" * 75)

    # Status
    status_str = "PASSED" if overall_pass else "FAILED"
    print(f"Validation Status: {status_str} (Benchmark tolerance: <= 25% rel or <= 0.03 abs f_L)\n")

    # Evaluate melt fronts at available benchmark times (2, 6, 10, 17 min)
    benchmark_front_times = [
        (120, 2),
        (360, 6),
        (600, 10),
        (1020, 17),
    ]

    front_results = {}
    for t_sec, t_min in benchmark_front_times:
        alpha_path = os.path.join(case_dir, str(int(t_sec)), "fluid", "phaseFraction")
        if os.path.exists(alpha_path):
            alpha_fld = parse_scalar_field(alpha_path)
            if alpha_fld is not None and len(alpha_fld) == nx * ny:
                a2d = alpha_fld.reshape((ny, nx))
                y_pts, x_pts = evaluate_melt_front(a2d, nx, ny, W, H)

                ref_front_csv = os.path.join(REF_DATA_DIR, f"brent1988_front_{t_min}min.csv")
                mean_diff, max_diff = None, None
                if os.path.exists(ref_front_csv):
                    import csv
                    b_y, b_x = [], []
                    with open(ref_front_csv, 'r') as f:
                        r = csv.DictReader(f)
                        for row in r:
                            b_y.append(float(row["y_m"]))
                            b_x.append(float(row["x_m"]))
                    b_x_interp = np.interp(y_pts, b_y, b_x)
                    dx_front = np.abs(x_pts - b_x_interp)
                    mean_diff = float(np.mean(dx_front) * 1000.0)
                    max_diff = float(np.max(dx_front) * 1000.0)
                    print("-" * 75)
                    print(f"  MELT FRONT PROFILE AT {t_min} MIN ({t_sec} s) vs BRENT ET AL. (1988)")
                    print("-" * 75)
                    print(f"  Mean front position difference: {mean_diff:.2f} mm (cavity width = {W*1000:.1f} mm)")
                    print(f"  Max front position difference:  {max_diff:.2f} mm")
                    print(f"  Grid cell width dx:             {W/nx*1000:.2f} mm")
                    print(f"  Front position accuracy:        Within {max_diff / (W/nx*1000.0):.1f} grid cell widths")
                    print("-" * 75 + "\n")

                front_results[t_sec] = {
                    "min": t_min,
                    "y": y_pts,
                    "x": x_pts,
                    "mean_diff": mean_diff,
                    "max_diff": max_diff,
                }

    # Write results summary
    res_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gallium_validation_results.txt")
    with open(res_path, "w") as f:
        f.write("=" * 75 + "\n")
        f.write("  Brent, Voller, Reid (1988) Gallium Melting Benchmark Validation\n")
        f.write("=" * 75 + "\n")
        f.write(f"Cavity: {W*100:.2f} cm x {H*100:.2f} cm\n")
        f.write(f"Mesh: {nx} x {ny}\n")
        f.write(f"Simulation end time: {end_time} s ({end_time/60.0:.1f} min)\n")
        f.write(f"Phase change model: enthalpyPorosity (Cu={cu:e}, deltaT={deltaT_melt} K)\n\n")
        f.write(f"{'Time (s)':<12} {'Time (min)':<12} {'Computed f_L':<15} {'Benchmark f_L':<15} {'Rel. Diff':<15}\n")
        f.write("-" * 75 + "\n")
        for t_val, f_L_computed, ref_str, diff_str in extracted_results:
            f.write(f"{t_val:<12.1f} {t_val/60.0:<12.2f} {f_L_computed:<15.4f} {ref_str:<15} {diff_str:<15}\n")
        f.write("-" * 75 + "\n")
        f.write(f"Validation Status: {status_str}\n\n")
        if front_results:
            f.write("Melt Front Profile Deviations vs Brent et al. (1988):\n")
            for t_sec in sorted(front_results.keys()):
                fd = front_results[t_sec]
                if fd["mean_diff"] is not None:
                    f.write(f"  t = {fd['min']} min ({t_sec} s): Mean diff = {fd['mean_diff']:.2f} mm, Max diff = {fd['max_diff']:.2f} mm (within {fd['max_diff']/(W/nx*1000.0):.1f} cells)\n")

    # Optional plotting
    if plot:
        try:
            # If both Enthalpy-Porosity and EHCDarcy cases exist, use unified comparative plotter
            comp_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis", "plot_comparison_both_models.py")
            ehc_case = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis", "case_ehcDarcy")
            if os.path.exists(comp_script) and os.path.exists(ehc_case):
                subprocess.run(f"python3 {comp_script}", shell=True)
                return overall_pass

            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

            # Panel 1: Liquid fraction curve
            times_ref = sorted(BENCHMARK_LIQUID_FRACTIONS.keys())
            fL_ref_vals = [BENCHMARK_LIQUID_FRACTIONS[t] for t in times_ref]
            ax1.plot([t/60.0 for t in times_ref], fL_ref_vals, 'rs--', markersize=7, lw=1.8, label='Brent et al. (1988) Numerical')
            times_comp = [r[0]/60.0 for r in extracted_results]
            fL_comp = [r[1] for r in extracted_results]
            ax1.plot(times_comp, fL_comp, 'b-o', markersize=4, lw=2.0, label='Computed (phaseChangeFoam)')
            ax1.set_xlabel('Time (min)', fontsize=11)
            ax1.set_ylabel('Liquid Fraction $f_L$', fontsize=11)
            ax1.set_title('Gallium Melting: Liquid Fraction vs Time', fontsize=12, fontweight='bold')
            ax1.set_xlim(0, max(max(times_comp) if times_comp else 2, 17) * 1.05)
            ax1.set_ylim(0, 0.65)
            ax1.grid(True, linestyle=':', alpha=0.6)
            ax1.legend(fontsize=10)

            # Panel 2: Melt front comparison at available times
            color_map = {
                120: '#1f77b4',   # 2 min: Blue
                360: '#2ca02c',   # 6 min: Green
                600: '#ff7f0e',   # 10 min: Orange
                1020: '#d62728',  # 17 min: Red
            }

            if front_results:
                # Plot reference benchmark fronts first (dashed)
                for t_sec, t_min in benchmark_front_times:
                    if t_sec in front_results:
                        ref_front_csv = os.path.join(REF_DATA_DIR, f"brent1988_front_{t_min}min.csv")
                        if os.path.exists(ref_front_csv):
                            import csv
                            b_y, b_x = [], []
                            with open(ref_front_csv, 'r') as f:
                                r = csv.DictReader(f)
                                for row in r:
                                    b_y.append(float(row["y_m"]))
                                    b_x.append(float(row["x_m"]))
                            c = color_map.get(t_sec, 'gray')
                            ax2.plot(np.array(b_x) * 1000.0, np.array(b_y) * 1000.0,
                                     color=c, linestyle='--', lw=1.8,
                                     label=f'Brent {t_min} min')

                # Plot computed fronts (solid)
                for t_sec in sorted(front_results.keys()):
                    fd = front_results[t_sec]
                    t_min = fd["min"]
                    c = color_map.get(t_sec, 'blue')
                    ax2.plot(fd["x"] * 1000.0, fd["y"] * 1000.0,
                             color=c, linestyle='-', lw=2.4,
                             label=f'Computed {t_min} min')

                ax2.set_xlim(0, W * 1000.0)
                ax2.set_ylim(0, H * 1000.0)
                ax2.set_xlabel('Width x (mm) [Hot Wall -> Cold Wall]', fontsize=11)
                ax2.set_ylabel('Height y (mm)', fontsize=11)
                ax2.set_title('Melt Front Profiles vs Brent et al. (1988)', fontsize=12, fontweight='bold')
                ax2.axvline(x=0, color='red', lw=1.5, alpha=0.5)
                ax2.axvline(x=W*1000.0, color='blue', lw=1.5, alpha=0.5)
                ax2.grid(True, linestyle=':', alpha=0.6)
                ax2.legend(loc='lower right', fontsize=8, ncol=2)
            else:
                ax2.axis('off')

            plt.tight_layout()
            plot_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gallium_melt_fraction_comparison.png")
            fig.savefig(plot_path, dpi=150)
            plt.close(fig)
            print(f"Saved validation comparison plot to {plot_path}")
        except Exception as e:
            print(f"Note: Plotting skipped ({e})")

    return overall_pass

def main():
    parser = argparse.ArgumentParser(description="Brent et al. (1988) Gallium Melting Benchmark")
    parser.add_argument("--nx", type=int, default=42, help="Grid points in x direction (default: 42)")
    parser.add_argument("--ny", type=int, default=32, help="Grid points in y direction (default: 32)")
    parser.add_argument("--endTime", type=int, default=120, help="Simulation duration in seconds (default: 120 s = 2 min)")
    parser.add_argument("--deltaT", type=float, default=0.5, help="Melting interval width in K (default: 0.5)")
    parser.add_argument("--cu", type=float, default=1.6e6, help="Darcy mushy constant (default: 1.6e6)")
    parser.add_argument("--plot", action="store_true", default=True, help="Generate comparison plot")
    parser.add_argument("--clean", action="store_true", help="Remove simulation folder after running")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.join(this_dir, "case")

    passed = run_gallium(
        case_dir=case_dir,
        nx=args.nx,
        ny=args.ny,
        end_time=args.endTime,
        deltaT_melt=args.deltaT,
        cu=args.cu,
        plot=args.plot
    )

    if args.clean and os.path.exists(case_dir):
        shutil.rmtree(case_dir)

    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
