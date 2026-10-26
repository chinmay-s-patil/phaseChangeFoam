#!/usr/bin/env python3
"""
Master runner for comprehensive numerical verification:
1. Parameter sweep across PIMPLE outer correctors (nOuter in [2, 3, 4, 5, 6, 7, 8, 10, 12, 16])
   for EP (lambda = 1.0, 0.5, 0.2) and EHCDarcy (nOuter in [2, 3, 4, 5, 8]) to t = 120s.
2. 17-minute full-benchmark runs with per-timestep wallHeatFlux logging.
3. Performs exact energy audits and front-deviation analyses.
4. Generates committed JSON and formatted TXT data files for report generation.
"""

import os
import sys
import time
import shutil
import subprocess
import json
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_energy_balance_audit import run_exact_audit, parse_openfoam_field

OF_ENV = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
BASE_CASE = os.path.abspath("validation/gallium2D/case")
SWEEP_DIR = os.path.abspath("validation/gallium2D/analysis/energy_balance/sweep_runs")
STUDY_17MIN_DIR = os.path.abspath("validation/gallium2D/analysis/study_17min")
REF_DIR = os.path.abspath("validation/gallium2D/reference_data")

W = 0.0889
H = 0.0635
D = 0.005

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing: {cmd}\nSTDERR: {res.stderr}\nSTDOUT: {res.stdout}")
        raise RuntimeError(f"Command failed: {cmd}")
    return res

def setup_case(base_dir, case_dir, model_type="enthalpyPorosity", lam=1.0, n_outer=8, end_time=120):
    os.makedirs(case_dir, exist_ok=True)
    for folder in ["0", "constant", "system"]:
        src = os.path.join(base_dir, folder)
        dst = os.path.join(case_dir, folder)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

    # 1. Modify controlDict
    ctrl_path = os.path.join(case_dir, "system", "controlDict")
    with open(ctrl_path, "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
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

functions
{{
    wallHeatFlux1
    {{
        type            wallHeatFlux;
        libs            (fieldFunctionObjects);
        region          fluid;
        writeControl    timeStep;
        writeInterval   1;
        patches         (hotWall coldWall);
    }}
}}
""")

    # 2. Modify fvSolution
    import re
    sol_paths = [os.path.join(case_dir, "system", "fluid", "fvSolution"), os.path.join(case_dir, "system", "fvSolution")]
    for sp in sol_paths:
        if os.path.exists(sp):
            with open(sp, "r") as f:
                content = f.read()
            content = re.sub(r"nOuterCorrectors\s+\d+;", f"nOuterCorrectors {n_outer};", content)
            with open(sp, "w") as f:
                f.write(content)

    # 3. Modify phaseChangeDict
    pc_path = os.path.join(case_dir, "constant", "fluid", "phaseChangeDict")
    if model_type == "EHCDarcy":
        with open(pc_path, "w") as f:
            f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }
active true;
phaseChange
{
    type EHCDarcy;
    active true;
    forward
    {
        T_lowerBound 302.68;
        T_upperBound 303.18;
        latentHeat 80160.0;
    }
    porosity
    {
        Cu 1.6e6;
        q 0.001;
    }
    thermophysical
    {
        mode custom;
        CpSolid 381.5;
        CpLiquid 381.5;
        kSolid 32.0;
        kLiquid 32.0;
    }
    density
    {
        rhoSolid 6095.0;
        rhoLiquid 6095.0;
    }
}
""")
    else:
        with open(pc_path, "w") as f:
            f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
active true;
phaseChange
{{
    type enthalpyPorosity;
    active true;
    forward
    {{
        T_lowerBound 302.68;
        T_upperBound 303.18;
        latentHeat 80160.0;
    }}
    lambda {lam};
    porosity
    {{
        Cu 1.6e6;
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

def launch_simulation(case_dir):
    cmd = f"bash -c '{OF_ENV}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    p = subprocess.Popen(cmd, shell=True)
    return p

def evaluate_melt_front(alpha_2d, nx=42, ny=32):
    dx = W / float(nx)
    dy = H / float(ny)
    y_coords = (np.arange(ny) + 0.5) * dy
    x_front = np.zeros(ny)

    for j in range(ny):
        row = alpha_2d[j, :]
        idx = np.where(row >= 0.5)[0]
        if len(idx) == 0:
            x_front[j] = row[0] * dx if row[0] > 0.0 else 0.0
        else:
            i_last = idx[-1]
            if i_last + 1 < nx:
                a1 = row[i_last]
                a2 = row[i_last + 1]
                frac = (0.5 - a1) / (a2 - a1) if abs(a2 - a1) > 1e-9 else 0.0
                x_front[j] = (i_last + 0.5 + frac) * dx
            else:
                x_front[j] = (i_last + 0.5) * dx
    return y_coords, x_front

def calculate_front_deviation(case_dir, t_sec, nx=42, ny=32):
    alpha_path = os.path.join(case_dir, str(t_sec), "fluid", "phaseFraction")
    alpha_fld = parse_openfoam_field(alpha_path)
    if alpha_fld is None or len(alpha_fld) != nx * ny:
        return None, None
    a2d = alpha_fld.reshape((ny, nx))
    y_pts, x_pts = evaluate_melt_front(a2d, nx, ny)

    t_min = int(round(t_sec / 60.0))
    ref_front_csv = os.path.join(REF_DIR, f"brent1988_front_{t_min}min.csv")
    if not os.path.exists(ref_front_csv):
        return None, None

    b_y, b_x = [], []
    with open(ref_front_csv, 'r') as f:
        r = csv.DictReader(f)
        for row in r:
            b_y.append(float(row["y_m"]))
            b_x.append(float(row["x_m"]))
    b_x_interp = np.interp(y_pts, b_y, b_x)
    dx_front = np.abs(x_pts - b_x_interp)
    mean_diff_mm = float(np.mean(dx_front) * 1000.0)
    max_diff_mm = float(np.max(dx_front) * 1000.0)
    return mean_diff_mm, max_diff_mm

def parse_execution_time(log_path):
    if not os.path.exists(log_path):
        return None
    with open(log_path, 'r') as f:
        lines = f.readlines()
        for line in reversed(lines[-50:]):
            if "ExecutionTime =" in line:
                try:
                    parts = line.split("ExecutionTime =")[1].split("s")[0].strip()
                    return float(parts)
                except:
                    pass
    return None

def main():
    print("=" * 80)
    print(" STARTING MASTER COMPREHENSIVE NUMERICAL BENCHMARK SUITE")
    print("=" * 80)

    # 1. Check/Run 17-minute EP tutorial baseline (lambda = 0.2, nOuter = 8) if not completed
    ep_n8_17m_dir = os.path.join(STUDY_17MIN_DIR, "audit_17min_ep_n8")
    if not os.path.exists(os.path.join(ep_n8_17m_dir, "1020", "fluid", "phaseFraction")):
        print(f"\nSetting up 17-minute EP (lambda=0.2, nOuter=8) in {ep_n8_17m_dir} ...")
        setup_case(BASE_CASE, ep_n8_17m_dir, model_type="enthalpyPorosity", lam=0.2, n_outer=8, end_time=1020)
        p_17m = launch_simulation(ep_n8_17m_dir)
        p_17m.wait()
        print("  --> 17-minute EP run completed.")

    # 2. Define all 2-minute sweep cases
    sweep_cases = [
        # Model, lambda, nOuter, case_id
        # EHCDarcy sweep
        ("EHCDarcy", None, 2, "ehc_n2"),
        ("EHCDarcy", None, 3, "ehc_n3"),
        ("EHCDarcy", None, 4, "ehc_n4"),
        ("EHCDarcy", None, 5, "ehc_n5"),
        ("EHCDarcy", None, 8, "ehc_n8"),
        # EP lambda = 1.0 sweep
        ("enthalpyPorosity", 1.0, 2, "ep_lam1.0_n2"),
        ("enthalpyPorosity", 1.0, 3, "ep_lam1.0_n3"),
        ("enthalpyPorosity", 1.0, 4, "ep_lam1.0_n4"),
        ("enthalpyPorosity", 1.0, 5, "ep_lam1.0_n5"),
        ("enthalpyPorosity", 1.0, 6, "ep_lam1.0_n6"),
        ("enthalpyPorosity", 1.0, 7, "ep_lam1.0_n7"),
        ("enthalpyPorosity", 1.0, 8, "ep_lam1.0_n8"),
        ("enthalpyPorosity", 1.0, 10, "ep_lam1.0_n10"),
        ("enthalpyPorosity", 1.0, 12, "ep_lam1.0_n12"),
        ("enthalpyPorosity", 1.0, 16, "ep_lam1.0_n16"),
        # EP lambda = 0.5 sweep
        ("enthalpyPorosity", 0.5, 2, "ep_lam0.5_n2"),
        ("enthalpyPorosity", 0.5, 3, "ep_lam0.5_n3"),
        ("enthalpyPorosity", 0.5, 5, "ep_lam0.5_n5"),
        ("enthalpyPorosity", 0.5, 8, "ep_lam0.5_n8"),
        # EP lambda = 0.2 sweep
        ("enthalpyPorosity", 0.2, 2, "ep_lam0.2_n2"),
        ("enthalpyPorosity", 0.2, 3, "ep_lam0.2_n3"),
        ("enthalpyPorosity", 0.2, 5, "ep_lam0.2_n5"),
        ("enthalpyPorosity", 0.2, 8, "ep_lam0.2_n8"),
        ("enthalpyPorosity", 0.2, 12, "ep_lam0.2_n12"),
    ]

    os.makedirs(SWEEP_DIR, exist_ok=True)
    sweep_results = []

    print(f"\nProcessing {len(sweep_cases)} sweep cases ...")
    for model_type, lam, n_outer, cid in sweep_cases:
        cdir = os.path.join(SWEEP_DIR, cid)
        whf_file = os.path.join(cdir, "postProcessing", "fluid", "wallHeatFlux1", "0", "wallHeatFlux.dat")
        log_file = os.path.join(cdir, "log.simulation")

        if not (os.path.exists(whf_file) and os.path.exists(log_file)):
            print(f"Running simulation for {cid} ...")
            setup_case(BASE_CASE, cdir, model_type=model_type, lam=lam if lam else 1.0, n_outer=n_outer, end_time=120)
            p = launch_simulation(cdir)
            p.wait()

        # Extract OpenFOAM internal CPU execution time
        cpu_time = parse_execution_time(log_file)
        if cpu_time is None:
            cpu_time = 0.0

        # Run exact audit (right-endpoint rule)
        audit_res = run_exact_audit(cdir)
        row_120 = None
        for r in audit_res:
            if abs(r["t_sec"] - 120.0) < 1e-3:
                row_120 = r
                break
        if row_120:
            row_120["model"] = model_type
            row_120["lambda"] = lam
            row_120["n_outer"] = n_outer
            row_120["case_id"] = cid
            row_120["clock_time"] = cpu_time
            row_120["cpu_time"] = cpu_time
            row_120["total_outer_passes"] = int(1200 * n_outer)
            sweep_results.append(row_120)

    # 3. Save Sweep Results
    sweep_json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ep_sweep_results.json")
    with open(sweep_json_path, "w") as f:
        json.dump(sweep_results, f, indent=2)

    sweep_txt_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ep_sweep_results.txt")
    with open(sweep_txt_path, "w") as f:
        f.write("=" * 125 + "\n")
        f.write(f"{'Model Formulation':<24} {'lambda':<8} {'nOuter':<8} {'Outer Passes':<14} {'f_L(2m)':<10} {'E_net(J)':<12} {'E_stored(J)':<13} {'Delta E(J)':<12} {'Closure(%)':<12} {'CPU Time':<10}\n")
        f.write("-" * 125 + "\n")
        for r in sweep_results:
            lam_str = f"{r['lambda']:.1f}" if r['lambda'] is not None else "—"
            f.write(f"{r['model']:<24} {lam_str:<8} {r['n_outer']:<8} {r['total_outer_passes']:<14} {r['f_L']:<10.4f} {r['E_net']:<12.2f} {r['E_stored']:<13.2f} {r['surplus']:<+12.2f} {r['rel_err']:<+11.2f}% {r['cpu_time']:<9.2f}s\n")
        f.write("=" * 125 + "\n")
    print(f"\nSaved sweep results to:\n  - {sweep_json_path}\n  - {sweep_txt_path}")

    # 4. Perform Audits on All Four 17-minute Cases
    scorecard_cases = [
        {"id": "ep_tutorial_n8", "name": "EP (Tutorial: lambda = 0.2, nOuter=8)", "dir": ep_n8_17m_dir},
        {"id": "ep_lam0.2_n12", "name": "EP (lambda = 0.2, nOuter=12)", "dir": os.path.join(STUDY_17MIN_DIR, "audit_17min_ep_n12")},
        {"id": "ep_lam1.0_n8", "name": "EP (Default: lambda = 1.0, nOuter=8)", "dir": os.path.join(STUDY_17MIN_DIR, "audit_17min_ep_lam1_n8")},
        {"id": "ehc_baseline_n8", "name": "EHCDarcy (Baseline: nOuter=8)", "dir": os.path.join(STUDY_17MIN_DIR, "audit_17min_ehc")},
    ]

    histories = {}
    scorecard_data = {}

    for sc in scorecard_cases:
        cid = sc["id"]
        cdir = sc["dir"]
        print(f"\nAuditing 17-minute case: {sc['name']} ...")
        audit_history = run_exact_audit(cdir)
        histories[cid] = audit_history

        # Extract values at target times: 2 min (120s), 6 min (360s), 10 min (600s), 17 min (1020s)
        def get_row(t_sec):
            for r in audit_history:
                if abs(r["t_sec"] - t_sec) < 1e-3:
                    return r
            return None

        r2 = get_row(120.0)
        r6 = get_row(360.0)
        r10 = get_row(600.0)
        r17 = get_row(1020.0)

        # Front deviations
        mean_2, max_2 = calculate_front_deviation(cdir, 120)
        mean_6, max_6 = calculate_front_deviation(cdir, 360)
        mean_10, max_10 = calculate_front_deviation(cdir, 600)
        mean_17, max_17 = calculate_front_deviation(cdir, 1020)

        scorecard_data[cid] = {
            "name": sc["name"],
            "closure_2m": r2["rel_err"] if r2 else None,
            "surplus_2m": r2["surplus"] if r2 else None,
            "closure_6m": r6["rel_err"] if r6 else None,
            "surplus_6m": r6["surplus"] if r6 else None,
            "closure_17m": r17["rel_err"] if r17 else None,
            "surplus_17m": r17["surplus"] if r17 else None,
            "mean_front_2m": mean_2,
            "mean_front_6m": mean_6,
            "mean_front_10m": mean_10,
            "mean_front_17m": mean_17,
            "max_front_17m": max_17,
            "f_L_17m": r17["f_L"] if r17 else None,
        }

    scorecard_json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scorecard_data.json")
    with open(scorecard_json_path, "w") as f:
        json.dump(scorecard_data, f, indent=2)

    history_json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "transient_audit_histories.json")
    with open(history_json_path, "w") as f:
        json.dump(histories, f, indent=2)

    print(f"\nSaved scorecard and transient data to:\n  - {scorecard_json_path}\n  - {history_json_path}")
    print("\nALL RUNS AND AUDITS COMPLETED SUCCESSFULLY!")

if __name__ == "__main__":
    main()
