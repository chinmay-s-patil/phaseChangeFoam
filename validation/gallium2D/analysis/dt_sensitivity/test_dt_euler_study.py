#!/usr/bin/env python3
"""
Tests the Brent et al. (1988) large-dt hypothesis:
"Did Brent's use of dt = 5-10s with first-order Euler ddt scheme cause low f_L(2 min) = 0.1102?"
Compares:
  - dt = 0.1s (baseline OpenFOAM)
  - dt = 1.0s
  - dt = 5.0s
  - dt = 5s (0-20s) -> 10s (20-120s) (Brent 1988 exact time-stepping schedule)
All using Euler ddtSchemes on the 42x32 mesh to t = 120s.
"""

import os
import sys
import shutil
import subprocess
import numpy as np

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing: {cmd}\nSTDERR: {res.stderr}\nSTDOUT: {res.stdout}")
        sys.exit(1)
    return res

def parse_phase_fraction(filepath):
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return None
    start_paren = content.find("(", idx)
    end_paren = content.rfind(")", 0, content.find(";", start_paren))
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def setup_dt_case(base_dir, case_dir, dt=0.1, end_time=120):
    os.makedirs(case_dir, exist_ok=True)
    for folder in ["0", "constant", "system"]:
        src = os.path.join(base_dir, folder)
        dst = os.path.join(case_dir, folder)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

    ctrl_path = os.path.join(case_dir, "system", "controlDict")
    with open(ctrl_path, "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT {dt};
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
""")

def run_case(case_dir):
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    run_cmd(cmd)

def run_brent_schedule(base_dir, case_dir):
    """Runs 0->20s with dt=5s, then 20->120s with dt=10s"""
    setup_dt_case(base_dir, case_dir, dt=5.0, end_time=20)
    # Set writeInterval to 20 so it writes at t=20
    ctrl_path = os.path.join(case_dir, "system", "controlDict")
    with open(ctrl_path, "r") as f:
        c = f.read()
    c = c.replace("writeInterval 60;", "writeInterval 20;")
    with open(ctrl_path, "w") as f:
        f.write(c)

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    # Run stage 1 (0 -> 20s)
    cmd1 = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.stage1 2>&1'"
    run_cmd(cmd1)

    # Modify controlDict for stage 2 (20 -> 120s with dt=10s)
    with open(ctrl_path, "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom latestTime;
stopAt endTime;
endTime 120;
deltaT 10;
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
""")
    # Run stage 2
    cmd2 = f"bash -c '{of_env}; cd {case_dir} && phaseChangeMultiRegionFoam >> log.stage2 2>&1'"
    run_cmd(cmd2)

if __name__ == "__main__":
    base = os.path.abspath("validation/gallium2D/case")
    out_dir = os.path.abspath("validation/gallium2D/analysis/dt_sensitivity")
    os.makedirs(out_dir, exist_ok=True)

    dt_list = [0.1, 1.0, 5.0]
    results = {}

    for dt in dt_list:
        c_dir = os.path.join(out_dir, f"case_dt_{dt}")
        print(f"Running dt = {dt} s to 120 s ...")
        setup_dt_case(base, c_dir, dt=dt, end_time=120)
        run_case(c_dir)
        pf = parse_phase_fraction(os.path.join(c_dir, "120", "fluid", "phaseFraction"))
        if pf is not None:
            results[f"dt = {dt}s"] = float(np.mean(pf))
            print(f"  --> f_L(120s) = {results[f'dt = {dt}s']:.4f}")

    # Brent schedule
    brent_dir = os.path.join(out_dir, "case_dt_brent_schedule")
    print("Running Brent schedule (dt=5s to 20s, then dt=10s to 120s) ...")
    run_brent_schedule(base, brent_dir)
    pf_brent = parse_phase_fraction(os.path.join(brent_dir, "120", "fluid", "phaseFraction"))
    if pf_brent is not None:
        results["Brent schedule (5s->10s)"] = float(np.mean(pf_brent))
        print(f"  --> f_L(120s) = {results['Brent schedule (5s->10s)']:.4f}")

    print("\n" + "=" * 65)
    print("      TIME-STEPPING & BRENT DT HYPOTHESIS TEST RESULTS")
    print("=" * 65)
    print(f"{'Condition':<30} {'f_L(2 min)':<15} {'Brent (1988)':<15}")
    print("-" * 65)
    for k, v in results.items():
        print(f"{k:<30} {v:<15.4f} {0.1102:<15.4f}")
    print("=" * 65)
