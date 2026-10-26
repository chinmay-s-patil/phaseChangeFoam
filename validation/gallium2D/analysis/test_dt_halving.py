#!/usr/bin/env python3
"""
Test script to investigate the EHCDarcy conservation deficit under dt halving:
Runs EHCDarcy at dt = 0.05s to t = 120s (2 min) and runs exact right-endpoint audit.
Compares with dt = 0.10s (baseline: Delta E = -6.43 J, closure = -0.34%).
"""

import os
import sys
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GALLIUM_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, GALLIUM_DIR)
from runTest import find_solver

sys.path.insert(0, os.path.join(SCRIPT_DIR, "energy_balance"))
from exact_energy_balance_audit import run_exact_audit

def run_test():
    case_dir = os.path.join(SCRIPT_DIR, "case_dt_halving_0.05")
    base_case = os.path.join(GALLIUM_DIR, "analysis", "energy_balance", "sweep_runs", "ehc_n8")
    
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir, exist_ok=True)

    for folder in ["0", "constant", "system"]:
        shutil.copytree(os.path.join(base_case, folder), os.path.join(case_dir, folder))

    # Update controlDict to dt = 0.05s
    ctrl_path = os.path.join(case_dir, "system", "controlDict")
    with open(ctrl_path, "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 120;
deltaT 0.05;
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
{
    wallHeatFlux1
    {
        type            wallHeatFlux;
        libs            (fieldFunctionObjects);
        region          fluid;
        writeControl    timeStep;
        writeInterval   1;
        patches         (hotWall coldWall);
    }
}
""")

    solver = find_solver()
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* && blockMesh && mkdir -p constant/fluid && cp -r constant/polyMesh constant/fluid/polyMesh && {solver} > log.simulation 2>&1'"
    
    print(f"Running EHCDarcy at dt = 0.05s to 120s (2400 steps)...")
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Simulation failed:\n{res.stderr}")
        return

    print("Simulation complete. Running exact right-endpoint energy audit...")
    results = run_exact_audit(case_dir)
    for r in results:
        if abs(r['t_sec'] - 120.0) < 1e-3:
            print("\n" + "="*70)
            print(f"EHCDarcy (dt = 0.05s, nOuter=8) at t = 120s (2 min):")
            print(f"  E_hot:     {r['E_in_hot']:.2f} J")
            print(f"  E_cold:    {r['E_out_cold']:.2f} J")
            print(f"  E_net:     {r['E_net']:.2f} J")
            print(f"  E_stored:  {r['E_stored']:.2f} J")
            print(f"  Delta E:   {r['surplus']:+.2f} J")
            print(f"  Closure:   {r['rel_err']:+.3f}%")
            print("="*70)

if __name__ == "__main__":
    run_test()
