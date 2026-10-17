#!/usr/bin/env python3
"""
Runs Enthalpy-Porosity with code-default lambda = 1.0 and nOuterCorrectors = 8
to full 17 min (1020 s) with per-step wallHeatFlux logging.
"""

import os
import sys
import shutil
import subprocess

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing: {cmd}\nSTDERR: {res.stderr}\nSTDOUT: {res.stdout}")
        sys.exit(1)
    return res

def main():
    base_dir = os.path.abspath("validation/gallium2D/case")
    case_dir = os.path.abspath("validation/gallium2D/analysis/study_17min/audit_17min_ep_lam1_n8")
    os.makedirs(case_dir, exist_ok=True)

    for folder in ["0", "constant", "system"]:
        src = os.path.join(base_dir, folder)
        dst = os.path.join(case_dir, folder)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

    # 1. Modify controlDict: end_time = 1020, writeInterval = 60
    ctrl_path = os.path.join(case_dir, "system", "controlDict")
    with open(ctrl_path, "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "system"; object controlDict; }
application phaseChangeMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1020;
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

    # 2. Modify fvSolution for nOuterCorrectors = 8
    import re
    sol_paths = [os.path.join(case_dir, "system", "fluid", "fvSolution"), os.path.join(case_dir, "system", "fvSolution")]
    for sp in sol_paths:
        if os.path.exists(sp):
            with open(sp, "r") as f:
                content = f.read()
            content = re.sub(r"nOuterCorrectors\s+\d+;", "nOuterCorrectors 8;", content)
            with open(sp, "w") as f:
                f.write(content)

    # 3. Modify phaseChangeDict: lambda = 1.0
    pc_path = os.path.join(case_dir, "constant", "fluid", "phaseChangeDict")
    with open(pc_path, "w") as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }
active true;
phaseChange
{
    type enthalpyPorosity;
    active true;
    forward
    {
        T_lowerBound 302.68;
        T_upperBound 303.18;
        latentHeat 80160.0;
    }
    lambda 1.0;
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

    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    print(f"Running 17-minute EP lambda=1.0, nOuter=8 in {case_dir} ...")
    run_cmd(cmd)
    print("Simulation complete. Running audit...")
    audit_script = os.path.abspath("validation/gallium2D/analysis/energy_balance/exact_energy_balance_audit.py")
    run_cmd(f"python3 {audit_script} {case_dir}")

if __name__ == "__main__":
    main()
