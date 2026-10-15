#!/usr/bin/env python3
"""
Runs exact energy conservation audit comparison for:
1. EHCDarcy (baseline)
2. EP with nOuterCorr = 2
3. EP with nOuterCorr = 5
4. EP with nOuterCorr = 8 (baseline)
5. EP with nOuterCorr = 12
All cases run to t = 120s with per-step wallHeatFlux logging.
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

def setup_case(base_dir, case_dir, model_type="enthalpyPorosity", n_outer=8, end_time=120):
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

    # 2. Modify fvSolution for nOuterCorrectors
    sol_paths = [os.path.join(case_dir, "system", "fluid", "fvSolution"), os.path.join(case_dir, "system", "fvSolution")]
    for sp in sol_paths:
        if os.path.exists(sp):
            with open(sp, "r") as f:
                content = f.read()
            import re
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
        # enthalpyPorosity
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
    lambda 0.2;
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

def run_case(case_dir):
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    print(f"Running simulation in {case_dir} ...")
    run_cmd(cmd)
    print("Done.")

if __name__ == "__main__":
    base = os.path.abspath("validation/gallium2D/case")
    audit_script = os.path.abspath("validation/gallium2D/analysis/energy_balance/exact_energy_balance_audit.py")

    cases = [
        ("case_audit_ehcDarcy", "EHCDarcy", 8),
        ("case_audit_ep_nOuter2", "enthalpyPorosity", 2),
        ("case_audit_ep_nOuter5", "enthalpyPorosity", 5),
        ("case_audit_ep_nOuter8", "enthalpyPorosity", 8),
        ("case_audit_ep_nOuter12", "enthalpyPorosity", 12),
    ]

    analysis_root = os.path.abspath("validation/gallium2D/analysis/energy_balance")

    for cname, mtype, nouter in cases:
        cdir = os.path.join(analysis_root, cname)
        print(f"\n---> Setting up {cname} (Model: {mtype}, nOuterCorr: {nouter})")
        setup_case(base, cdir, model_type=mtype, n_outer=nouter, end_time=120)
        run_case(cdir)
        print(f"---> Running exact audit for {cname}:")
        run_cmd(f"python3 {audit_script} {cdir}")
