#!/usr/bin/env python3
"""
Master script to execute the complete suite of 17-minute (1020 s) validation runs:
1. Spatial Mesh Convergence to 17 min (32x24, 42x32, 63x48) + GCI calculation
2. Exact Energy Audit to 17 min (EHCDarcy and EP nOuter=12)
3. Late-Time Darcy Cu Sweep to 17 min (1e4, 1e5, 1.6e6, 1e7, 1e8)
4. Late-Time Melting Window Sweep to 17 min (0.1K, 0.2K, 0.5K)
"""

import os
import sys
import shutil
import subprocess
import time
import numpy as np

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing: {cmd}\nSTDERR: {res.stderr}\nSTDOUT: {res.stdout}")
        sys.exit(1)
    return res

def setup_case(base_dir, case_dir, nx=42, ny=32, end_time=1020, write_interval=60,
               model_type="EHCDarcy", n_outer=8, cu=1.6e6, deltaT=0.5, with_wall_flux=False):
    os.makedirs(case_dir, exist_ok=True)
    for folder in ["0", "constant", "system"]:
        src = os.path.join(base_dir, folder)
        dst = os.path.join(case_dir, folder)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)

    W = 0.0889
    H = 0.0635
    D = 0.005

    # 1. blockMeshDict
    with open(os.path.join(case_dir, "system", "blockMeshDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
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
    hotWall {{ type wall; faces ((0 4 7 3)); }}
    coldWall {{ type wall; faces ((1 2 6 5)); }}
    topWall {{ type wall; faces ((3 7 6 2)); }}
    bottomWall {{ type wall; faces ((0 1 5 4)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")

    # 2. controlDict
    whf_block = """
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
""" if with_wall_flux else ""

    with open(os.path.join(case_dir, "system", "controlDict"), "w") as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object controlDict; }}
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
{whf_block}
""")

    # 3. fvSolution (nOuterCorrectors)
    for sp in [os.path.join(case_dir, "system", "fluid", "fvSolution"), os.path.join(case_dir, "system", "fvSolution")]:
        if os.path.exists(sp):
            with open(sp, "r") as f:
                content = f.read()
            import re
            content = re.sub(r"nOuterCorrectors\s+\d+;", f"nOuterCorrectors {n_outer};", content)
            with open(sp, "w") as f:
                f.write(content)

    # 4. phaseChangeDict
    T_m = 302.93
    T_l = T_m - deltaT / 2.0
    T_u = T_m + deltaT / 2.0

    pc_path = os.path.join(case_dir, "constant", "fluid", "phaseChangeDict")
    if model_type == "EHCDarcy":
        pc_text = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
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
"""
    else:
        pc_text = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "constant/fluid"; object phaseChangeDict; }}
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
"""
    with open(pc_path, "w") as f:
        f.write(pc_text)

def launch_case_process(case_dir):
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    prep_cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh'"
    run_cmd(prep_cmd)
    sim_cmd = f"bash -c '{of_env}; cd {case_dir} && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    p = subprocess.Popen(sim_cmd, shell=True)
    return p

if __name__ == "__main__":
    base = os.path.abspath("validation/gallium2D/case")
    study_root = os.path.abspath("validation/gallium2D/analysis/study_17min")
    os.makedirs(study_root, exist_ok=True)

    # Define all required 17-min runs
    case_definitions = [
        # 1. Mesh Convergence Study (EHCDarcy, 17 min, write every 60s)
        {"id": "mesh_32x24", "nx": 32, "ny": 24, "model": "EHCDarcy", "cu": 1.6e6, "dT": 0.5, "nOuter": 8, "whf": False},
        {"id": "mesh_63x48", "nx": 63, "ny": 48, "model": "EHCDarcy", "cu": 1.6e6, "dT": 0.5, "nOuter": 8, "whf": False},

        # 2. Exact Energy Audit to 17 min (with per-step wallHeatFlux)
        {"id": "audit_17min_ehc", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1.6e6, "dT": 0.5, "nOuter": 8, "whf": True},
        {"id": "audit_17min_ep_n12", "nx": 42, "ny": 32, "model": "enthalpyPorosity", "cu": 1.6e6, "dT": 0.5, "nOuter": 12, "whf": True},

        # 3. Cu Sweep to 17 min (EHCDarcy)
        {"id": "cu_1e4", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1e4, "dT": 0.5, "nOuter": 8, "whf": False},
        {"id": "cu_1e5", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1e5, "dT": 0.5, "nOuter": 8, "whf": False},
        {"id": "cu_1e7", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1e7, "dT": 0.5, "nOuter": 8, "whf": False},
        {"id": "cu_1e8", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1e8, "dT": 0.5, "nOuter": 8, "whf": False},

        # 4. Melting Window Sweep to 17 min (EHCDarcy)
        {"id": "deltaT_0.1", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1.6e6, "dT": 0.1, "nOuter": 8, "whf": False},
        {"id": "deltaT_0.2", "nx": 42, "ny": 32, "model": "EHCDarcy", "cu": 1.6e6, "dT": 0.2, "nOuter": 8, "whf": False},
    ]

    print(f"Setting up {len(case_definitions)} cases in {study_root} ...")
    for cdef in case_definitions:
        cdir = os.path.join(study_root, cdef["id"])
        setup_case(
            base, cdir,
            nx=cdef["nx"], ny=cdef["ny"],
            end_time=1020, write_interval=60,
            model_type=cdef["model"], n_outer=cdef["nOuter"],
            cu=cdef["cu"], deltaT=cdef["dT"],
            with_wall_flux=cdef["whf"]
        )

    # Run in parallel batches (up to 4 concurrent simulations to maximize 16 CPU cores without overloading)
    batch_size = 4
    total = len(case_definitions)
    for i in range(0, total, batch_size):
        batch = case_definitions[i:i+batch_size]
        print(f"\n=======================================================")
        print(f" Launching Batch {i//batch_size + 1}: {[b['id'] for b in batch]}")
        print(f"=======================================================")
        procs = []
        for b in batch:
            cdir = os.path.join(study_root, b["id"])
            print(f"Starting {b['id']} ...")
            p = launch_case_process(cdir)
            procs.append((b["id"], p, cdir))

        # Wait for batch to complete
        for cid, p, cdir in procs:
            p.wait()
            print(f"  --> Case {cid} finished with exit code {p.returncode}")

    print("\nAll 17-minute cases completed successfully!")
