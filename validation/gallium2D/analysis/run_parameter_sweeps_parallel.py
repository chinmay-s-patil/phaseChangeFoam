#!/usr/bin/env python3
"""
Runs Cu sweep and deltaT sweep cases in parallel using available CPU cores.
"""

import os
import sys
import subprocess

def launch(case_dir):
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    return subprocess.Popen(cmd, shell=True)

if __name__ == "__main__":
    study_root = os.path.abspath("validation/gallium2D/analysis/study_17min")
    cases = ["cu_1e4", "cu_1e5", "cu_1e7", "cu_1e8", "deltaT_0.1", "deltaT_0.2"]
    procs = []
    for c in cases:
        cdir = os.path.join(study_root, c)
        print(f"Launching {c} ...")
        p = launch(cdir)
        procs.append((c, p))

    for c, p in procs:
        p.wait()
        print(f"Case {c} finished with exit code {p.returncode}")
