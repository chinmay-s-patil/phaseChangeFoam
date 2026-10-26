#!/usr/bin/env python3
"""
Study of Enthalpy-Porosity Newton Relaxation Factor lambda and nOuterCorrectors.
Evaluates lambda in {1.0, 0.5, 0.2} across nOuterCorrectors in {2, 3, 5, 8} at t = 120s (2 min).
Computes exact step-resolved energy closure (E_net = E_hot - E_cold_out),
liquid fraction f_L, stored enthalpy, and run time.
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

def parse_openfoam_field(filepath):
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
    end_paren = content.rfind(")", 0, content.find(";", start_paren))
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def parse_wall_heat_flux(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"wallHeatFlux file not found at: {filepath}")
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                t = float(parts[0])
                patch = parts[1]
                integral = float(parts[4])
                if patch not in data:
                    data[patch] = []
                data[patch].append((t, integral))
    for p in data:
        data[p].sort(key=lambda x: x[0])
    return data

def integrate_flux(time_flux_list):
    if not time_flux_list:
        return np.array([]), np.array([])
    times = np.array([x[0] for x in time_flux_list])
    fluxes = np.array([x[1] for x in time_flux_list])
    cum_energy = np.zeros(len(times))
    dt0 = times[0]
    cum_energy[0] = fluxes[0] * dt0
    for i in range(1, len(times)):
        dt = times[i] - times[i-1]
        cum_energy[i] = cum_energy[i-1] + 0.5 * (fluxes[i-1] + fluxes[i]) * dt
    return times, cum_energy

def setup_case(base_dir, case_dir, lam=1.0, n_outer=2, end_time=120):
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

def run_case(case_dir):
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    cmd = f"bash -c '{of_env}; cd {case_dir} && rm -rf [1-9]* 0.* postProcessing && blockMesh && mkdir -p constant/fluid && rm -rf constant/fluid/polyMesh && cp -r constant/polyMesh constant/fluid/polyMesh && phaseChangeMultiRegionFoam > log.simulation 2>&1'"
    t0 = time.time()
    run_cmd(cmd)
    t1 = time.time()
    return t1 - t0

def audit_case(case_dir, target_t=120.0, nx=42, ny=32):
    W = 0.0889
    H = 0.0635
    D = 0.005
    dx = W / nx
    dy = H / ny
    V_cell = dx * dy * D
    rho = 6095.0
    Cp = 381.5
    L = 80160.0
    T_init = 301.45

    whf_file = os.path.join(case_dir, "postProcessing", "fluid", "wallHeatFlux1", "0", "wallHeatFlux.dat")
    if not os.path.exists(whf_file):
        whf_file = os.path.join(case_dir, "postProcessing", "wallHeatFlux1", "0", "wallHeatFlux.dat")
    
    flux_data = parse_wall_heat_flux(whf_file)
    hot_times, E_in_hot = integrate_flux(flux_data.get('hotWall', []))
    cold_times, E_cold_raw = integrate_flux(flux_data.get('coldWall', []))

    # Outward cold loss magnitude is -E_cold_raw
    E_out_cold = -E_cold_raw

    # Field storage at target_t
    t_str = f"{int(target_t)}" if target_t.is_integer() else f"{target_t}"
    T_path = os.path.join(case_dir, t_str, "fluid", "T")
    alpha_path = os.path.join(case_dir, t_str, "fluid", "phaseFraction")

    T = parse_openfoam_field(T_path)
    alpha = parse_openfoam_field(alpha_path)

    E_sens = float(np.sum(rho * V_cell * Cp * (T - T_init)))
    E_lat = float(np.sum(rho * V_cell * L * alpha))
    E_stored = E_sens + E_lat
    f_L = float(np.mean(alpha))

    idx = np.argmin(np.abs(hot_times - target_t))
    e_hot = E_in_hot[idx]
    e_cold_loss = E_out_cold[idx]
    e_net = e_hot - e_cold_loss

    delta_E = E_stored - e_net
    rel_err = (delta_E / e_net) * 100.0

    return {
        "t": target_t,
        "E_hot": e_hot,
        "E_cold_loss": e_cold_loss,
        "E_net": e_net,
        "E_sens": E_sens,
        "E_lat": E_lat,
        "E_stored": E_stored,
        "delta_E": delta_E,
        "rel_err": rel_err,
        "f_L": f_L
    }

def main():
    base_dir = os.path.abspath("validation/gallium2D/case")
    study_root = os.path.abspath("validation/gallium2D/analysis/energy_balance/ep_relaxation_study")
    os.makedirs(study_root, exist_ok=True)

    test_matrix = [
        # (lambda, n_outer)
        (1.0, 2),
        (1.0, 3),
        (1.0, 5),
        (0.5, 2),
        (0.5, 3),
        (0.5, 5),
        (0.2, 2),
        (0.2, 3),
        (0.2, 5),
        (0.2, 8),
        (1.0, 8),
    ]

    print("=" * 115)
    print(" ENTHALPY-POROSITY NEWTON RELAXATION (lambda) & OUTER-ITERATION (nOuter) STUDY (t = 120s)")
    print("=" * 115)

    results = []
    for lam, n_outer in test_matrix:
        case_name = f"ep_lam{lam:.1f}_nOuter{n_outer}"
        case_dir = os.path.join(study_root, case_name)
        print(f"\n---> Running {case_name}: lambda={lam}, nOuterCorrectors={n_outer}")
        setup_case(base_dir, case_dir, lam=lam, n_outer=n_outer, end_time=120)
        wall_time = run_case(case_dir)
        audit = audit_case(case_dir, target_t=120.0)
        audit["lam"] = lam
        audit["n_outer"] = n_outer
        audit["wall_time"] = wall_time
        results.append(audit)
        print(f"     Done in {wall_time:.1f}s: f_L={audit['f_L']:.4f}, E_net={audit['E_net']:.2f}J, E_stored={audit['E_stored']:.2f}J, Error={audit['rel_err']:+.2f}%")

    print("\n" + "=" * 115)
    print(" SUMMARY TABLE: ENTHALPY-POROSITY RELAXATION & OUTER CORRECTOR SWEEP (t = 120s / 2 min)")
    print("=" * 115)
    header = f"{'Model':<18} {'lambda':<8} {'nOuter':<8} {'f_L(2m)':<10} {'E_hot(J)':<12} {'E_cold(J)':<11} {'E_net(J)':<12} {'E_stored(J)':<13} {'Delta_E(J)':<12} {'Closure Err':<12} {'Clock Time':<10}"
    print(header)
    print("-" * 115)

    # Reference EHCDarcy numbers:
    # E_net = 1888.43 J (with startup 8.12J), E_stored = 1878.37 J, delta = -10.06 J (-0.53%)
    # Or without startup in E_net: 1880.31 J, delta = -1.94 J (-0.10%)
    lines = [header, "-" * 115]
    for r in results:
        line = f"{'enthalpyPorosity':<18} {r['lam']:<8.1f} {r['n_outer']:<8} {r['f_L']:<10.4f} {r['E_hot']:<12.2f} {r['E_cold_loss']:<11.2f} {r['E_net']:<12.2f} {r['E_stored']:<13.2f} {r['delta_E']:<+12.2f} {r['rel_err']:<+10.2f}% {r['wall_time']:<8.1f}s"
        print(line)
        lines.append(line)

    out_file = os.path.join(study_root, "ep_relaxation_study_results.txt")
    with open(out_file, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nResults written to: {out_file}")

if __name__ == "__main__":
    main()
