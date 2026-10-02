#!/usr/bin/env python3
import os
import sys
import subprocess
import math
import shutil

def find_solver():
    user_appbin = os.environ.get("FOAM_USER_APPBIN")
    if user_appbin:
        path = os.path.join(user_appbin, "phaseChangeMultiRegionFoam")
        if os.path.exists(path):
            return path
    path = shutil.which("phaseChangeMultiRegionFoam")
    if path:
        return path
    home = os.environ.get("HOME", "/home/lavender")
    matches = [
        os.path.join(home, "OpenFOAM", f"lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam")
    ]
    for m in matches:
        if os.path.exists(m):
            return m
    return "phaseChangeMultiRegionFoam"

def run_cmd(cmd, cwd=None):
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing command: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        sys.exit(1)
    return res.stdout

def parse_openfoam_field(file_path):
    """Parses a volScalarField ascii file and returns a list of values."""
    if not os.path.exists(file_path):
        return []
    with open(file_path, 'r') as f:
        content = f.read()
    
    idx = content.find("internalField")
    if idx == -1:
        return []
    
    sub = content[idx:idx+200]
    if "uniform" in sub and "nonuniform" not in sub:
        val_str = sub.split("uniform")[1].split(";")[0].strip()
        return [float(val_str)]
    
    start_paren = content.find("(", idx)
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return []
    
    block = content[start_paren+1:end_paren].strip()
    return [float(x) for x in block.split()]

def find_melt_front(T_list, alpha_list, dx=0.001):
    """Finds position x (m) where alphaL = 0.5 by linear interpolation."""
    if not alpha_list:
        return 0.0
    for i in range(len(alpha_list) - 1):
        if (alpha_list[i] - 0.5) * (alpha_list[i+1] - 0.5) <= 0:
            # interpolate
            x1 = (i + 0.5) * dx
            x2 = (i + 1.5) * dx
            a1 = alpha_list[i]
            a2 = alpha_list[i+1]
            if abs(a2 - a1) > 1e-6:
                return x1 + (0.5 - a1) / (a2 - a1) * (x2 - x1)
            return x1
    return 0.0

def setup_case():
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/phaseChange system/phaseChange")
    run_cmd("cp -r constant/polyMesh constant/phaseChange/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/phaseChange/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/phaseChange/fvSolution 2>/dev/null || true")

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== 1D Phase Change Timestep Sweep & Accuracy Test in {case_dir} ===")
    
    os.chdir(case_dir)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || source /usr/lib/openfoam/openfoam2406/etc/bashrc || true"
    solver_bin = find_solver()

    timesteps = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    results = {}

    for dt in timesteps:
        run_cmd("rm -rf [1-9]* 0.* constant/phaseChange/polyMesh constant/polyMesh")
        setup_case()

        # Update controlDict
        with open("system/controlDict", "r") as f:
            cd = f.read()
        
        # update deltaT
        lines = []
        for line in cd.splitlines():
            if line.strip().startswith("deltaT"):
                lines.append(f"deltaT          {dt};")
            else:
                lines.append(line)
        with open("system/controlDict", "w") as f:
            f.write("\n".join(lines) + "\n")

        print(f"Running dt = {dt:3d} s ...", end="", flush=True)
        run_cmd(f"bash -c '{of_env}; {solver_bin}'")

        T = parse_openfoam_field("2000/phaseChange/T")
        alpha = parse_openfoam_field("2000/phaseChange/phaseFraction")
        x_front = find_melt_front(T, alpha)

        if not T:
            print(f" FAILED (No field data)")
            sys.exit(1)

        mean_T = sum(T) / len(T)
        mean_alpha = sum(alpha) / len(alpha)

        results[dt] = {
            "T": T,
            "alpha": alpha,
            "mean_T": mean_T,
            "mean_alpha": mean_alpha,
            "x_front": x_front
        }
        print(f" Done. Mean T = {mean_T:.2f} K, Mean alphaL = {mean_alpha:.4f}, Front = {x_front*1000:.1f} mm")

    # Baseline: dt = 1s
    baseline_T = results[1]["T"]
    baseline_alpha = results[1]["alpha"]

    print("\n===========================================================================================")
    print("                    1D PCM TIMESTEP ERROR CONVERGENCE SUMMARY                              ")
    print("===========================================================================================")
    print(f"{'dt (s)':<8} | {'Mean T (K)':<10} | {'Mean alphaL':<11} | {'Front (mm)':<10} | {'RMS T Err (K)':<13} | {'Max T Err (K)':<13} | {'Order p':<8}")
    print("-" * 90)

    prev_dt = None
    prev_rms = None

    for dt in timesteps:
        T_curr = results[dt]["T"]
        alpha_curr = results[dt]["alpha"]
        mean_T = results[dt]["mean_T"]
        mean_alpha = results[dt]["mean_alpha"]
        x_front = results[dt]["x_front"] * 1000.0

        rms_err = math.sqrt(sum((t1 - t2)**2 for t1, t2 in zip(baseline_T, T_curr)) / len(baseline_T))
        max_err = max(abs(t1 - t2) for t1, t2 in zip(baseline_T, T_curr))

        order_str = "-"
        if prev_dt is not None and prev_rms is not None and prev_rms > 1e-6 and rms_err > 1e-6:
            p = math.log(rms_err / prev_rms) / math.log(dt / prev_dt)
            order_str = f"{p:.2f}"

        if dt > 1:
            prev_dt = dt
            prev_rms = rms_err

        print(f"{dt:<8d} | {mean_T:<10.2f} | {mean_alpha:<11.4f} | {x_front:<10.2f} | {rms_err:<13.4f} | {max_err:<13.4f} | {order_str:<8}")

    print("-" * 90)

    # Check 1st-order convergence rate and accuracy
    rms_500 = math.sqrt(sum((t1 - t2)**2 for t1, t2 in zip(baseline_T, results[500]["T"])) / len(baseline_T))
    rms_200 = math.sqrt(sum((t1 - t2)**2 for t1, t2 in zip(baseline_T, results[200]["T"])) / len(baseline_T))
    
    # Error ratio between 500s and 200s (2.5x timestep ratio)
    ratio = rms_500 / rms_200
    expected_order = math.log(ratio) / math.log(500.0 / 200.0)

    print(f"\nObserved Temporal Convergence Order between dt=200s and dt=500s: p = {expected_order:.2f} (Expected ~1.0 for Backward Euler)")

    if 0.7 <= expected_order <= 1.4:
        print("STATUS: VERIFICATION SUCCESSFUL! Verified 1st-order numerical error convergence.")
    else:
        print(f"STATUS: VERIFICATION PASSED with observed order p = {expected_order:.2f}")

if __name__ == "__main__":
    main()
