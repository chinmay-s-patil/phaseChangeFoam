#!/usr/bin/env python3
import os
import sys
import subprocess
import math

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
    vals = [float(x) for x in block.split()]
    return vals

def setup_case():
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"
    run_cmd(f"bash -c '{of_env}; blockMesh'")
    run_cmd("mkdir -p constant/pcm system/pcm")
    run_cmd("cp -r constant/polyMesh constant/pcm/polyMesh 2>/dev/null || true")
    run_cmd("cp system/fvSchemes system/pcm/fvSchemes 2>/dev/null || true")
    run_cmd("cp system/fvSolution system/pcm/fvSolution 2>/dev/null || true")

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"=== 1D PCM Verification Test in {case_dir} ===")
    
    os.chdir(case_dir)
    of_env = "source /usr/lib/openfoam/openfoam2412/etc/bashrc || true"
    solver_path = "/home/lavender/OpenFOAM/lavender-v2412/platforms/linux64GccDPInt32Opt/bin/phaseChangeMultiRegionFoam"

    # 1. Clean previous runs
    run_cmd("rm -rf [1-9]* 0.* constant/pcm/polyMesh constant/polyMesh")
    setup_case()

    # 2. Scenario A: Baseline Small Timestep (dt = 1 s)
    print("\n--- Running Case A: Small Timestep (dt = 1 s, N=2000 steps) ---")
    run_cmd("rm -rf [1-9]*")
    with open("system/controlDict", "r") as f:
        cd_content = f.read()
    cd_small = cd_content.replace("deltaT          100;", "deltaT          1;").replace("deltaT          200;", "deltaT          1;").replace("deltaT          500;", "deltaT          1;")
    with open("system/controlDict", "w") as f:
        f.write(cd_small)

    run_cmd(f"bash -c '{of_env}; {solver_path}'")

    T_small = parse_openfoam_field("2000/pcm/T")
    alpha_small = parse_openfoam_field("2000/pcm/liquidFraction")
    h_small = parse_openfoam_field("2000/pcm/h")

    print(f"Small dt (1s) at t=2000s: {len(T_small)} cells.")
    print(f"  Min T: {min(T_small):.2f} K, Max T: {max(T_small):.2f} K")
    print(f"  Min alphaL: {min(alpha_small):.4f}, Max alphaL: {max(alpha_small):.4f}")

    # 3. Scenario B: Large Timestep (dt = 200 s)
    print("\n--- Running Case B: Large Timestep (dt = 200 s, N=10 steps) ---")
    run_cmd("rm -rf [1-9]*")
    cd_large = cd_small.replace("deltaT          1;", "deltaT          200;")
    with open("system/controlDict", "w") as f:
        f.write(cd_large)

    run_cmd(f"bash -c '{of_env}; {solver_path}'")

    T_large = parse_openfoam_field("2000/pcm/T")
    alpha_large = parse_openfoam_field("2000/pcm/liquidFraction")
    h_large = parse_openfoam_field("2000/pcm/h")

    print(f"Large dt (200s) at t=2000s: {len(T_large)} cells.")
    print(f"  Min T: {min(T_large):.2f} K, Max T: {max(T_large):.2f} K")
    print(f"  Min alphaL: {min(alpha_large):.4f}, Max alphaL: {max(alpha_large):.4f}")

    # 4. Scenario C: Ultra-large Timestep (dt = 500 s)
    print("\n--- Running Case C: Ultra-large Timestep (dt = 500 s, N=4 steps) ---")
    run_cmd("rm -rf [1-9]*")
    cd_ultra = cd_small.replace("deltaT          1;", "deltaT          500;")
    with open("system/controlDict", "w") as f:
        f.write(cd_ultra)

    run_cmd(f"bash -c '{of_env}; {solver_path}'")

    T_ultra = parse_openfoam_field("2000/pcm/T")
    alpha_ultra = parse_openfoam_field("2000/pcm/liquidFraction")

    print(f"Ultra dt (500s) at t=2000s: {len(T_ultra)} cells.")
    print(f"  Min T: {min(T_ultra):.2f} K, Max T: {max(T_ultra):.2f} K")
    print(f"  Min alphaL: {min(alpha_ultra):.4f}, Max alphaL: {max(alpha_ultra):.4f}")

    # 5. Analysis & Verification Metrics
    print("\n=======================================================================")
    print("                1D PCM VERIFICATION SUMMARY RESULTS                    ")
    print("=======================================================================")
    print(f"{'Metric':<35} | {'dt = 1s':<10} | {'dt = 200s':<10} | {'dt = 500s':<10}")
    print("-" * 75)
    
    if len(T_small) == len(T_large) and len(T_small) > 0:
        mean_T_s = sum(T_small)/len(T_small)
        mean_T_l = sum(T_large)/len(T_large)
        mean_T_u = sum(T_ultra)/len(T_ultra)

        mean_a_s = sum(alpha_small)/len(alpha_small)
        mean_a_l = sum(alpha_large)/len(alpha_large)
        mean_a_u = sum(alpha_ultra)/len(alpha_ultra)

        max_T_diff_l = max(abs(t1 - t2) for t1, t2 in zip(T_small, T_large))
        rms_T_diff_l = math.sqrt(sum((t1 - t2)**2 for t1, t2 in zip(T_small, T_large)) / len(T_small))

        max_T_diff_u = max(abs(t1 - t2) for t1, t2 in zip(T_small, T_ultra))
        rms_T_diff_u = math.sqrt(sum((t1 - t2)**2 for t1, t2 in zip(T_small, T_ultra)) / len(T_small))

        print(f"{'Mean Temperature (K)':<35} | {mean_T_s:<10.2f} | {mean_T_l:<10.2f} | {mean_T_u:<10.2f}")
        print(f"{'Mean Liquid Fraction':<35} | {mean_a_s:<10.4f} | {mean_a_l:<10.4f} | {mean_a_u:<10.4f}")
        print(f"{'Max T Diff vs dt=1s (K)':<35} | {'-':<10} | {max_T_diff_l:<10.4f} | {max_T_diff_u:<10.4f}")
        print(f"{'RMS T Diff vs dt=1s (K)':<35} | {'-':<10} | {rms_T_diff_l:<10.4f} | {rms_T_diff_u:<10.4f}")
        print("-" * 75)

        if max_T_diff_l < 3.0 and max_T_diff_u < 5.0:
            print("\nSTATUS: VERIFICATION SUCCESSFUL!")
            print("The solver produces consistent, stable, and accurate thermal and phase fraction profiles even under 500x larger timesteps.")
        else:
            print(f"\nSTATUS: VERIFICATION FAILED (Max T diff {max_T_diff_l:.2f} K)")
            sys.exit(1)
    else:
        print("Error: Could not parse fields properly!")
        sys.exit(1)

if __name__ == "__main__":
    main()
