#!/usr/bin/env python3
"""
Master verification gate runner for phaseChangeFoam repository.
Runs all verification test suites under verification/ and exits non-zero if any test suite fails.
"""

import os
import sys
import time
import subprocess

SUITES = [
    "testModelFeatures/runTest.py",
    "phaseChangeCoupled/runTest.py",
    "singleCellOvershoot/runTest.py",
    "leeModel/runTest.py",
    "constantSource/runTest.py",
    "phaseChange1D_fluidEP/runTest.py",
    "phaseChange1D/runTest.py",
    "phaseChange1D_fluxBC/runTest.py",
    "phaseChange1D_heatingCooling/runTest.py",
    "phaseChange1D_porosity/runTest.py",
    "phaseChange1D_unequalCp/runTest.py",
    "phaseChange1D_variableK/runTest.py",
]

def main():
    verification_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(verification_dir)

    print("=" * 80)
    print("      PHASECHANGEFOAM MASTER VERIFICATION GATE RUNNER")
    print("=" * 80)
    print(f"Repository Root: {repo_dir}")
    print(f"Total Suites Registered: {len(SUITES)}\n")

    start_all_time = time.time()
    results = []

    for idx, suite_rel_path in enumerate(SUITES, 1):
        suite_abs_path = os.path.join(verification_dir, suite_rel_path)
        suite_dir = os.path.dirname(suite_abs_path)
        suite_name = os.path.basename(suite_dir)

        if not os.path.exists(suite_abs_path):
            print(f"[{idx}/{len(SUITES)}] ERROR: Test suite not found at {suite_abs_path}")
            results.append((suite_name, "MISSING", 0.0))
            continue

        print("-" * 80)
        print(f"[{idx}/{len(SUITES)}] Running suite: {suite_name} ({suite_rel_path})")
        print("-" * 80)

        suite_start = time.time()
        res = subprocess.run(
            [sys.executable, suite_abs_path],
            cwd=suite_dir,
            stdout=sys.stdout,
            stderr=sys.stderr,
            text=True
        )
        duration = time.time() - suite_start

        if res.returncode == 0:
            print(f"\n---> {suite_name}: PASSED ({duration:.2f}s)\n")
            results.append((suite_name, "PASSED", duration))
        else:
            print(f"\n---> {suite_name}: FAILED (exit code {res.returncode}, {duration:.2f}s)\n")
            results.append((suite_name, "FAILED", duration))

    total_duration = time.time() - start_all_time

    print("=" * 80)
    print("      VERIFICATION GATE SUMMARY REPORT")
    print("=" * 80)
    all_passed = True
    for suite_name, status, duration in results:
        status_str = f"PASSED ({duration:.2f}s)" if status == "PASSED" else f"FAILED ({duration:.2f}s)"
        print(f"  * {suite_name:<35} : {status_str}")
        if status != "PASSED":
            all_passed = False

    print("-" * 80)
    print(f"Total Verification Execution Time: {total_duration:.2f}s")

    if all_passed:
        print("\nALL VERIFICATION SUITES PASSED SUCCESSFULLY! (GATE GREEN)")
        sys.exit(0)
    else:
        print("\nONE OR MORE VERIFICATION SUITES FAILED! (GATE RED)")
        sys.exit(1)

if __name__ == "__main__":
    main()
