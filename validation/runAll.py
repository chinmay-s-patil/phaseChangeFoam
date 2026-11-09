#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Chinmay S Patil
"""
Master validation ladder runner for phaseChangeFoam repository.
Runs validation suites under validation/ and exits non-zero if any benchmark fails.

Usage:
    python3 validation/runAll.py          # Runs full validation ladder (All 5 rungs, ~3.5 min)
    python3 validation/runAll.py --fast   # Runs fast validation gate (Rungs 1 & 2: de Vahl Davis ~2s, EHC hysteresis ~2s)
    python3 validation/runAll.py --clean  # Cleans case directories after running
"""

import os
import sys
import time
import argparse
import subprocess

VALIDATION_SUITES = [
    {
        "name": "deVahlDavis2D",
        "title": "de Vahl Davis (1983) Natural Convection Cavity (Benchmark)",
        "script": "deVahlDavis2D/runTest.py",
        "fast": True,
        "description": "Steady natural convection cavity at Ra = 10^5, Nu_mean benchmark target 4.519 (+/- 5%) (~3s)"
    },
    {
        "name": "pcmHysteresis",
        "title": "EHC Thermal Hysteresis Cyclic Path Logic (Model Verification)",
        "script": "pcmHysteresis/runTest.py",
        "fast": True,
        "description": "Direction-dependent melting and freezing loop with 5 K hysteresis window (~2s)"
    },
    {
        "name": "tin2D",
        "title": "Bertrand / Hannoun Pure Tin Melting (Conduction Verified)",
        "script": "tin2D/runTest.py",
        "fast": False,
        "description": "Square cavity tin melting (Ra = 2.50e5, Pr = 0.02, Ste = 0.01) vs 1D Stefan solution (~15s)"
    },
    {
        "name": "gallium2D",
        "title": "Brent et al. (1988) / Ye & Arıcı (2023) Gallium Melting (Benchmark)",
        "script": "gallium2D/runTest.py",
        "fast": False,
        "description": "2-minute regression gate: f_L(2m) = 0.1303 +/- 0.002, |energy closure| <= 0.5% (~35s)"
    },
    {
        "name": "gallium3D_experimental",
        "title": "Gau & Viskanta (1986) 3D Gallium Melting (Validation in Progress)",
        "script": "gallium3D_experimental/runTest.py",
        "fast": False,
        "description": "3D cavity melting with physical no-slip sidewall braking and front position check (~2.5m)"
    },
]

def main():
    parser = argparse.ArgumentParser(description="Master Validation Ladder Runner for phaseChangeFoam")
    parser.add_argument("--fast", action="store_true", help="Run only fast benchmarks (Rungs 1 & 2: de Vahl Davis & EHC hysteresis)")
    parser.add_argument("--clean", action="store_true", help="Remove simulation directories after completion")
    args = parser.parse_args()

    validation_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(validation_dir)

    print("=" * 80)
    print("      PHASECHANGEFOAM NATURAL CONVECTION & PHASE CHANGE VALIDATION LADDER")
    print("=" * 80)
    print(f"Repository Root: {repo_dir}")
    print(f"Validation Root: {validation_dir}")
    mode_str = "FAST GATE (Rungs 1 & 2)" if args.fast else "FULL VALIDATION LADDER (All 5 rungs)"
    print(f"Execution Mode:  {mode_str}\n")

    suites_to_run = [s for s in VALIDATION_SUITES if not args.fast or s["fast"]]

    start_all_time = time.time()
    results = []

    for idx, suite in enumerate(suites_to_run, 1):
        suite_abs_path = os.path.join(validation_dir, suite["script"])
        suite_dir = os.path.dirname(suite_abs_path)
        suite_name = suite["name"]

        if not os.path.exists(suite_abs_path):
            print(f"[{idx}/{len(suites_to_run)}] ERROR: Suite script missing at {suite_abs_path}")
            results.append((suite_name, suite["title"], "MISSING", 0.0))
            continue

        print("-" * 80)
        print(f"[{idx}/{len(suites_to_run)}] Rung {idx}: {suite['title']}")
        print(f"       Description: {suite['description']}")
        print(f"       Running:     {suite['script']}")
        print("-" * 80)

        cmd = [sys.executable, suite_abs_path]
        if args.clean:
            cmd.append("--clean")

        suite_start = time.time()
        res = subprocess.run(
            cmd,
            cwd=suite_dir,
            stdout=sys.stdout,
            stderr=sys.stderr,
            text=True
        )
        duration = time.time() - suite_start

        if res.returncode == 0:
            print(f"\n---> {suite_name}: PASSED ({duration:.2f}s)\n")
            results.append((suite_name, suite["title"], "PASSED", duration))
        else:
            print(f"\n---> {suite_name}: FAILED (exit code {res.returncode}, {duration:.2f}s)\n")
            results.append((suite_name, suite["title"], "FAILED", duration))

    total_duration = time.time() - start_all_time

    print("=" * 80)
    print("      VALIDATION LADDER SUMMARY REPORT")
    print("=" * 80)
    all_passed = True
    for suite_name, title, status, duration in results:
        status_str = f"PASSED ({duration:.2f}s)" if status == "PASSED" else f"FAILED ({duration:.2f}s)"
        print(f"  * {title:<55} : {status_str}")
        if status != "PASSED":
            all_passed = False

    print("-" * 80)
    print(f"Total Validation Execution Time: {total_duration:.2f}s")

    passed_count = sum(1 for _, _, status, _ in results if status == "PASSED")
    total_count = len(results)
    if all_passed:
        print(f"\nAll {total_count} validation benchmarks passed.")
        sys.exit(0)
    else:
        print(f"\n{total_count - passed_count} of {total_count} validation benchmark(s) failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
