#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Master Pipeline Runner for Unified ND-i-009 SAXS Analysis Suite.
----------------------------------------------------------------
Orchestrates end-to-end fitting, plotting, trajectory comparisons, and dashboard generation.

Usage:
  /home/tomek/anaconda3/bin/python unified_fitting_ND-i-009/scripts/run_unified_pipeline.py [--all] [--step STEP]
"""

import os
import sys
import time
import argparse
import subprocess

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PYTHON_BIN = sys.executable

STEPS = [
    ("Step 0: Fit Merged Composite Baseline", os.path.join(PROJECT_DIR, 'scripts', '00_fit_merged_baseline.py')),
    ("Step 1: Evolutionary Fit Short-Distance Models", os.path.join(PROJECT_DIR, 'scripts', '01_fit_short_distance_models.py')),
    ("Step 2: Fit Long-Distance Multi-Models", os.path.join(PROJECT_DIR, 'scripts', '02_fit_long_distance_models.py')),
    ("Step 3: Generate All Frame Diagnostic Plots", os.path.join(PROJECT_DIR, 'scripts', '03_generate_all_frame_plots.py')),
    ("Step 4: Parameter Trajectory Analysis & Comparisons", os.path.join(PROJECT_DIR, 'scripts', '04_compare_parameter_trajectories.py')),
    ("Step 5: Compile Master Interactive Dashboard", os.path.join(PROJECT_DIR, 'scripts', '05_build_unified_dashboard.py'))
]


def run_command(desc, script_path):
    print("\n" + "=" * 80)
    print(f"  RUNNING: {desc}")
    print(f"  Script: {script_path}")
    print("=" * 80)
    start_t = time.time()
    res = subprocess.run([PYTHON_BIN, script_path], cwd=PROJECT_DIR)
    elapsed = time.time() - start_t
    if res.returncode != 0:
        print(f"\n[ERROR] {desc} failed with exit code {res.returncode}")
        sys.exit(res.returncode)
    print(f"[SUCCESS] {desc} completed in {elapsed:.1f} s.")


def main():
    parser = argparse.ArgumentParser(description="Unified ND-i-009 SAXS Pipeline Runner")
    parser.add_argument('--all', action='store_true', help="Run all pipeline steps end-to-end")
    parser.add_argument('--step', type=int, choices=[0, 1, 2, 3, 4, 5], help="Run a specific step (0-5)")
    args = parser.parse_args()

    print("*" * 80)
    print("  UNIFIED FITTING SUITE FOR DATASET ND-i-009 (CaCO3 40 mM, pH 10.0)")
    print(f"  Python runtime: {PYTHON_BIN}")
    print(f"  Project root:   {PROJECT_DIR}")
    print("*" * 80)

    start_total = time.time()

    if args.step is not None:
        desc, script = STEPS[args.step]
        run_command(desc, script)
    else:
        # Default or --all: Run trajectory generation and dashboard compilation
        for desc, script in STEPS[3:]:  # steps 3, 4, 5
            run_command(desc, script)

    total_time = time.time() - start_total
    print("\n" + "*" * 80)
    print(f"  PIPELINE FINISHED IN {total_time:.1f} SECONDS!")
    print(f"  Master Dashboard: file://{os.path.join(PROJECT_DIR, 'output', 'dashboard', 'index.html')}")
    print("*" * 80)


if __name__ == '__main__':
    main()
