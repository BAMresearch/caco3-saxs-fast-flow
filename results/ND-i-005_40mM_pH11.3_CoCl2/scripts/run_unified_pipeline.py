#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Master Pipeline Runner for Unified ND-i-005 SAXS Analysis Suite.
----------------------------------------------------------------
Orchestrates end-to-end data preparation, baseline composite fitting,
time-resolved short/long series modelling, diagnostic frame rendering,
cross-scale trajectory synthesis, and interactive dashboard compilation.

Usage:
  /home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py [--all] [--step STEP]
"""

import os
import sys
import time
import argparse
import subprocess

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PYTHON_BIN = '/home/tomek/anaconda3/bin/python' if os.path.exists('/home/tomek/anaconda3/bin/python') else sys.executable

STEPS = [
    ("Step 0: Prepare Data and Stitched Composite Anchor", os.path.join(PROJECT_DIR, 'scripts', 'analysis', '00_prepare_data_and_composite.py')),
    ("Step 1: Fit Multiscale Composite Baseline", os.path.join(PROJECT_DIR, 'scripts', 'analysis', '00_fit_composite_baseline.py')),
    ("Step 2: Fit Short-Distance Time Series (3.0 m, 50 frames)", os.path.join(PROJECT_DIR, 'scripts', 'analysis', '01_fit_short_distance_models.py')),
    ("Step 3: Fit Long-Distance Time Series (31.0 m, 50 frames)", os.path.join(PROJECT_DIR, 'scripts', 'analysis', '02_fit_long_distance_models.py')),
    ("Step 4: Generate All Universal Diagnostic Frame Plots (100 plots)", os.path.join(PROJECT_DIR, 'scripts', 'plotting', '03_generate_all_frame_plots.py')),
    ("Step 5: Cross-Scale Trajectory Synthesis", os.path.join(PROJECT_DIR, 'scripts', 'plotting', '04_compare_parameter_trajectories.py')),
    ("Step 6: Build Master Interactive Dashboard", os.path.join(PROJECT_DIR, 'scripts', 'dashboard', '05_build_unified_dashboard.py'))
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
    parser = argparse.ArgumentParser(description="Unified ND-i-005 SAXS Pipeline Runner")
    parser.add_argument('--all', action='store_true', help="Run all pipeline steps end-to-end")
    parser.add_argument('--fit', '--analysis', dest='fit', action='store_true', help="Run data preparation and fitting models")
    parser.add_argument('--plot', action='store_true', help="Run all diagnostic frame and trajectory plots")
    parser.add_argument('--dashboard', action='store_true', help="Compile unified interactive HTML dashboard")
    parser.add_argument('--step', type=int, choices=range(len(STEPS)), help="Run a specific step (0 to 6)")
    args = parser.parse_args()

    print("*" * 80)
    print("  UNIFIED FITTING SUITE FOR DATASET ND-i-005 (CaCO3 40 mM, pH 11.3, Cobalt)")
    print(f"  Python runtime: {PYTHON_BIN}")
    print(f"  Project root:   {PROJECT_DIR}")
    print("*" * 80)

    start_total = time.time()

    if args.step is not None:
        desc, script = STEPS[args.step]
        run_command(desc, script)
    elif args.all:
        for desc, script in STEPS:
            run_command(desc, script)
    else:
        # Default: verify synthesis and dashboard compilation
        for desc, script in STEPS[4:]:
            run_command(desc, script)

    total_time = time.time() - start_total
    print("\n" + "*" * 80)
    print(f"  PIPELINE FINISHED IN {total_time:.1f} SECONDS!")
    print(f"  Master Dashboard: file://{os.path.join(PROJECT_DIR, 'output', 'dashboard', 'index.html')}")
    print("*" * 80)


if __name__ == '__main__':
    main()
