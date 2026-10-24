#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collective Pipeline Orchestrator for Multi-Scale SAXS Synthesis
--------------------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Orchestrates collective analysis, plotting, reporting, and dashboard compilation
across all experimental series (ND-i-001 through ND-i-009).

Allows independent execution of each functional responsibility:
  python run_collective_pipeline.py --fit       # (or --analysis) Run data synthesis and metrics
  python run_collective_pipeline.py --plot      # Run all collective cross-series plots and diagrams
  python run_collective_pipeline.py --report    # Build VPS Markdown verification report from narratives
  python run_collective_pipeline.py --dashboard # Build master interactive HTML dashboard
  python run_collective_pipeline.py --all       # Run end-to-end collective pipeline
"""

import os
import sys
import time
import argparse
import subprocess

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
PYTHON_BIN = '/home/tomek/anaconda3/bin/python' if os.path.exists('/home/tomek/anaconda3/bin/python') else sys.executable

STEPS = {
    'analysis': [
        ("Compute VPS Metrics & Parameters Synthesis", os.path.join(SCRIPT_DIR, 'analysis', 'compute_vps_metrics.py')),
        ("Extract Cross-Series Trajectories", os.path.join(SCRIPT_DIR, 'analysis', 'extract_cross_series_trajectories.py'))
    ],
    'plotting': [
        ("Plot VPS Synthesis 4-Panel Figure", os.path.join(SCRIPT_DIR, 'plotting', 'plot_vps_framework.py')),
        ("Plot Cross-Series Trajectories Comparison", os.path.join(SCRIPT_DIR, 'plotting', 'plot_cross_series_trajectories.py')),
        ("Generate Model Transformation Diagram (SVG)", os.path.join(SCRIPT_DIR, 'plotting', 'generate_model_transformation_diagram.py')),
        ("Plot Standardized Composite Baseline Fits", os.path.join(SCRIPT_DIR, 'plotting', 'plot_standardized_composite_baselines.py'))
    ],
    'reporting': [
        ("Generate VPS Verification Report from Standalone Narratives", os.path.join(SCRIPT_DIR, 'reporting', 'generate_vps_verification_report.py'))
    ],
    'dashboard': [
        ("Build Master Synthesis Dashboard (HTML)", os.path.join(SCRIPT_DIR, 'dashboard', 'build_master_dashboard.py'))
    ]
}


def run_command(desc, script_path):
    print("\n" + "=" * 80)
    print(f"  RUNNING: {desc}")
    print(f"  Script:  {script_path}")
    print("=" * 80)
    start_t = time.time()
    res = subprocess.run([PYTHON_BIN, script_path], cwd=RESULTS_DIR)
    elapsed = time.time() - start_t
    if res.returncode != 0:
        print(f"\n[ERROR] {desc} failed with exit code {res.returncode}")
        sys.exit(res.returncode)
    print(f"[SUCCESS] {desc} completed in {elapsed:.1f} s.")


def main():
    parser = argparse.ArgumentParser(description="Collective Multi-Scale SAXS Synthesis Pipeline Orchestrator")
    parser.add_argument('--all', action='store_true', help="Execute entire collective pipeline (analysis, plotting, reporting, dashboard)")
    parser.add_argument('--fit', '--analysis', dest='analysis', action='store_true', help="Execute data analysis and trajectory extraction")
    parser.add_argument('--plot', action='store_true', help="Execute all collective visualization and plotting scripts")
    parser.add_argument('--report', action='store_true', help="Compile Markdown verification reports from standalone narrative texts")
    parser.add_argument('--dashboard', action='store_true', help="Compile interactive master HTML synthesis dashboard")
    args = parser.parse_args()

    # If no specific flags passed, show help or run all
    if not (args.all or args.analysis or args.plot or args.report or args.dashboard):
        parser.print_help()
        print("\n[Notice] No execution mode specified. Use --all to run everything, or select --fit, --plot, --report, or --dashboard.")
        return

    print("*" * 80)
    print("  COLLECTIVE SAXS SYNTHESIS PIPELINE")
    print(f"  Python runtime: {PYTHON_BIN}")
    print(f"  Results root:   {RESULTS_DIR}")
    print("*" * 80)

    start_total = time.time()

    if args.all or args.analysis:
        for desc, script in STEPS['analysis']:
            run_command(desc, script)

    if args.all or args.plot:
        for desc, script in STEPS['plotting']:
            run_command(desc, script)

    if args.all or args.report:
        for desc, script in STEPS['reporting']:
            run_command(desc, script)

    if args.all or args.dashboard:
        for desc, script in STEPS['dashboard']:
            run_command(desc, script)

    total_time = time.time() - start_total
    print("\n" + "*" * 80)
    print(f"  COLLECTIVE PIPELINE COMPLETED IN {total_time:.1f} SECONDS!")
    print(f"  Verification Report: file://{os.path.join(RESULTS_DIR, 'output', 'vps_framework_verification_report.md')}")
    print(f"  Master Dashboard:    file://{os.path.join(RESULTS_DIR, 'index.html')}")
    print("*" * 80)


if __name__ == '__main__':
    main()
