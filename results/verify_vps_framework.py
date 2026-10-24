#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification of Viscoelastic Phase Separation (VPS) Framework in CaCO3 SAXS
---------------------------------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

NOTE: This script is maintained for backward compatibility.
The architecture has been refactored into modular functional sub-packages:
- Analysis:  results/scripts/analysis/compute_vps_metrics.py
- Plotting:  results/scripts/plotting/plot_vps_framework.py
- Reporting: results/scripts/reporting/generate_vps_verification_report.py
- Narrative: results/docs/narratives/ (vps_theory.md, regime_descriptions.md, etc.)

Orchestrated via: results/scripts/run_collective_pipeline.py
"""

import os
import sys

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'analysis'))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'plotting'))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'reporting'))

from compute_vps_metrics import compute_vps_metrics, load_all_parameters, compute_parameter_summary_table
from plot_vps_framework import generate_vps_synthesis_plot
from generate_vps_verification_report import generate_report


def main():
    print("================================================================================")
    print("  VERIFYING VISCOELASTIC PHASE SEPARATION (VPS) SCATTERING FRAMEWORK")
    print("  (Delegating to modular analysis, plotting, and reporting engines)")
    print("================================================================================")
    
    # 1. Compute metrics & parameters synthesis
    metrics = compute_vps_metrics()
    
    # 2. Generate VPS Markdown verification report from standalone narratives
    generate_report()
    
    # 3. Generate synthesis figure
    generate_vps_synthesis_plot()
    
    print("================================================================================")
    print("  VPS FRAMEWORK VERIFICATION COMPLETED SUCCESSFULLY")
    print("================================================================================")


if __name__ == '__main__':
    main()
