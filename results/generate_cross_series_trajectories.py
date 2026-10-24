#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-Series Trajectory Synthesis and Comparison
------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

NOTE: This script is maintained for backward compatibility.
The architecture has been refactored into modular functional sub-packages:
- Analysis: results/scripts/analysis/extract_cross_series_trajectories.py
- Plotting: results/scripts/plotting/plot_cross_series_trajectories.py

Orchestrated via: results/scripts/run_collective_pipeline.py
"""

import os
import sys

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'analysis'))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'plotting'))

from extract_cross_series_trajectories import load_series_trajectories, export_trajectories_json
from plot_cross_series_trajectories import generate_synthesis_figure


def main():
    print("================================================================================")
    print("  CROSS-SERIES TRAJECTORY SYNTHESIS")
    print("  (Delegating to modular analysis and plotting engines)")
    print("================================================================================")
    trajectories = load_series_trajectories()
    export_trajectories_json(trajectories)
    generate_synthesis_figure(trajectories)
    print("================================================================================")


if __name__ == '__main__':
    main()
