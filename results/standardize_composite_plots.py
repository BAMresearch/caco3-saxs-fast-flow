#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standardized Composite Baseline Plot Generator
----------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

NOTE: This script is maintained for backward compatibility.
The architecture has been refactored into modular functional sub-packages:
- Plotting: results/scripts/plotting/plot_standardized_composite_baselines.py

Orchestrated via: results/scripts/run_collective_pipeline.py --plot
"""

import os
import sys

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'plotting'))

from plot_standardized_composite_baselines import main as generate_standardized_plots


if __name__ == '__main__':
    generate_standardized_plots()
