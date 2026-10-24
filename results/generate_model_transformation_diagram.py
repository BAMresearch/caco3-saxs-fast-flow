#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model Transformation Diagram Generator
--------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

NOTE: This script is maintained for backward compatibility.
The architecture has been refactored into modular functional sub-packages:
- Plotting: results/scripts/plotting/generate_model_transformation_diagram.py

Orchestrated via: results/scripts/run_collective_pipeline.py --plot
"""

import os
import sys

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'plotting'))

from generate_model_transformation_diagram import main as generate_diagram


if __name__ == '__main__':
    generate_diagram()
