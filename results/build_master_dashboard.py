#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master Multi-Scale SAXS Synthesis Dashboard Generator
-----------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

NOTE: This script is maintained for backward compatibility.
The architecture has been refactored into modular functional sub-packages:
- Dashboard Engine: results/scripts/dashboard/build_master_dashboard.py

Orchestrated via: results/scripts/run_collective_pipeline.py --dashboard
"""

import os
import sys

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'dashboard'))

from build_master_dashboard import build_master_dashboard


def main():
    print("================================================================================")
    print("  BUILDING MASTER MULTI-SCALE SAXS SYNTHESIS DASHBOARD")
    print("  (Delegating to results/scripts/dashboard/build_master_dashboard.py)")
    print("================================================================================")
    build_master_dashboard()
    print("================================================================================")


if __name__ == '__main__':
    main()
