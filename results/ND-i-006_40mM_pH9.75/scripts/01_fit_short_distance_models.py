#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Backward-compatible wrapper for 01_fit_short_distance_models.py.
Delegates execution to scripts/analysis/01_fit_short_distance_models.py.
"""
import os
import sys
import subprocess

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
TARGET = os.path.join(SCRIPT_DIR, 'analysis', '01_fit_short_distance_models.py')
PYTHON_BIN = '/home/tomek/anaconda3/bin/python' if os.path.exists('/home/tomek/anaconda3/bin/python') else sys.executable

if __name__ == '__main__':
    res = subprocess.run([PYTHON_BIN, TARGET] + sys.argv[1:], cwd=os.path.dirname(SCRIPT_DIR))
    sys.exit(res.returncode)
