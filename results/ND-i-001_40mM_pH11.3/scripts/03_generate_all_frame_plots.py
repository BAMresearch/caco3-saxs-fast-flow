#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Backward-compatible wrapper for 03_generate_all_frame_plots.py.
Delegates execution to scripts/plotting/03_generate_all_frame_plots.py.
"""
import os
import sys
import subprocess

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
TARGET = os.path.join(SCRIPT_DIR, 'plotting', '03_generate_all_frame_plots.py')
PYTHON_BIN = '/home/tomek/anaconda3/bin/python' if os.path.exists('/home/tomek/anaconda3/bin/python') else sys.executable

if __name__ == '__main__':
    res = subprocess.run([PYTHON_BIN, TARGET] + sys.argv[1:], cwd=os.path.dirname(SCRIPT_DIR))
    sys.exit(res.returncode)
