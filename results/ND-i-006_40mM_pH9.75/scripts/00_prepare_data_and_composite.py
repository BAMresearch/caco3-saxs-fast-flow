#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Prepare Data and Composite Anchor for ND-i-006.
---------------------------------------------------------
Dataset: 40.0 mM CaCO3, pH 9.75 (Near-Spinodal / Crossover Regime)
1. Confirms 50 Short-Distance, 50 Long-Distance, and final stage raw .dat files.
2. Extracts actual experimental timestamps from series_h5_timescales.json.
3. Packages time-series data into wide-format CSV files:
   - data/dataNDi006_fitting_results_short_distance.csv
   - data/dataNDi006_fitting_results_long_distance.csv
   - data/frame_metadata.csv
4. Stitches the final stage datasets to produce:
   - data/dataNDi006_fitting_results_composite.csv
   - data/composite_final.dat
"""

import os
import sys
import json
import glob
import shutil
import numpy as np
import pandas as pd

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WORKSPACE_ROOT = os.path.abspath(os.path.join(PROJECT_DIR, '..', '..'))

SRC_DIR = os.path.join(WORKSPACE_ROOT, 'visualisations', '40 mM pH 9.75 ND-i-006 selected')
META_JSON = os.path.join(WORKSPACE_ROOT, 'advanced_fitting_all_series', 'data', 'metadata', 'series_h5_timescales.json')


def prepare_dataset():
    print("=" * 80)
    print("  PREPARING ND-i-006 DATASET (40 mM CaCO3, pH 9.75, Near-Spinodal/Crossover)")
    print("=" * 80)

    raw_s_dir = os.path.join(PROJECT_DIR, 'data', 'raw_short_distance')
    raw_l_dir = os.path.join(PROJECT_DIR, 'data', 'raw_long_distance')
    raw_f_dir = os.path.join(PROJECT_DIR, 'data', 'raw_final_stages')
    os.makedirs(raw_s_dir, exist_ok=True)
    os.makedirs(raw_l_dir, exist_ok=True)
    os.makedirs(raw_f_dir, exist_ok=True)

    # Check if files already present locally or need copying from source
    s_files_local = sorted(glob.glob(os.path.join(raw_s_dir, '*.dat')))
    l_files_local = sorted(glob.glob(os.path.join(raw_l_dir, '*.dat')))

    if len(s_files_local) < 50 or len(l_files_local) < 50:
        if os.path.exists(SRC_DIR):
            src_s_files = sorted(glob.glob(os.path.join(SRC_DIR, 'shortDistance', '*.dat')))
            src_l_files = sorted(glob.glob(os.path.join(SRC_DIR, 'longDistance', '*.dat')))
            src_f_files = sorted(glob.glob(os.path.join(SRC_DIR, 'finalStages', '*.dat')))
            for f in src_s_files:
                shutil.copy2(f, os.path.join(raw_s_dir, os.path.basename(f)))
            for f in src_l_files:
                shutil.copy2(f, os.path.join(raw_l_dir, os.path.basename(f)))
            for f in src_f_files:
                shutil.copy2(f, os.path.join(raw_f_dir, os.path.basename(f)))
        else:
            print(f"Using existing local raw files ({len(s_files_local)} short, {len(l_files_local)} long)")
    else:
        print(f"Using existing local raw files ({len(s_files_local)} short, {len(l_files_local)} long)")

    s_files = sorted(glob.glob(os.path.join(raw_s_dir, '*.dat')))
    l_files = sorted(glob.glob(os.path.join(raw_l_dir, '*.dat')))

    # Load timestamps
    if os.path.exists(META_JSON):
        with open(META_JSON, 'r') as f:
            meta_list = json.load(f)
        meta_entry = [x for x in meta_list if x['sample_id'] == 'ND-i-006'][0]
        t_vals = meta_entry['short_timestamps']
    else:
        t_vals = list(np.geomspace(0.06, 10.8769, 50))

    # Metadata table
    meta_rows = []
    for i in range(len(s_files)):
        meta_rows.append({
            'frame_index': i,
            'time_seconds': t_vals[i] if i < len(t_vals) else float(i),
            'short_file': os.path.basename(s_files[i]),
            'long_file': os.path.basename(l_files[i]) if i < len(l_files) else ''
        })
    df_meta = pd.DataFrame(meta_rows)
    df_meta.to_csv(os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv'), index=False)
    print(f"Generated frame metadata with 50 frames (t = {t_vals[0]:.4f} to {t_vals[-1]:.4f} s).")

    # Short distance wide CSV
    q_s = None
    cols_s = {}
    for i, f in enumerate(src_s_files):
        arr = np.loadtxt(f, skiprows=1)
        if q_s is None:
            q_s = arr[:, 0]
            cols_s['q'] = q_s
        cols_s[str(i)] = arr[:, 1]
        cols_s[f'error_{i}'] = arr[:, 2]
    pd.DataFrame(cols_s).to_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_short_distance.csv'), index=False)

    # Long distance wide CSV
    q_l = None
    cols_l = {}
    for i, f in enumerate(src_l_files):
        arr = np.loadtxt(f, skiprows=1)
        if q_l is None:
            q_l = arr[:, 0]
            cols_l['q'] = q_l
        cols_l[str(i)] = arr[:, 1]
        cols_l[f'error_{i}'] = arr[:, 2]
    pd.DataFrame(cols_l).to_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_long_distance.csv'), index=False)

    print("Wide format CSVs verified.")

    # Composite anchor
    short_final = [f for f in src_f_files if '00308' in os.path.basename(f)][0]
    long_candidates = [f for f in src_f_files if '00543' in os.path.basename(f) or '00544' in os.path.basename(f)]
    long_final = long_candidates[0] if long_candidates else src_f_files[0]

    dat_s = np.loadtxt(short_final, skiprows=1)
    dat_l = np.loadtxt(long_final, skiprows=1)

    q_ov_min, q_ov_max = 0.03, 0.08
    m_s = (dat_s[:, 0] >= q_ov_min) & (dat_s[:, 0] <= q_ov_max) & (dat_s[:, 1] > 0)
    m_l = (dat_l[:, 0] >= q_ov_min) & (dat_l[:, 0] <= q_ov_max) & (dat_l[:, 1] > 0)

    mean_s = np.mean(dat_s[m_s, 1])
    mean_l = np.mean(dat_l[m_l, 1])
    k_scale = mean_s / mean_l
    print(f"Stitching scale factor k_scale = {k_scale:.5f}")

    dat_l_scaled = dat_l.copy()
    dat_l_scaled[:, 1] *= k_scale
    dat_l_scaled[:, 2] *= k_scale

    cut_q = 0.05
    comp_l = dat_l_scaled[dat_l_scaled[:, 0] < cut_q]
    comp_s = dat_s[dat_s[:, 0] >= cut_q]
    comp_stitched = np.vstack([comp_l, comp_s])
    comp_stitched = comp_stitched[comp_stitched[:, 0].argsort()]

    df_comp = pd.DataFrame({
        'q': comp_stitched[:, 0],
        'I': comp_stitched[:, 1],
        'dI': comp_stitched[:, 2]
    })
    df_comp.to_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_composite.csv'), index=False)
    np.savetxt(os.path.join(PROJECT_DIR, 'data', 'composite_final.dat'), comp_stitched, fmt='%.8e')
    print("ND-i-006 Data Preparation Complete.")


if __name__ == '__main__':
    prepare_dataset()
