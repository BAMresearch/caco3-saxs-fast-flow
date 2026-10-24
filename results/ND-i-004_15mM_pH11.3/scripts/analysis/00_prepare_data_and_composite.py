#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Prepare Data and Composite Anchor for ND-i-004.
---------------------------------------------------------
Dataset: 15.0 mM CaCO3, pH 11.3 (Transitional Supersaturation Regime)
1. Confirms 50 Short-Distance (00204), 50 Long-Distance (00625), and final stage raw .dat files.
2. Extracts actual experimental timestamps from series_h5_timescales.json.
3. Packages time-series data into wide-format CSV files:
   - data/dataNDi004_fitting_results_short_distance.csv
   - data/dataNDi004_fitting_results_long_distance.csv
   - data/frame_metadata.csv
4. Stitches the final stage datasets to produce:
   - data/dataNDi004_fitting_results_composite.csv
   - data/composite_final.dat
"""

import os
import sys
import json
import glob
import shutil
import numpy as np
import pandas as pd

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
WORKSPACE_ROOT = os.path.abspath(os.path.join(PROJECT_DIR, '..'))

SRC_DIR = os.path.join(WORKSPACE_ROOT, 'visualisations', '15 mM pH 11.3 ND-i-004 selected')
META_JSON = os.path.join(WORKSPACE_ROOT, 'advanced_fitting_all_series', 'data', 'metadata', 'series_h5_timescales.json')


def prepare_dataset():
    print("=" * 80)
    print("  PREPARING ND-i-004 DATASET (15 mM CaCO3, pH 11.3, Pure)")
    print("=" * 80)

    # 1. Setup local raw directories
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

    # 2. Load Timestamps
    if os.path.exists(META_JSON):
        with open(META_JSON, 'r') as f:
            meta = json.load(f)
        ndi004_meta = next(item for item in meta if item.get('sample_id') == 'ND-i-004')
        timestamps = ndi004_meta['short_timestamps']
        print(f"Loaded {len(timestamps)} timestamps from metadata JSON (t: {timestamps[0]:.4f}s to {timestamps[-1]:.4f}s)")
    else:
        timestamps = list(np.geomspace(0.06, 10.8769, 50))
        print("Warning: Used geometric fallback timestamps.")

    # 3. Short Distance Wide CSV
    s_files = sorted(glob.glob(os.path.join(raw_s_dir, '*.dat')))
    s_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']) for f in s_files]
    q_s = s_dfs[0]['q'].values

    short_dict = {'q': q_s}
    for f_idx, df in enumerate(s_dfs):
        short_dict[str(f_idx)] = df['I'].values
        short_dict['error_' + str(f_idx)] = df['dI'].values

    df_short = pd.DataFrame(short_dict)
    short_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_short_distance.csv')
    df_short.to_csv(short_csv, index=False)
    print(f"Saved short-distance wide CSV: {short_csv} ({df_short.shape})")

    # 4. Long Distance Wide CSV
    l_files = sorted(glob.glob(os.path.join(raw_l_dir, '*.dat')))
    l_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']) for f in l_files]
    q_l = l_dfs[0]['q'].values

    long_dict = {'q': q_l}
    for f_idx, df in enumerate(l_dfs):
        long_dict[str(f_idx)] = df['I'].values
        long_dict['error_' + str(f_idx)] = df['dI'].values

    df_long = pd.DataFrame(long_dict)
    long_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_long_distance.csv')
    df_long.to_csv(long_csv, index=False)
    print(f"Saved long-distance wide CSV: {long_csv} ({df_long.shape})")

    # 5. Frame Metadata
    meta_records = []
    for f_idx in range(50):
        meta_records.append({
            'frame_index': f_idx,
            'time_seconds': float(timestamps[f_idx]),
            'short_file': os.path.basename(s_files[f_idx]),
            'long_file': os.path.basename(l_files[f_idx])
        })
    df_meta = pd.DataFrame(meta_records)
    meta_csv = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
    df_meta.to_csv(meta_csv, index=False)
    print(f"Saved frame metadata table: {meta_csv}")

    # 6. Composite Stitched SAXS Curve (Long Distance 00627 + Short Distance 00221)
    final_s = os.path.join(raw_f_dir, 'ACC_eiger2_00221_sub.dat')
    final_l = os.path.join(raw_f_dir, 'ACC_eiger2_00627_sub.dat')
    if not os.path.exists(final_s):
        final_s = os.path.join(raw_f_dir, 'ACC_eiger2_00205_sub.dat')
    if not os.path.exists(final_l):
        final_l = l_files[-1]

    df_fs = pd.read_csv(final_s, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna()
    df_fl = pd.read_csv(final_l, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna()

    # Overlap analysis across q in [0.03, 0.06] nm^-1
    q_overlap_min, q_overlap_max = 0.03, 0.06
    sel_fs = (df_fs['q'] >= q_overlap_min) & (df_fs['q'] <= q_overlap_max)
    sel_fl = (df_fl['q'] >= q_overlap_min) & (df_fl['q'] <= q_overlap_max)

    q_common = df_fs.loc[sel_fs, 'q'].values
    I_fs = df_fs.loc[sel_fs, 'I'].values
    I_fl_interp = np.interp(q_common, df_fl['q'].values, df_fl['I'].values)
    k_scale = float(np.median(I_fl_interp / np.maximum(I_fs, 1e-4)))
    print(f"Final stage scale factor (Long / Short overlap): k_scale = {k_scale:.4f}")

    # Construct composite scattering curve:
    # Long distance up to q = 0.045 nm^-1
    # Short distance from q >= 0.045 nm^-1 (scaled by k_scale ~ 1.135 for exact seam continuity)
    q_stitch = 0.045
    sel_long_part = df_fl['q'] < q_stitch
    q_comp_long = df_fl.loc[sel_long_part, 'q'].values
    I_comp_long = df_fl.loc[sel_long_part, 'I'].values
    dI_comp_long = df_fl.loc[sel_long_part, 'dI'].values

    sel_short_part = df_fs['q'] >= q_stitch
    q_comp_short = df_fs.loc[sel_short_part, 'q'].values
    I_comp_short = df_fs.loc[sel_short_part, 'I'].values * k_scale
    dI_comp_short = df_fs.loc[sel_short_part, 'dI'].values * k_scale

    q_comp = np.concatenate([q_comp_long, q_comp_short])
    I_comp = np.concatenate([I_comp_long, I_comp_short])
    dI_comp = np.concatenate([dI_comp_long, dI_comp_short])

    sort_idx = np.argsort(q_comp)
    df_composite = pd.DataFrame({
        'q': q_comp[sort_idx],
        'I': I_comp[sort_idx],
        'dI': dI_comp[sort_idx]
    })
    comp_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_composite.csv')
    df_composite.to_csv(comp_csv, index=False)

    comp_dat = os.path.join(PROJECT_DIR, 'data', 'composite_final.dat')
    with open(comp_dat, 'w') as f_out:
        f_out.write("# q [nm^-1]\tI(q) [cm^-1]\terror [cm^-1]\n")
        for _, r in df_composite.iterrows():
            f_out.write(f"{r['q']:.8e}\t{r['I']:.8e}\t{r['dI']:.8e}\n")
    print(f"Saved composite scattering curve: {comp_csv} & {comp_dat} ({len(df_composite)} points, q: [{df_composite['q'].min():.5f}, {df_composite['q'].max():.5f}] nm^-1)")

    print("=" * 80)
    print("  DATA PREPARATION FOR ND-i-004 COMPLETED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == '__main__':
    prepare_dataset()
