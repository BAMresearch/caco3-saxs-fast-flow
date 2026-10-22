#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Prepare Data and Composite Anchor for ND-i-005.
---------------------------------------------------------
1. Reads 50 Short-Distance and 50 Long-Distance raw .dat files.
2. Extracts actual experimental timestamps from series_h5_timescales.json.
3. Packages time-series data into wide-format CSV files:
   - data/dataNDi005_fitting_results_short_distance.csv
   - data/dataNDi005_fitting_results_long_distance.csv
   - data/frame_metadata.csv
4. Stitches the final stage datasets (00721 long + 00253 short) to produce:
   - data/dataNDi005_fitting_results_composite.csv
   - data/composite_final.dat
"""

import os
import sys
import json
import glob
import numpy as np
import pandas as pd

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WORKSPACE_ROOT = os.path.abspath(os.path.join(PROJECT_DIR, '..'))

META_JSON = os.path.join(WORKSPACE_ROOT, 'advanced_fitting_all_series', 'data', 'metadata', 'series_h5_timescales.json')

def prepare_dataset():
    print("=" * 80)
    print("  PREPARING ND-i-005 DATASET (40 mM CaCO3, pH 11.3, Cobalt)")
    print("=" * 80)

    # 1. Load Timestamps
    if os.path.exists(META_JSON):
        with open(META_JSON, 'r') as f:
            meta = json.load(f)
        ndi005_meta = next(item for item in meta if item.get('sample_id') == 'ND-i-005')
        timestamps = ndi005_meta['short_timestamps']
        print(f"Loaded {len(timestamps)} timestamps from metadata JSON (t: {timestamps[0]:.4f}s to {timestamps[-1]:.4f}s)")
    else:
        # Fallback if json not found
        timestamps = list(np.geomspace(0.06, 10.8769, 50))
        print("Warning: Used geometric fallback timestamps.")

    # 2. Short Distance
    s_files = sorted(glob.glob(os.path.join(PROJECT_DIR, 'data', 'raw_short_distance', '*.dat')))
    print(f"Found {len(s_files)} short-distance raw files.")
    assert len(s_files) == 50, f"Expected 50 short files, found {len(s_files)}"

    s_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']) for f in s_files]
    q_s = s_dfs[0]['q'].values

    short_dict = {'q': q_s}
    for f_idx, df in enumerate(s_dfs):
        short_dict[str(f_idx)] = df['I'].values
        short_dict['error_' + str(f_idx)] = df['dI'].values

    df_short = pd.DataFrame(short_dict)
    short_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_short_distance.csv')
    df_short.to_csv(short_csv, index=False)
    print(f"Saved short-distance wide CSV: {short_csv} ({df_short.shape})")

    # 3. Long Distance
    l_files = sorted(glob.glob(os.path.join(PROJECT_DIR, 'data', 'raw_long_distance', '*.dat')))
    print(f"Found {len(l_files)} long-distance raw files.")
    assert len(l_files) == 50, f"Expected 50 long files, found {len(l_files)}"

    l_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']) for f in l_files]
    q_l = l_dfs[0]['q'].values

    long_dict = {'q': q_l}
    for f_idx, df in enumerate(l_dfs):
        long_dict[str(f_idx)] = df['I'].values
        long_dict['error_' + str(f_idx)] = df['dI'].values

    df_long = pd.DataFrame(long_dict)
    long_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_long_distance.csv')
    df_long.to_csv(long_csv, index=False)
    print(f"Saved long-distance wide CSV: {long_csv} ({df_long.shape})")

    # 4. Frame Metadata
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

    # 5. Composite Final Stages
    final_s = os.path.join(PROJECT_DIR, 'data', 'raw_final_stages', 'ACC_eiger2_00253_sub.dat')
    final_l = os.path.join(PROJECT_DIR, 'data', 'raw_final_stages', 'ACC_eiger2_00721_sub.dat')
    if os.path.exists(final_s) and os.path.exists(final_l):
        df_fs = pd.read_csv(final_s, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI'])
        df_fl = pd.read_csv(final_l, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI'])

        # Overlap analysis
        q_overlap_min, q_overlap_max = 0.0225, 0.20
        sel_fs = (df_fs['q'] >= q_overlap_min) & (df_fs['q'] <= q_overlap_max)
        sel_fl = (df_fl['q'] >= q_overlap_min) & (df_fl['q'] <= q_overlap_max)
        
        q_common = df_fs.loc[sel_fs, 'q'].values
        I_fs = df_fs.loc[sel_fs, 'I'].values
        I_fl_interp = np.interp(q_common, df_fl['q'].values, df_fl['I'].values)
        k_scale = float(np.median(I_fs / I_fl_interp))
        print(f"Final stage scale factor (Short / Long): k_scale = {k_scale:.4f}")

        # Construct composite scattering curve:
        # Long distance scaled by k_scale up to q = 0.04 nm^-1
        # Short distance from q >= 0.04 nm^-1
        sel_long_part = df_fl['q'] < 0.04
        q_comp_long = df_fl.loc[sel_long_part, 'q'].values
        I_comp_long = df_fl.loc[sel_long_part, 'I'].values * k_scale
        dI_comp_long = df_fl.loc[sel_long_part, 'dI'].values * k_scale

        sel_short_part = df_fs['q'] >= 0.04
        q_comp_short = df_fs.loc[sel_short_part, 'q'].values
        I_comp_short = df_fs.loc[sel_short_part, 'I'].values
        dI_comp_short = df_fs.loc[sel_short_part, 'dI'].values

        q_comp = np.concatenate([q_comp_long, q_comp_short])
        I_comp = np.concatenate([I_comp_long, I_comp_short])
        dI_comp = np.concatenate([dI_comp_long, dI_comp_short])

        sort_idx = np.argsort(q_comp)
        df_composite = pd.DataFrame({
            'q': q_comp[sort_idx],
            'I': I_comp[sort_idx],
            'dI': dI_comp[sort_idx]
        })
        comp_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_composite.csv')
        df_composite.to_csv(comp_csv, index=False)

        comp_dat = os.path.join(PROJECT_DIR, 'data', 'composite_final.dat')
        with open(comp_dat, 'w') as f_out:
            f_out.write("# q [nm^-1]\tI(q) [cm^-1]\terror [cm^-1]\n")
            for _, r in df_composite.iterrows():
                f_out.write(f"{r['q']:.8e}\t{r['I']:.8e}\t{r['dI']:.8e}\n")
        print(f"Saved composite scattering curve: {comp_csv} & {comp_dat} ({len(df_composite)} points, q: [{df_composite['q'].min():.5f}, {df_composite['q'].max():.5f}] nm^-1)")

    print("=" * 80)
    print("  DATA PREPARATION COMPLETED SUCCESSFULLY!")
    print("=" * 80)

if __name__ == '__main__':
    prepare_dataset()
