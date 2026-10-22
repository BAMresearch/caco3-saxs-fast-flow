#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Prepare Data and Composite Anchor for ND-i-003.
---------------------------------------------------------
1. Copies 50 Short-Distance, 50 Long-Distance, and final stage raw .dat files.
2. Extracts actual experimental timestamps from series_h5_timescales.json.
3. Packages time-series data into wide-format CSV files:
   - data/dataNDi003_fitting_results_short_distance.csv
   - data/dataNDi003_fitting_results_long_distance.csv
   - data/frame_metadata.csv
4. Stitches the final stage datasets (00661 long + 00177 short) to produce:
   - data/dataNDi003_fitting_results_composite.csv
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
WORKSPACE_ROOT = os.path.abspath(os.path.join(PROJECT_DIR, '..'))

SRC_DIR = os.path.join(WORKSPACE_ROOT, 'visualisations', '10 mM pH 11.3 ND-i-003 selected')
META_JSON = os.path.join(WORKSPACE_ROOT, 'advanced_fitting_all_series', 'data', 'metadata', 'series_h5_timescales.json')


def prepare_dataset():
    print("=" * 80)
    print("  PREPARING ND-i-003 DATASET (10 mM CaCO3, pH 11.3, Low Supersaturation)")
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
        ndi003_meta = next(item for item in meta if item.get('sample_id') == 'ND-i-003')
        timestamps = ndi003_meta['short_timestamps']
        print(f"Loaded {len(timestamps)} timestamps from metadata JSON (t: {timestamps[0]:.4f}s to {timestamps[-1]:.4f}s)")
    else:
        timestamps = list(np.geomspace(0.06, 10.8769, 50))
        print("Warning: Used geometric fallback timestamps.")

    # 3. Short Distance Wide CSV
    s_files = sorted(glob.glob(os.path.join(raw_s_dir, '*.dat')))
    s_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna() for f in s_files]
    q_s = s_dfs[0]['q'].values

    short_dict = {'q': q_s}
    for f_idx, df in enumerate(s_dfs):
        # Align to common q_s if minor point count differences exist
        if len(df) == len(q_s) and np.allclose(df['q'].values, q_s, atol=1e-5):
            short_dict[str(f_idx)] = df['I'].values
            short_dict['error_' + str(f_idx)] = df['dI'].values
        else:
            interp_I = np.interp(q_s, df['q'].values, df['I'].values)
            interp_dI = np.interp(q_s, df['q'].values, df['dI'].values)
            short_dict[str(f_idx)] = interp_I
            short_dict['error_' + str(f_idx)] = interp_dI

    df_short = pd.DataFrame(short_dict)
    short_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_short_distance.csv')
    df_short.to_csv(short_csv, index=False)
    print(f"Saved short-distance wide CSV: {short_csv} ({df_short.shape})")

    # 4. Long Distance Wide CSV
    l_files = sorted(glob.glob(os.path.join(raw_l_dir, '*.dat')))
    l_dfs = [pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna() for f in l_files]
    q_l = l_dfs[0]['q'].values

    long_dict = {'q': q_l}
    for f_idx, df in enumerate(l_dfs):
        if len(df) == len(q_l) and np.allclose(df['q'].values, q_l, atol=1e-5):
            long_dict[str(f_idx)] = df['I'].values
            long_dict['error_' + str(f_idx)] = df['dI'].values
        else:
            interp_I = np.interp(q_l, df['q'].values, df['I'].values)
            interp_dI = np.interp(q_l, df['q'].values, df['dI'].values)
            long_dict[str(f_idx)] = interp_I
            long_dict['error_' + str(f_idx)] = interp_dI

    df_long = pd.DataFrame(long_dict)
    long_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_long_distance.csv')
    df_long.to_csv(long_csv, index=False)
    print(f"Saved long-distance wide CSV: {long_csv} ({df_long.shape})")

    # 5. Metadata CSV
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
    print(f"Saved frame metadata CSV: {meta_csv} ({len(df_meta)} records)")

    # 6. Build Stitched Composite Curve from Final Stages
    # Long final: ACC_eiger2_00661_sub.dat; Short final: ACC_eiger2_00177_sub.dat
    long_final_path = os.path.join(raw_f_dir, 'ACC_eiger2_00661_sub.dat')
    short_final_path = os.path.join(raw_f_dir, 'ACC_eiger2_00177_sub.dat')

    df_lf = pd.read_csv(long_final_path, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna()
    df_sf = pd.read_csv(short_final_path, sep=r'\s+', skiprows=1, header=None, names=['q', 'I', 'dI']).dropna()

    # Seamless stitching threshold at q ~ 0.18 nm^-1 where long distance signal meets baseline noise floor
    q_stitch = 0.180
    part_long = df_lf[df_lf['q'] < q_stitch].copy()
    part_short = df_sf[df_sf['q'] >= q_stitch].copy()

    df_composite = pd.concat([part_long, part_short], ignore_index=True)
    df_composite = df_composite.sort_values('q').reset_index(drop=True)

    # Save ASCII composite
    comp_dat = os.path.join(PROJECT_DIR, 'data', 'composite_final.dat')
    with open(comp_dat, 'w') as f:
        f.write("# q [nm^-1]   I(q) [cm^-1]   dI [cm^-1]\n")
        for _, row in df_composite.iterrows():
            f.write(f"{row['q']:.8e}  {row['I']:.8e}  {row['dI']:.8e}\n")

    # Save wide composite CSV
    comp_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_composite.csv')
    df_comp_wide = pd.DataFrame({
        'q': df_composite['q'].values,
        '0': df_composite['I'].values,
        'error_0': df_composite['dI'].values
    })
    df_comp_wide.to_csv(comp_csv, index=False)
    print(f"Stitched composite scattering curve: {comp_dat} and {comp_csv} ({len(df_composite)} points, q: {df_composite['q'].min():.4e} to {df_composite['q'].max():.4e} nm^-1)")
    print("Data preparation complete!\n")


if __name__ == '__main__':
    prepare_dataset()
