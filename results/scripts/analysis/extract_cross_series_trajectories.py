#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-Series Trajectory Extraction and Normalization Engine
-----------------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Extracts, aligns, and normalises time-resolved structural trajectories
across all calcium carbonate reaction series (ND-i-001 to ND-i-009).

Pure computation module:
- No matplotlib plotting logic.
- Standalone execution or modular import.

Output:
  results/output/cross_series_trajectories.json
"""

import os
import sys
import json
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# Series metadata and theoretical maximum volume fractions
SERIES_CONFIG = {
    'ND-i-003': {
        'name': 'ND-i-003 (10 mM, pure)',
        'dir': 'ND-i-003_10mM_pH11.3',
        'conc_mM': 10.0,
        'pH': 11.3,
        'phi_max': 5.0e-4,
        'colour': '#2563eb',
        'marker': 's',
        'regime': 'Sub-spinodal nucleation'
    },
    'ND-i-004': {
        'name': 'ND-i-004 (15 mM, pure)',
        'dir': 'ND-i-004_15mM_pH11.3',
        'conc_mM': 15.0,
        'pH': 11.3,
        'phi_max': 7.5e-4,
        'colour': '#0891b2',
        'marker': '^',
        'regime': 'Crossover boundary'
    },
    'ND-i-002': {
        'name': 'ND-i-002 (20 mM, pure)',
        'dir': 'ND-i-002_20mM_pH11.3',
        'conc_mM': 20.0,
        'pH': 11.3,
        'phi_max': 1.0e-3,
        'colour': '#16a34a',
        'marker': 'D',
        'regime': 'Spinodal gelation'
    },
    'ND-i-009': {
        'name': 'ND-i-009 (40 mM, pure benchmark)',
        'dir': 'ND-i-009_40mM_pH10.0',
        'conc_mM': 40.0,
        'pH': 10.0,
        'phi_max': 2.0e-3,
        'colour': '#ea580c',
        'marker': 'o',
        'regime': 'Primary pure 40 mM benchmark'
    },
    'ND-i-001': {
        'name': 'ND-i-001 (40 mM, pure, pH 11.3)',
        'dir': 'ND-i-001_40mM_pH11.3',
        'conc_mM': 40.0,
        'pH': 11.3,
        'phi_max': 2.0e-3,
        'colour': '#dc2626',
        'marker': 'p',
        'regime': 'Rapid quench spinodal'
    },
    'ND-i-005': {
        'name': 'ND-i-005 (40 mM, +Co2+ additive)',
        'dir': 'ND-i-005_40mM_pH11.3_CoCl2',
        'conc_mM': 40.0,
        'pH': 11.3,
        'phi_max': 2.0e-3,
        'colour': '#9333ea',
        'marker': 'v',
        'regime': 'Dopant additive branch'
    }
}


def load_series_trajectories():
    """Load and process time series parameters for all series."""
    trajectories = {}

    for sid, cfg in SERIES_CONFIG.items():
        p_ld = os.path.join(RESULTS_DIR, cfg['dir'], 'output', 'long_distance', 'parameters.csv')
        p_sd = os.path.join(RESULTS_DIR, cfg['dir'], 'output', 'short_distance', 'parameters.csv')

        df_ld = pd.read_csv(p_ld) if os.path.exists(p_ld) else pd.DataFrame()
        df_sd = pd.read_csv(p_sd) if os.path.exists(p_sd) else pd.DataFrame()

        t_col = 'time_seconds' if 'time_seconds' in df_ld.columns else 'time'
        t_vals = df_ld[t_col].values if not df_ld.empty else np.array([])

        phi_vals = np.zeros_like(t_vals)
        xi_fluc = np.nan
        if sid == 'ND-i-003':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['Rmean'].values
            r_meso_label = r'ACC radius $\langle R_{\rm ACC} \rangle$'
            r_prim = np.nan
            xi_fluc = np.nan
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 2167.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 5.1e-5)
        elif sid == 'ND-i-004':
            phi_vals = (df_ld['phi_net'].fillna(0.0) + df_ld['phi_acc'].fillna(0.0)).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 4.94
            xi_fluc = np.nan
            r_drop = np.full_like(t_vals, 2800.0)
            phi_drop = np.full_like(t_vals, 4.5e-5)
        elif sid == 'ND-i-002':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 11.20
            xi_fluc = 7.63
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 7163.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 8.4e-5)
        elif sid == 'ND-i-009':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 8.35
            xi_fluc = 2.72
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 2904.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 5.9e-4)
        elif sid == 'ND-i-001':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 6.24
            xi_fluc = 3.05
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 3015.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 4.98e-5)
        elif sid == 'ND-i-005':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values if 'xi' in df_ld.columns else df_ld['xi_fluc'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 7.03
            xi_fluc = 8.05
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 3530.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 2.92e-5)

        chi2_col = 'redchi' if 'redchi' in df_ld.columns else ('chi2_red' if 'chi2_red' in df_ld.columns else None)
        chi2_vals = df_ld[chi2_col].values if chi2_col else np.ones_like(t_vals)

        alpha_vals = phi_vals / cfg['phi_max'] if cfg['phi_max'] > 0 else np.zeros_like(phi_vals)
        alpha_vals = np.clip(alpha_vals, 0.0, 1.2)

        trajectories[sid] = {
            'config': cfg,
            'time': t_vals,
            'phi': phi_vals,
            'alpha': alpha_vals,
            'r_meso': r_meso,
            'r_meso_label': r_meso_label,
            'r_prim': r_prim,
            'xi_fluc': xi_fluc,
            'r_drop': r_drop,
            'phi_drop': phi_drop,
            'chi2': chi2_vals
        }

    return trajectories


def export_trajectories_json(trajectories=None):
    """Save serializable trajectory telemetry to JSON."""
    if trajectories is None:
        trajectories = load_series_trajectories()

    export_data = {}
    for sid, d in trajectories.items():
        cfg = d['config']
        t = d['time']
        step = 1 if len(t) <= 60 else (len(t) // 60)
        export_data[sid] = {
            'name': cfg['name'],
            'conc_mM': cfg['conc_mM'],
            'pH': cfg['pH'],
            'regime': cfg['regime'],
            'colour': cfg['colour'],
            'time': [round(float(v), 3) for v in t[::step]],
            'alpha': [round(float(v), 4) if np.isfinite(v) else None for v in d['alpha'][::step]],
            'phi': [float(f"{v:.3e}") if np.isfinite(v) else None for v in d['phi'][::step]],
            'r_meso': [round(float(v), 1) if np.isfinite(v) else None for v in d['r_meso'][::step]],
            'r_meso_label': d['r_meso_label'],
            'chi2': [round(float(v), 3) if np.isfinite(v) else None for v in d['chi2'][::step]]
        }

    json_path = os.path.join(OUT_DIR, 'cross_series_trajectories.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(export_data, f, indent=2)
    print(f"[Analysis] Saved cross-series trajectories JSON to: {json_path}")
    return export_data


def main():
    print("=" * 80)
    print("  RUNNING: Cross-Series Trajectory Extraction")
    print("=" * 80)
    trajectories = load_series_trajectories()
    export_trajectories_json(trajectories)
    print(f"[Success] Extracted trajectories for {len(trajectories)} series.")


if __name__ == '__main__':
    main()
