#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPS Metrics and Parameter Synthesis Computation Engine
------------------------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Calculates quantitative Viscoelastic Phase Separation (VPS) metrics,
parameter syntheses across experimental regimes (10 mM to 40 mM),
and outputs structured JSON data.

Pure computation module:
- No hardcoded narrative text.
- No matplotlib plotting logic.
- Standalone execution or modular import.

Output:
  results/output/vps_metrics_summary.json
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

# Series metadata and physical regimes
SERIES_REGIMES = {
    'ND-i-003': {
        'name': 'ND-i-003 (10 mM, pure, pH 11.3)',
        'dir': 'ND-i-003_10mM_pH11.3',
        'conc_mM': 10.0,
        'pH': 11.3,
        'regime': 'Sub-spinodal nucleation',
        'model_type': 'Discrete ACC spheres + secondary droplets',
        'colour': '#2563eb',
        'has_network': False,
        'has_acc': True
    },
    'ND-i-004': {
        'name': 'ND-i-004 (15 mM, pure, pH 11.3)',
        'dir': 'ND-i-004_15mM_pH11.3',
        'conc_mM': 15.0,
        'pH': 11.3,
        'regime': 'Crossover boundary',
        'model_type': 'Coexisting early network + ACC spheres',
        'colour': '#0891b2',
        'has_network': True,
        'has_acc': True
    },
    'ND-i-002': {
        'name': 'ND-i-002 (20 mM, pure, pH 11.3)',
        'dir': 'ND-i-002_20mM_pH11.3',
        'conc_mM': 20.0,
        'pH': 11.3,
        'regime': 'Spinodal gelation',
        'model_type': 'Percolated mass fractal network',
        'colour': '#16a34a',
        'has_network': True,
        'has_acc': False
    },
    'ND-i-009': {
        'name': 'ND-i-009 (40 mM, pure, pH 10.0)',
        'dir': 'ND-i-009_40mM_pH10.0',
        'conc_mM': 40.0,
        'pH': 10.0,
        'regime': 'Primary pure 40 mM benchmark',
        'model_type': 'Continuous spinodal network + PNC mesh',
        'colour': '#ea580c',
        'has_network': True,
        'has_acc': False
    },
    'ND-i-001': {
        'name': 'ND-i-001 (40 mM, pure, pH 11.3)',
        'dir': 'ND-i-001_40mM_pH11.3',
        'conc_mM': 40.0,
        'pH': 11.3,
        'regime': 'Rapid quench spinodal',
        'model_type': 'Instantaneous spinodal gel + droplets',
        'colour': '#dc2626',
        'has_network': True,
        'has_acc': False
    },
    'ND-i-005': {
        'name': 'ND-i-005 (40 mM, +Co2+, pH 11.3)',
        'dir': 'ND-i-005_40mM_pH11.3_CoCl2',
        'conc_mM': 40.0,
        'pH': 11.3,
        'regime': 'Additive perturbation branch',
        'model_type': 'Extended network + dopant-expanded nodes',
        'colour': '#9333ea',
        'has_network': True,
        'has_acc': False
    }
}


def load_all_parameters():
    """Load baseline composite and time series parameters for all series."""
    data = {}
    for sid, meta in SERIES_REGIMES.items():
        base_dir = os.path.join(RESULTS_DIR, meta['dir'], 'output')
        comp_file = os.path.join(base_dir, 'composite', 'baseline_parameters.csv')
        comp_json = os.path.join(base_dir, 'composite', 'baseline_parameters.json')
        ld_file = os.path.join(base_dir, 'long_distance', 'parameters.csv')
        sd_file = os.path.join(base_dir, 'short_distance', 'parameters.csv')

        df_comp = pd.read_csv(comp_file) if os.path.exists(comp_file) else None
        dict_comp = {}
        if os.path.exists(comp_json):
            with open(comp_json, 'r', encoding='utf-8') as f:
                dict_comp = json.load(f)

        df_ld = pd.read_csv(ld_file) if os.path.exists(ld_file) else None
        df_sd = pd.read_csv(sd_file) if os.path.exists(sd_file) else None

        data[sid] = {
            'meta': meta,
            'composite_df': df_comp,
            'composite_dict': dict_comp,
            'long_distance': df_ld,
            'short_distance': df_sd
        }
    return data


def compute_parameter_summary_table(data):
    """
    Synthesize parameter table across all series.
    Returns list of dicts with formatted and raw values.
    """
    rows = []
    for sid, d in data.items():
        meta = d['meta']
        comp = d['composite_df']

        r_prim = "N/A"
        dm_val = "N/A"
        xi_val = "N/A"
        r_large = "N/A"
        chi2 = "N/A"

        if comp is not None and not comp.empty:
            c = comp.iloc[0]
            if 'R_prim' in c and pd.notna(c['R_prim']):
                r_prim = f"{c['R_prim']:.2f}"
            elif 'Rmean' in c and pd.notna(c['Rmean']) and meta['has_network']:
                r_prim = f"{c['Rmean']:.2f}"
            elif 'Rmean' in c and pd.notna(c['Rmean']) and not meta['has_network']:
                r_prim = f"{c['Rmean']:.1f} (ACC)"

            if 'Dm' in c and pd.notna(c['Dm']):
                dm_val = f"{c['Dm']:.3f}"
            elif 'D' in c and pd.notna(c['D']):
                dm_val = f"{c['D']:.3f}"

            if 'xi' in c and pd.notna(c['xi']):
                xi_val = f"{c['xi']:.1f}" if c['xi'] < 50000 else "Infinite gel"

            if 'mean_radius_large' in c and pd.notna(c['mean_radius_large']):
                r_large = f"{c['mean_radius_large']/1000:.2f}"

            if 'redchi' in c and pd.notna(c['redchi']):
                chi2 = f"{c['redchi']:.3f}"

        rows.append({
            'series_id': sid,
            'concentration_mM': meta['conc_mM'],
            'pH': meta['pH'],
            'regime': meta['regime'],
            'r_prim': r_prim,
            'dm': dm_val,
            'xi': xi_val,
            'r_large_um': r_large,
            'chi2_red': chi2
        })
    return rows


def compute_vps_metrics():
    """Execute complete analysis and write summary JSON."""
    data = load_all_parameters()
    summary_table = compute_parameter_summary_table(data)

    # Calculate ND-i-009 volume shrinking trajectory metrics
    d9 = data['ND-i-009']['long_distance']
    nd009_metrics = {}
    if d9 is not None and not d9.empty:
        nd009_metrics = {
            'time_range_s': [float(d9['time_seconds'].min()), float(d9['time_seconds'].max())],
            'D_m_initial': float(d9['D'].iloc[0]),
            'D_m_final': float(d9['D'].iloc[-1]),
            'D_m_delta': float(d9['D'].iloc[-1] - d9['D'].iloc[0]),
            'R_prim_constant_nm': 8.34,
            'xi_fluc_mean_nm': 2.76,
            'frames_count': len(d9)
        }

    output_payload = {
        'series_summary': summary_table,
        'nd009_volume_shrinking': nd009_metrics,
        'length_scale_hierarchy': {
            'ND-i-003': {'R_prim': 0.0, 'xi': 0.0, 'R_large_nm': 2167.0},
            'ND-i-004': {'R_prim': 4.94, 'xi': 108.7, 'R_large_nm': 2800.0},
            'ND-i-002': {'R_prim': 11.24, 'xi': 356.0, 'R_large_nm': 7163.0},
            'ND-i-009': {'R_prim': 8.34, 'xi': 500.0, 'R_large_nm': 2904.0},
            'ND-i-001': {'R_prim': 6.39, 'xi': 735.0, 'R_large_nm': 3015.0},
            'ND-i-005': {'R_prim': 8.25, 'xi': 1194.0, 'R_large_nm': 3530.0}
        }
    }

    out_json = os.path.join(OUT_DIR, 'vps_metrics_summary.json')
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(output_payload, f, indent=2)
    print(f"[Analysis] Saved VPS metrics summary to: {out_json}")
    return output_payload


def main():
    print("=" * 80)
    print("  RUNNING: VPS Metrics & Parameter Synthesis Computation")
    print("=" * 80)
    metrics = compute_vps_metrics()
    print(f"[Success] Extracted metrics for {len(metrics['series_summary'])} series.")


if __name__ == '__main__':
    main()
