#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-Series Trajectory Synthesis and Comparison
------------------------------------------------
Author: Tomasz Stawski, tomasz.stawski@bam.de

This script extracts, aligns, and visualises time-resolved structural
trajectories across the six calcium carbonate reaction series (ND-i-001,
ND-i-002, ND-i-003, ND-i-004, ND-i-005, ND-i-009). ND-i-009 serves as
the primary 40 mM pure CaCO3 benchmark for viscoelastic phase separation.
Outputs:
  results/output/cross_series_trajectories_comparison.png
  results/output/cross_series_trajectories.json
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, RESULTS_DIR)
from saxs_colors import SAXS_COLORS

OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# Series metadata and theoretical maximum volume fractions
# phi_max = [Ca2+]0 * M(CaCO3) / (rho_ACC * 1000)
# 10 mM -> 5.0e-4; 15 mM -> 7.5e-4; 20 mM -> 1.0e-3; 40 mM -> 2.0e-3
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
        if sid == 'ND-i-003':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['Rmean'].values
            r_meso_label = r'ACC radius $\langle R_{\rm ACC} \rangle$'
            r_prim = np.nan
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 2167.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 5.1e-5)
        elif sid == 'ND-i-004':
            phi_vals = (df_ld['phi_net'].fillna(0.0) + df_ld['phi_acc'].fillna(0.0)).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 4.94
            r_drop = np.full_like(t_vals, 2800.0)
            phi_drop = np.full_like(t_vals, 4.5e-5)
        elif sid == 'ND-i-002':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 11.20
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 7163.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 8.4e-5)
        elif sid == 'ND-i-009':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 2.80
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 2904.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 5.9e-4)
        elif sid == 'ND-i-001':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 3.50
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 3015.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 4.98e-5)
        elif sid == 'ND-i-005':
            phi_vals = df_ld['phi'].fillna(0.0).values
            r_meso = df_ld['xi'].values if 'xi' in df_ld.columns else df_ld['xi_fluc'].values
            r_meso_label = r'Correlation cutoff $\xi$'
            r_prim = 6.59
            r_drop = df_ld['mean_radius_large'].values if 'mean_radius_large' in df_ld.columns else np.full_like(t_vals, 3530.0)
            phi_drop = df_ld['total_volume_fraction_large'].values if 'total_volume_fraction_large' in df_ld.columns else np.full_like(t_vals, 2.92e-5)

        chi2_col = 'redchi' if 'redchi' in df_ld.columns else ('chi2_red' if 'chi2_red' in df_ld.columns else None)
        chi2_vals = df_ld[chi2_col].values if chi2_col else np.ones_like(t_vals)

        alpha_vals = phi_vals / cfg['phi_max']
        alpha_vals = np.clip(alpha_vals, 0.0, 1.2)

        trajectories[sid] = {
            'config': cfg,
            'time': t_vals,
            'phi': phi_vals,
            'alpha': alpha_vals,
            'r_meso': r_meso,
            'r_meso_label': r_meso_label,
            'r_prim': r_prim,
            'r_drop': r_drop,
            'phi_drop': phi_drop,
            'chi2': chi2_vals
        }

    return trajectories


def generate_synthesis_figure(trajectories):
    """
    Generate 4-panel publication figure showing cross-series kinetic evolution.
    Panel A: Precipitation extent alpha(t) vs time.
    Panel B: Characteristic meso-scale length scale (ACC radius vs network cutoff xi).
    Panel C: Structural bifurcation across concentration (spheres vs network).
    Panel D: Macroscopic phase inversion droplet population.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5), dpi=250)
    fig.patch.set_facecolor('#ffffff')

    for ax in axes.flat:
        ax.set_facecolor('#fafafa')
        ax.tick_params(direction='in', which='both', top=True, right=True)

    # -------------------------------------------------------------------------
    # Panel A: Precipitation Extent alpha(t) = phi(t) / phi_max
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    for sid, d in trajectories.items():
        cfg = d['config']
        t = d['time']
        alpha = d['alpha']
        step = 1 if len(t) <= 60 else (len(t) // 60)
        ax.plot(
            t[::step], alpha[::step],
            marker=cfg['marker'], color=cfg['colour'], linewidth=2.0, markersize=5.0,
            label=f"{cfg['name']} ({cfg['regime']})"
        )

    ax.set_xscale('log')
    ax.set_xlim(0.04, 12.0)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("Reaction time $t$ [s]", fontsize=11)
    ax.set_ylabel(r"Conversion $\alpha(t) = \phi(t) / \phi_{\rm max}$", fontsize=11)
    ax.set_title(r"(a) Precipitation extent: induction and growth kinetics", fontsize=11, fontweight='bold')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend(fontsize=8.5, loc='upper left', framealpha=0.92)

    # -------------------------------------------------------------------------
    # Panel B: Dominant Meso-Scale Length Scale (R_acc vs xi)
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    for sid, d in trajectories.items():
        cfg = d['config']
        t = d['time']
        r = d['r_meso']
        step = 1 if len(t) <= 60 else (len(t) // 60)
        valid = np.isfinite(r) & (r > 0)
        if np.any(valid):
            ax.plot(
                t[valid][::step], r[valid][::step],
                marker=cfg['marker'], color=cfg['colour'], linewidth=2.0, markersize=5.0,
                label=f"{sid}: {d['r_meso_label']}"
            )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.04, 12.0)
    ax.set_ylim(15, 2500)
    ax.set_xlabel("Reaction time $t$ [s]", fontsize=11)
    ax.set_ylabel("Characteristic meso-scale length [nm]", fontsize=11)
    ax.set_title(r"(b) Meso-scale length scale: ACC radius vs network cutoff $\xi$", fontsize=11, fontweight='bold')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend(fontsize=8.5, loc='lower right', framealpha=0.92)

    # -------------------------------------------------------------------------
    # Panel C: Structural Bifurcation Across Initial Concentration [Ca2+]0
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    ax.axvspan(12.0, 15.5, color='#fed7aa', alpha=0.45, label=r'Spinodal boundary $C^* \approx 12 - 15$ mM')

    # Pure CaCO3 series
    ax.plot([10, 15], [107.0, 47.7], 's-', color='#2563eb', linewidth=2.2, markersize=8, label=r'Discrete ACC nanospheres $\langle R_{\rm ACC} \rangle$')
    ax.plot([15, 20, 40], [4.94, 11.20, 2.80], 'o-', color='#16a34a', linewidth=2.2, markersize=8, label=r'Primary building blocks $R_{\rm prim}$ (pure $\mathrm{CaCO}_3$)')
    ax.plot([15, 20, 40], [108.7, 356.0, 735.0], '^--', color='#dc2626', linewidth=2.0, markersize=8, label=r'Network correlation cutoff $\xi$ (pure $\mathrm{CaCO}_3$)')

    # Additive Dopant annotations
    ax.plot(40.0, 6.59, 'v', color='#9333ea', markersize=9, label=r'Co$^{2+}$ additive: enlarged node ($6.59$ nm)')
    ax.plot(40.0, 1194.0, 'v', color='#9333ea', markersize=9, label=r'Co$^{2+}$ additive: expanded cutoff $\xi$ ($1,194$ nm)')

    ax.set_yscale('log')
    ax.set_xlim(8, 44)
    ax.set_ylim(2, 2500)
    ax.set_xlabel(r"Initial concentration $[\mathrm{Ca}^{2+}]_0$ [mM]", fontsize=11)
    ax.set_ylabel("Length scale [nm]", fontsize=11)
    ax.set_title(r"(c) Structural bifurcation across concentration: spheres vs network", fontsize=11, fontweight='bold')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend(fontsize=8.0, loc='lower left', framealpha=0.92)

    # Inset: Stacked relative phase proportions across concentration regimes
    ax_in = ax.inset_axes([0.54, 0.50, 0.43, 0.45])
    concs = ['10 mM', '15 mM', '20 mM', '40 mM']
    acc_pct = np.array([100.0, 63.1, 0.0, 0.0])
    net_pct = np.array([0.0, 36.9, 75.9, 72.9])
    pnc_pct = np.array([0.0, 0.0, 24.1, 25.9])
    drp_pct = np.array([0.0, 0.0, 0.0, 1.2])
    
    x_in = np.arange(len(concs))
    w_in = 0.55
    ax_in.bar(x_in, acc_pct, w_in, color=SAXS_COLORS['acc_spheres'], label='ACC Spheres')
    ax_in.bar(x_in, net_pct, w_in, bottom=acc_pct, color=SAXS_COLORS['fractal_network'], label='Fractal Gel')
    ax_in.bar(x_in, pnc_pct, w_in, bottom=acc_pct + net_pct, color=SAXS_COLORS['pnc_fluctuations'], label='Tanaka PNC Mesh')
    ax_in.bar(x_in, drp_pct, w_in, bottom=acc_pct + net_pct + pnc_pct, color=SAXS_COLORS['droplets_porod'], label='Droplets / Porod')
    
    ax_in.set_xticks(x_in)
    ax_in.set_xticklabels(['10', '15', '20', '40'], fontsize=7.5)
    ax_in.set_xlabel(r'$[\mathrm{Ca}^{2+}]_0$ [mM]', fontsize=7.5)
    ax_in.set_ylim(0, 105)
    ax_in.set_ylabel('Share [%]', fontsize=7.5)
    ax_in.set_title('Phase Dominance', fontsize=8.0, fontweight='bold')
    ax_in.tick_params(labelsize=7.5)
    ax_in.legend(fontsize=6.0, loc='upper right', framealpha=0.88)

    # -------------------------------------------------------------------------
    # Panel D: Universal Phase Inversion Droplet Population Across Series
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    series_names = ['ND-i-003', 'ND-i-004', 'ND-i-002', 'ND-i-009', 'ND-i-001', 'ND-i-005']
    r_drop_mean = []
    phi_drop_mean = []

    for sid in series_names:
        d = trajectories[sid]
        r_drop_mean.append(float(np.nanmedian(d['r_drop'])))
        phi_drop_mean.append(float(np.nanmedian(d['phi_drop'])) * 1e5)

    x_pos = np.arange(len(series_names))
    width = 0.35

    color_bars = [trajectories[s]['config']['colour'] for s in series_names]
    ax.bar(x_pos - width/2, r_drop_mean, width, color=color_bars, alpha=0.85, edgecolor='black', label=r'Droplet radius $\langle R_{\rm large} \rangle$ [nm]')
    
    ax2 = ax.twinx()
    ax2.bar(x_pos + width/2, phi_drop_mean, width, color='#64748b', alpha=0.55, hatch='//', edgecolor='black', label=r'Volume fraction $\phi_{\rm large} \times 10^5$')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(series_names, fontsize=9.5)
    ax.set_ylabel(r'Droplet radius $\langle R_{\rm large} \rangle$ [nm]', fontsize=11, color='#1e293b')
    ax2.set_ylabel(r'Volume fraction $\phi_{\rm large} \times 10^5$', fontsize=11, color='#475569')
    ax.set_ylim(0, 8000)
    ax2.set_ylim(0, 150)
    ax.set_title(r"(d) Stationary phase inversion droplets across series", fontsize=11, fontweight='bold')
    ax.grid(True, axis='y', linestyle=':', alpha=0.6)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8.5, framealpha=0.92)

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'cross_series_trajectories_comparison.png')
    plt.savefig(plot_path, dpi=250)
    plt.close(fig)
    print(f"Saved cross-series trajectory comparison plot to: {plot_path}")


def export_trajectories_json(trajectories):
    """Export aligned trajectory data to JSON format."""
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
    print(f"Saved cross-series trajectories JSON to: {json_path}")


def main():
    trajectories = load_series_trajectories()
    generate_synthesis_figure(trajectories)
    export_trajectories_json(trajectories)


if __name__ == '__main__':
    main()
