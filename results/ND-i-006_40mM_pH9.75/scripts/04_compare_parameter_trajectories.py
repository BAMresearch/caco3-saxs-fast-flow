#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 04: Unified Cross-Scale Kinetic Trajectories for ND-i-006.
----------------------------------------------------------------
Synthesizes structural evolution across both camera geometries
(Short Distance 3.0 m & Long Distance 31.0 m) from t = 0.060 s to 10.877 s.

Panels:
(A) Structural Length Scales (Primary Building Blocks, Large Droplets, Network Correlation Cutoff)
(B) Mesoscale Spinodal / Mass Fractal Dimension Dm(t) +/- 1 sigma
(C) Volume Fraction Kinetics (phi_net, phi_large) & Porod Surface Growth
(D) Cross-Scale Goodness-of-Fit Metric chi^2_red(t)

Output:
  output/unified_trajectories/unified_cross_scale_trajectories.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'unified_trajectories')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
LONG_CSV = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')


def generate_unified_trajectories():
    print("=" * 80)
    print("  GENERATING UNIFIED CROSS-SCALE TRAJECTORY SYNTHESIS (ND-i-006)")
    print("  Model: Spinodal Mass Fractal Network + Large Droplets (No ACC Particles)")
    print("=" * 80)
    os.makedirs(OUT_DIR, exist_ok=True)

    df_sd = pd.read_csv(SHORT_CSV)
    df_ld = pd.read_csv(LONG_CSV)

    t_sd = df_sd['time_seconds'].values
    t_ld = df_ld['time_seconds'].values

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("ND-i-006 Multiscale Kinetic Synthesis: Spinodal Gelation & Large Droplets\n(40.0 mM CaCO$_3$, pH 9.75, ESRF ID02 Time-Resolved SAXS)", fontsize=13, fontweight='bold', y=0.98)

    # Panel A: Length Scales
    ax1 = axes[0, 0]
    r_prim_val = float(df_sd['R_prim'].iloc[0]) if 'R_prim' in df_sd else 1.50
    std_prim_val = float(df_sd['std_prim'].iloc[0]) if 'std_prim' in df_sd else 3.03
    ax1.axhline(r_prim_val, color='#16a34a', linestyle='--', linewidth=1.6, label=f'Primary Building Blocks ($R_{{\\rm prim}} = {r_prim_val:.1f} \\pm {std_prim_val:.1f}$ nm)')
    ax1.axhspan(max(0.5, r_prim_val - std_prim_val), r_prim_val + std_prim_val, color='#16a34a', alpha=0.15)

    if 'mean_radius_large' in df_ld:
        ax1.plot(t_ld, df_ld['mean_radius_large'], 'd-', color='#0284c7', markersize=3.5, linewidth=1.8, label=r'Large Droplets $\langle R_{\mathrm{large}} \rangle \pm 1\sigma$')
        err_large = df_ld.get('std_err_mean_radius_large', pd.Series(np.zeros_like(t_ld))).values
        ax1.fill_between(t_ld, np.maximum(500.0, df_ld['mean_radius_large'] - err_large), df_ld['mean_radius_large'] + err_large, color='#38bdf8', alpha=0.25)

    ax1.plot(t_ld, df_ld['xi'], '^--', color='#2563eb', markersize=3.5, linewidth=1.6, label=r'Network Cutoff $\xi(t) \pm 1\sigma$')
    err_xi = df_ld.get('std_err_xi', pd.Series(np.zeros_like(t_ld))).values
    ax1.fill_between(t_ld, np.maximum(10.0, df_ld['xi'] - err_xi), df_ld['xi'] + err_xi, color='#2563eb', alpha=0.18)

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylim(1.0, 10000)
    ax1.set_ylabel("Length Scale [nm]", fontsize=11, fontweight='bold')
    ax1.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax1.set_title(r"(A) Hierarchy of Structural Length Scales ($0.060 - 10.88$ s)", fontsize=11, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=8.5, framealpha=0.92)

    # Panel B: Dm
    ax2 = axes[0, 1]
    ax2.plot(t_ld, df_ld['Dm'], 'o-', color='#d97706', markersize=3.5, linewidth=1.8, label=r'Long Distance ($31.0$ m, $D_m \pm 1\sigma$)')
    ax2.fill_between(t_ld, df_ld['Dm'] - df_ld['std_err_Dm'], df_ld['Dm'] + df_ld['std_err_Dm'], color='#d97706', alpha=0.25)
    ax2.plot(t_sd, df_sd['Dm'], 's--', color='#ea580c', markersize=3.0, linewidth=1.4, alpha=0.7, label=r'Short Distance ($3.0$ m, $D_m$)')
    ax2.axhline(1.78, color='#64748b', linestyle=':', label='Diffusion-Limited / DLCA Limit ($D_m \\approx 1.78$)')
    ax2.axhline(2.0, color='#94a3b8', linestyle='--', alpha=0.7, label='Percolated Gel Limit ($D_m = 2.0$)')

    ax2.set_ylabel(r"Mass Fractal Dimension $D_m$", fontsize=11, fontweight='bold')
    ax2.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax2.set_xscale('log')
    ax2.set_ylim(1.1, 2.85)
    ax2.set_title(r"(B) Mesoscale Network Fractal Dimension ($D_m \pm 1\sigma$)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=8.5)

    # Panel C: Volume Fractions & Porod
    ax3 = axes[1, 0]
    ax3_twin = ax3.twinx()
    l1 = ax3.plot(t_ld, df_ld['phi_net'] * 1e4, 's-', color='#16a34a', markersize=3.0, linewidth=1.8, label=r'Network $\phi_{\mathrm{net}} \times 10^4$')
    ax3.fill_between(t_ld, (df_ld['phi_net'] - df_ld['std_err_phi_net']) * 1e4, (df_ld['phi_net'] + df_ld['std_err_phi_net']) * 1e4, color='#16a34a', alpha=0.20)

    if 'total_volume_fraction_large' in df_ld:
        l2 = ax3.plot(t_ld, df_ld['total_volume_fraction_large'] * 1e4, 'o-', color='#0284c7', markersize=3.0, linewidth=1.8, label=r'Large Droplets $\phi_{\mathrm{large}} \times 10^4$')
        err_phi_l = df_ld.get('std_err_total_volume_fraction_large', pd.Series(np.zeros_like(t_ld))).values
        ax3.fill_between(t_ld, (df_ld['total_volume_fraction_large'] - err_phi_l) * 1e4, (df_ld['total_volume_fraction_large'] + err_phi_l) * 1e4, color='#0284c7', alpha=0.22)
    else:
        l2 = []

    l3 = ax3_twin.plot(t_sd, df_sd['A_porod'] * 1e6, 'd:', color='#d97706', markersize=3.0, linewidth=1.6, label=r'Porod Tail $A_{\mathrm{porod}} \times 10^6$')
    ax3_twin.fill_between(t_sd, (df_sd['A_porod'] - df_sd['std_err_A_porod']) * 1e6, (df_sd['A_porod'] + df_sd['std_err_A_porod']) * 1e6, color='#d97706', alpha=0.18)

    ax3.set_xscale('log')
    ax3.set_ylabel(r"Volume Fraction $\phi \times 10^4$", fontsize=11, fontweight='bold', color='#1e293b')
    ax3_twin.set_ylabel(r"Porod Amplitude $A_{\mathrm{porod}} \times 10^6$", fontsize=11, fontweight='bold', color='#d97706')
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax3.set_title(r"(C) Phase Volume Fraction & Porod Surface Kinetics", fontsize=11, fontweight='bold')
    ax3.grid(True, linestyle='--', alpha=0.5)

    all_lines = l1 + (l2 if isinstance(l2, list) else [l2]) + l3
    all_labels = [line.get_label() for line in all_lines]
    ax3.legend(all_lines, all_labels, loc='upper left', fontsize=8.5, framealpha=0.92)

    # Panel D: Chi2
    ax4 = axes[1, 1]
    ax4.plot(t_sd, df_sd['chi2_red'], 'o-', color='#3b82f6', markersize=3.5, linewidth=1.6, label=r'Short-Distance ($3.0$ m, $\chi^2_{\mathrm{red}}$)')
    ax4.plot(t_ld, df_ld['chi2_red'], 's-', color='#ef4444', markersize=3.5, linewidth=1.6, label=r'Long-Distance ($31.0$ m, $\chi^2_{\mathrm{red}}$)')
    ax4.axhline(1.0, color='#64748b', linestyle=':', label=r'Ideal Statistic ($\chi^2_{\mathrm{red}} = 1.0$)')

    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylim(0.08, 15.0)
    ax4.set_ylabel(r"Reduced $\chi^2_{\mathrm{red}}$", fontsize=11, fontweight='bold')
    ax4.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax4.set_title(r"(D) Cross-Scale Goodness-of-Fit Metric ($\chi^2_{\mathrm{red}}$)", fontsize=11, fontweight='bold')
    ax4.grid(True, which='both', linestyle='--', alpha=0.5)
    ax4.legend(loc='upper right', fontsize=8.5, framealpha=0.92)

    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, 'unified_cross_scale_trajectories.png')
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved cross-scale trajectories plot to: {out_png}")


if __name__ == '__main__':
    generate_unified_trajectories()
