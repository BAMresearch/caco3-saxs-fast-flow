#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 04: Unified Cross-Scale Kinetic Trajectories for ND-i-001.
----------------------------------------------------------------
Synthesizes the multiscale structural evolution across both camera geometries
(Short Distance 3.0 m & Long Distance 31.0 m) from t = 0.060 s to 10.877 s.

Panels:
(A) Structural Hierarchy Across 4 Orders of Magnitude (0.5 nm to 4.0 um)
(B) Cross-Scale Network Fractal Dimension Agreement D_f(t) +/- 1 sigma
(C) Phase Volume Fraction Evolution (phi_prim, phi_large, A_fluc) +/- 1 sigma
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
    print("  GENERATING UNIFIED CROSS-SCALE TRAJECTORY SYNTHESIS (ND-i-001)")
    print("=" * 80)
    os.makedirs(OUT_DIR, exist_ok=True)

    df_sd = pd.read_csv(SHORT_CSV)
    df_ld = pd.read_csv(LONG_CSV)

    t_sd = df_sd['time_seconds'].values
    t_ld = df_ld['time_seconds'].values

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.28, wspace=0.24)
    fig.suptitle("ND-i-001 Multiscale Kinetic Synthesis: Time-Resolved Phase Separation\n(40 mM CaCO$_3$, pH 11.3, Pure System)", fontsize=13, fontweight='bold', y=0.98)

    # Panel A: Hierarchy of Structural Length Scales Across 4 Orders of Magnitude
    ax1 = fig.add_subplot(gs[0, :])
    
    # 1. Viscoelastic fluctuation mesh correlation length
    ax1.plot(t_sd, df_sd['xi_fluc'], 's-', color='#e63946', markersize=3.5, linewidth=1.8, label=r'Viscoelastic Fluctuation Mesh $\xi_{\mathrm{fluc}} \pm 1\sigma$')
    ax1.fill_between(t_sd, df_sd['xi_fluc'] - df_sd['std_err_xi_fluc'], df_sd['xi_fluc'] + df_sd['std_err_xi_fluc'], color='#e63946', alpha=0.25)

    # 2. Primary building blocks
    ax1.plot(t_sd, df_sd['Rmean'], 'o-', color='#2a9d8f', markersize=3.5, linewidth=1.8, label=r'Primary Building Block $\langle R_{\mathrm{prim}} \rangle \pm 1\sigma$')
    ax1.fill_between(t_sd, df_sd['Rmean'] - df_sd['std_err_Rmean'], df_sd['Rmean'] + df_sd['std_err_Rmean'], color='#2a9d8f', alpha=0.25)

    # 3. Mesoscale network correlation length
    ax1.plot(t_ld, df_ld['xi'], '^--', color='#9d4edd', markersize=3.5, linewidth=1.6, label=r'Network Cutoff $\xi(t) \pm 1\sigma$')
    ax1.fill_between(t_ld, np.maximum(10.0, df_ld['xi'] - df_ld['std_err_xi']), df_ld['xi'] + df_ld['std_err_xi'], color='#9d4edd', alpha=0.18)

    # 4. Large droplet radius
    ax1.plot(t_ld, df_ld['mean_radius_large'], 'D-', color='#1d3557', markersize=3.5, linewidth=1.8, label=r'Secondary Droplets $\langle R_{\mathrm{large}} \rangle \pm 1\sigma$')
    ax1.fill_between(t_ld, df_ld['mean_radius_large'] - df_ld['std_err_mean_radius_large'], df_ld['mean_radius_large'] + df_ld['std_err_mean_radius_large'], color='#1d3557', alpha=0.25)

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylim(0.8, 8000)
    ax1.set_ylabel("Length Scale [nm]", fontsize=11, fontweight='bold')
    ax1.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax1.set_title(r"(A) Hierarchy of Structural Length Scales Across Real Time ($0.060 - 10.88$ s) $\pm 1\sigma$", fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='center right', fontsize=9.5, framealpha=0.92)

    # Panel B: Cross-Scale Network Dimension Agreement
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(t_sd, df_sd['D'], 'o-', color='#d62728', markersize=3.0, linewidth=1.8, label=r'Short Distance ($3.0$ m, local mesh)')
    ax2.fill_between(t_sd, df_sd['D'] - df_sd['std_err_D'], df_sd['D'] + df_sd['std_err_D'], color='#d62728', alpha=0.25, label=r'SD $\pm 1\sigma$ band')

    ax2.plot(t_ld, df_ld['D'], 's-', color='#1f77b4', markersize=3.0, linewidth=1.8, label=r'Long Distance ($31.0$ m, mesoscale mesh)')
    ax2.fill_between(t_ld, df_ld['D'] - df_ld['std_err_D'], df_ld['D'] + df_ld['std_err_D'], color='#1f77b4', alpha=0.20, label=r'LD $\pm 1\sigma$ band')

    ax2.set_ylabel(r"Mass Fractal Dimension $D_f$", fontsize=11, fontweight='bold')
    ax2.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax2.set_xscale('log')
    ax2.set_title(r"(B) Cross-Scale Network Dimension Agreement ($D_f \pm 1\sigma$)", fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='best', fontsize=9.0)

    # Panel C: Volume Fraction & Fluctuation Evolution
    ax3 = fig.add_subplot(gs[1, 1])
    ax3_twin = ax3.twinx()
    l1 = ax3.plot(t_sd, df_sd['phi'] * 1e4, 'o-', color='#2a9d8f', markersize=3.0, linewidth=1.8, label=r'Primary $\phi_{\mathrm{prim}} \times 10^4$')
    ax3.fill_between(t_sd, (df_sd['phi'] - df_sd['std_err_phi']) * 1e4, (df_sd['phi'] + df_sd['std_err_phi']) * 1e4, color='#2a9d8f', alpha=0.22)

    l2 = ax3.plot(t_ld, df_ld['total_volume_fraction_large'] * 1e5, 'D-', color='#1d3557', markersize=3.0, linewidth=1.8, label=r'Droplet $\phi_{\mathrm{large}} \times 10^5$')
    ax3.fill_between(t_ld, (df_ld['total_volume_fraction_large'] - df_ld['std_err_total_volume_fraction_large']) * 1e5,
                     (df_ld['total_volume_fraction_large'] + df_ld['std_err_total_volume_fraction_large']) * 1e5, color='#1d3557', alpha=0.22)

    l3 = ax3_twin.plot(t_sd, df_sd['A_fluc'], '^--', color='#e76f51', markersize=3.0, linewidth=1.6, label=r'Fluct. Amp. $A_{\mathrm{fluc}}$')
    ax3_twin.fill_between(t_sd, df_sd['A_fluc'] - df_sd['std_err_A_fluc'], df_sd['A_fluc'] + df_sd['std_err_A_fluc'], color='#e76f51', alpha=0.20)

    ax3.set_ylabel(r"Volume Fraction Scale", fontsize=11, fontweight='bold')
    ax3_twin.set_ylabel(r"$A_{\mathrm{fluc}}$ [cm$^{-1}$]", fontsize=11, fontweight='bold', color='#e76f51')
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    ax3.set_xscale('log')
    ax3.set_title(r"(C) Phase Volume Fraction & Fluctuation Kinetics $\pm 1\sigma$", fontsize=12, fontweight='bold')
    ax3.grid(True, linestyle='--', alpha=0.5)
    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='best', fontsize=8.5)

    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, 'unified_cross_scale_trajectories.png')
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved cross-scale trajectories figure: {out_png}")
    print("=" * 80)


if __name__ == '__main__':
    generate_unified_trajectories()
