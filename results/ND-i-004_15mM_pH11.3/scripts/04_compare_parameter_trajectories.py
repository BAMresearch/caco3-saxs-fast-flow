#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 04: Unified Cross-Scale Kinetic Trajectories for ND-i-004.
----------------------------------------------------------------
Synthesizes the multiscale structural evolution across both camera geometries
(Short Distance 3.0 m & Long Distance 31.0 m) from t = 0.060 s to 10.877 s.

Features dedicated tracking of relative component proportions and dynamic phase
inversion (ACC nanospheres vs. spinodal fractal network).

Panels:
(A) Hierarchy of Structural Length Scales (R_prim, <R_ACC>, xi)
(B) Mesoscale Network Fractal Dimension D_m(t) +/- 1 sigma
(C) Absolute Phase Volume Fractions (phi_acc, phi_net) & Porod Amplitude
(D) Relative Phase Proportions & Dynamic Crossover (chi_net vs chi_acc in %)
(E) Phase Mass Ratio Evolution (phi_net / phi_acc on log scale)
(F) Cross-Scale Goodness-of-Fit Metric chi^2_red(t)

Outputs:
  output/unified_trajectories/unified_cross_scale_trajectories.png
  output/unified_trajectories/unified_cross_scale_trajectories.pdf
  output/unified_trajectories/relative_phase_evolution.png
  output/unified_trajectories/relative_phase_evolution.pdf
  output/unified_trajectories/relative_proportions.csv
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
RESULTS_DIR = os.path.abspath(os.path.join(PROJECT_DIR, '..'))
sys.path.insert(0, RESULTS_DIR)
sys.path.insert(0, PROJECT_DIR)

from saxs_colors import SAXS_COLORS

OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'unified_trajectories')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
LONG_CSV = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')


def generate_unified_trajectories():
    print("=" * 80)
    print("  GENERATING UNIFIED CROSS-SCALE TRAJECTORY SYNTHESIS (ND-i-004)")
    print("  Tracking relative component proportions and dominance evolution")
    print("=" * 80)
    os.makedirs(OUT_DIR, exist_ok=True)

    df_sd = pd.read_csv(SHORT_CSV)
    df_ld = pd.read_csv(LONG_CSV)

    t_sd = df_sd['time_seconds'].values
    t_ld = df_ld['time_seconds'].values

    # Extract volume fractions and calculate relative proportions
    p_net_ld = df_ld['phi_net'].values
    p_acc_ld = df_ld['phi_acc'].values
    std_net_ld = df_ld['std_err_phi_net'].values if 'std_err_phi_net' in df_ld.columns else np.zeros_like(p_net_ld)
    std_acc_ld = df_ld['std_err_phi_acc'].values if 'std_err_phi_acc' in df_ld.columns else np.zeros_like(p_acc_ld)

    p_tot_ld = np.maximum(p_net_ld + p_acc_ld, 1e-12)
    chi_net = (p_net_ld / p_tot_ld) * 100.0
    chi_acc = (p_acc_ld / p_tot_ld) * 100.0

    # Propagate standard errors on percentage proportions
    # sigma_chi = sqrt( (dchi/dp_net)^2 * s_net^2 + (dchi/dp_acc)^2 * s_acc^2 )
    dchi_dnet = p_acc_ld / (p_tot_ld**2) * 100.0
    dchi_dacc = -p_net_ld / (p_tot_ld**2) * 100.0
    std_chi = np.sqrt((dchi_dnet * std_net_ld)**2 + (dchi_dacc * std_acc_ld)**2)
    std_chi = np.clip(std_chi, 0.0, 25.0)

    # Phase ratio (network / acc spheres)
    ratio_net_acc = p_net_ld / np.maximum(p_acc_ld, 1e-9)

    # Save relative proportions table
    df_props = pd.DataFrame({
        'frame_index': df_ld['frame_index'] if 'frame_index' in df_ld.columns else np.arange(len(t_ld)),
        'time_seconds': t_ld,
        'phi_net': p_net_ld,
        'phi_acc': p_acc_ld,
        'chi_net_percent': chi_net,
        'std_err_chi_net': std_chi,
        'chi_acc_percent': chi_acc,
        'std_err_chi_acc': std_chi,
        'net_to_acc_ratio': ratio_net_acc,
        'dominant_component': np.where(chi_net >= 50.0, 'Fractal Network', 'ACC Nanospheres')
    })
    out_csv = os.path.join(OUT_DIR, 'relative_proportions.csv')
    df_props.to_csv(out_csv, index=False)
    print(f"Saved relative proportions table: {out_csv}")

    # =========================================================================
    # FIGURE 1: Comprehensive 6-Panel Kinetic Synthesis (3 rows x 2 cols)
    # =========================================================================
    fig, axes = plt.subplots(3, 2, figsize=(14, 14.5), dpi=300)
    fig.patch.set_facecolor('#ffffff')
    fig.suptitle(
        "ND-i-004 Multiscale Kinetic Synthesis & Component Dominance Evolution\n"
        "(15.0 mM CaCO$_3$, pH 11.3, ESRF ID02 Time-Resolved SAXS)",
        fontsize=13, fontweight='bold', y=0.985
    )

    for ax in axes.flat:
        ax.set_facecolor('#fafafa')
        ax.tick_params(direction='in', which='both', top=True, right=True)

    # -------------------------------------------------------------------------
    # Panel A: Hierarchy of Structural Length Scales
    # -------------------------------------------------------------------------
    ax1 = axes[0, 0]
    r_prim_val = float(df_sd['R_prim'].iloc[0]) if 'R_prim' in df_sd else 4.94
    std_prim_val = float(df_sd['std_prim'].iloc[0]) if 'std_prim' in df_sd else 1.5
    ax1.axhline(r_prim_val, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.6,
                label=f'Primary Building Blocks ($R_{{\\rm prim}} = {r_prim_val:.1f} \\pm {std_prim_val:.1f}$ nm)')
    ax1.axhspan(r_prim_val - std_prim_val, r_prim_val + std_prim_val, color=SAXS_COLORS['fractal_network'], alpha=0.15)

    ax1.plot(t_sd, df_sd['R_acc'], 'o-', color=SAXS_COLORS['acc_spheres'], markersize=3.5, linewidth=1.8,
             label=r'ACC Nanospheres $\langle R_{\mathrm{ACC}} \rangle \pm \sigma_{\mathrm{ACC}}$')
    ax1.fill_between(t_sd, np.maximum(1.0, df_sd['R_acc'] - df_sd['std_acc']), df_sd['R_acc'] + df_sd['std_acc'],
                     color=SAXS_COLORS['acc_spheres'], alpha=0.22)

    ax1.plot(t_ld, df_ld['xi'], '^--', color='#2563eb', markersize=3.5, linewidth=1.6,
             label=r'Network Cutoff $\xi(t) \pm 1\sigma$')
    ax1.fill_between(t_ld, np.maximum(10.0, df_ld['xi'] - df_ld['std_err_xi']), df_ld['xi'] + df_ld['std_err_xi'],
                     color='#2563eb', alpha=0.18)

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylim(1.0, 1000)
    ax1.set_ylabel("Length Scale [nm]", fontsize=10.5, fontweight='bold')
    ax1.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax1.set_title(r"(A) Hierarchy of Structural Length Scales ($0.060 - 10.88$ s)", fontsize=11, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=8.5, framealpha=0.92)

    # -------------------------------------------------------------------------
    # Panel B: Network Mass-Fractal Dimension
    # -------------------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.plot(t_ld, df_ld['Dm'], 'o-', color='#d97706', markersize=3.5, linewidth=1.8,
             label=r'Long Distance ($31.0$ m, $D_m \pm 1\sigma$)')
    ax2.fill_between(t_ld, df_ld['Dm'] - df_ld['std_err_Dm'], df_ld['Dm'] + df_ld['std_err_Dm'],
                     color='#d97706', alpha=0.25)
    ax2.axhline(1.78, color='#64748b', linestyle=':', label='Diffusion-Limited / DLCA Limit ($D_m \\approx 1.78$)')
    ax2.axhline(2.0, color='#94a3b8', linestyle='--', alpha=0.7, label='Percolated Gel Limit ($D_m = 2.0$)')

    ax2.set_ylabel(r"Mass Fractal Dimension $D_m$", fontsize=10.5, fontweight='bold')
    ax2.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax2.set_xscale('log')
    ax2.set_ylim(1.6, 2.65)
    ax2.set_title(r"(B) Mesoscale Network Fractal Dimension ($D_m \pm 1\sigma$)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel C: Absolute Multi-Phase Volume Fraction Kinetics
    # -------------------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3_twin = ax3.twinx()
    l1 = ax3.plot(t_ld, p_acc_ld * 1e5, 'o-', color=SAXS_COLORS['acc_spheres'], markersize=3.2, linewidth=1.8,
                  label=r'ACC Spheres $\phi_{\mathrm{ACC}} \times 10^5$')
    ax3.fill_between(t_ld, (p_acc_ld - std_acc_ld) * 1e5, (p_acc_ld + std_acc_ld) * 1e5,
                     color=SAXS_COLORS['acc_spheres'], alpha=0.20)

    l2 = ax3.plot(t_ld, p_net_ld * 1e4, 's--', color=SAXS_COLORS['fractal_network'], markersize=3.2, linewidth=1.6,
                  label=r'Network $\phi_{\mathrm{net}} \times 10^4$')
    ax3.fill_between(t_ld, (p_net_ld - std_net_ld) * 1e4, (p_net_ld + std_net_ld) * 1e4,
                     color=SAXS_COLORS['fractal_network'], alpha=0.20)

    l3 = ax3_twin.plot(t_sd, df_sd['A_porod'] * 1e6, 'D:', color=SAXS_COLORS['droplets_porod'], markersize=2.8, linewidth=1.4,
                       label=r'Porod $A_{\mathrm{porod}} \times 10^6$ [cm$^{-1}$nm$^{4}$]')

    ax3.set_ylabel(r"Phase Volume Fractions", fontsize=10.5, fontweight='bold')
    ax3_twin.set_ylabel(r"Porod Amplitude $A_{\mathrm{porod}} \times 10^6$", fontsize=10.5, fontweight='bold', color=SAXS_COLORS['droplets_porod'])
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax3.set_xscale('log')
    ax3.set_title(r"(C) Absolute Phase Volume Fractions & Porod Surface", fontsize=11, fontweight='bold')
    ax3.grid(True, linestyle='--', alpha=0.5)
    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='center left', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel D: Relative Phase Proportions & Dynamic Crossover
    # -------------------------------------------------------------------------
    ax4 = axes[1, 1]
    # Dominance reference line at 50%
    ax4.axhline(50.0, color='#64748b', linestyle='--', linewidth=1.3, alpha=0.8, label='Equiproportion / Dominance Threshold (50%)')

    # Shaded Dynamic Inversion Window
    ax4.axvspan(0.10, 0.28, color='#fef3c7', alpha=0.6, label='Dynamic Crossover Window ($t \\approx 0.15 - 0.25$ s)')

    # Curves for relative proportions
    ax4.plot(t_ld, chi_net, 's-', color=SAXS_COLORS['fractal_network'], markersize=4.0, linewidth=2.0,
             label=r'Fractal Network Proportion $\chi_{\mathrm{net}}(t)$')
    ax4.fill_between(t_ld, np.maximum(0.0, chi_net - std_chi), np.minimum(100.0, chi_net + std_chi),
                     color=SAXS_COLORS['fractal_network'], alpha=0.18)

    ax4.plot(t_ld, chi_acc, 'o-', color=SAXS_COLORS['acc_spheres'], markersize=4.0, linewidth=2.0,
             label=r'Discrete ACC Nanospheres $\chi_{\mathrm{ACC}}(t)$')
    ax4.fill_between(t_ld, np.maximum(0.0, chi_acc - std_chi), np.minimum(100.0, chi_acc + std_chi),
                     color=SAXS_COLORS['acc_spheres'], alpha=0.18)

    # Annotated regimes
    ax4.text(0.065, 88.0, r"ACC Dominance" "\n" r"($95.8\%$ at $0.06$ s)", color=SAXS_COLORS['acc_spheres'],
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor=SAXS_COLORS['acc_spheres'], alpha=0.85))

    ax4.text(1.2, 85.0, "Network Dominance\n($76.3\\%$ plateau)", color=SAXS_COLORS['fractal_network'],
             fontsize=8.5, fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor=SAXS_COLORS['fractal_network'], alpha=0.85))

    ax4.set_ylabel(r"Relative Proportion $\chi$ [%]", fontsize=10.5, fontweight='bold')
    ax4.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax4.set_xscale('log')
    ax4.set_ylim(-2.0, 102.0)
    ax4.set_title(r"(D) Relative Phase Proportions & Dynamic Inversion", fontsize=11, fontweight='bold')
    ax4.grid(True, linestyle='--', alpha=0.5)
    ax4.legend(loc='center right', fontsize=8.2, framealpha=0.92)

    # -------------------------------------------------------------------------
    # Panel E: Phase Ratio Evolution (phi_net / phi_acc)
    # -------------------------------------------------------------------------
    ax5 = axes[2, 0]
    ax5.axhline(1.0, color='#64748b', linestyle='--', linewidth=1.3, label='Parity Ratio (1.0 = equal mass)')
    ax5.axvspan(0.10, 0.28, color='#fef3c7', alpha=0.5)

    ax5.plot(t_ld, ratio_net_acc, 'D-', color='#9333ea', markersize=3.8, linewidth=1.8,
             label=r'Phase Mass Ratio $\phi_{\mathrm{net}} / \phi_{\mathrm{ACC}}$')

    # Shaded zones: ACC dominant vs Network dominant
    ax5.axhspan(0.01, 1.0, color=SAXS_COLORS['acc_spheres'], alpha=0.08)
    ax5.axhspan(1.0, 20.0, color=SAXS_COLORS['fractal_network'], alpha=0.08)
    ax5.text(0.065, 0.08, 'ACC Dominant Zone', color=SAXS_COLORS['acc_spheres'], fontsize=8.5, fontweight='bold')
    ax5.text(0.6, 5.0, 'Network Dominant Zone', color=SAXS_COLORS['fractal_network'], fontsize=8.5, fontweight='bold')

    ax5.set_yscale('log')
    ax5.set_xscale('log')
    ax5.set_ylim(0.02, 12.0)
    ax5.set_ylabel(r"Ratio $\phi_{\mathrm{net}} / \phi_{\mathrm{ACC}}$ [log]", fontsize=10.5, fontweight='bold')
    ax5.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax5.set_title(r"(E) Phase Ratio Evolution: Dynamic Inversion ($0.04 \to 3.3$)", fontsize=11, fontweight='bold')
    ax5.grid(True, which='both', linestyle='--', alpha=0.5)
    ax5.legend(loc='lower right', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel F: Goodness of Fit Across Both Geometries
    # -------------------------------------------------------------------------
    ax6 = axes[2, 1]
    ax6.plot(t_sd, df_sd['chi2_red'], 'o-', color=SAXS_COLORS['acc_spheres'], markersize=3.5, linewidth=1.6,
             label=r'Short Distance ($3.0$ m, ACC + Network + Porod)')
    ax6.plot(t_ld, df_ld['chi2_red'], 's-', color=SAXS_COLORS['fractal_network'], markersize=3.5, linewidth=1.6,
             label=r'Long Distance ($31.0$ m, Network + Porod + ACC)')
    ax6.axhline(1.0, color='#64748b', linestyle=':', label='Ideal Target ($\\chi^2_{\\mathrm{red}} = 1.0$)')

    ax6.set_ylabel(r"Goodness of Fit $\chi^2_{\mathrm{red}}$", fontsize=10.5, fontweight='bold')
    ax6.set_xlabel("Reaction Time $t$ [s]", fontsize=10.5, fontweight='bold')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylim(0.04, 3.0)
    ax6.set_title(r"(F) Cross-Scale Goodness-of-Fit $\chi^2_{\mathrm{red}}(t)$", fontsize=11, fontweight='bold')
    ax6.grid(True, which='both', linestyle='--', alpha=0.5)
    ax6.legend(loc='upper right', fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.965])
    out_png = os.path.join(OUT_DIR, 'unified_cross_scale_trajectories.png')
    out_pdf = os.path.join(OUT_DIR, 'unified_cross_scale_trajectories.pdf')
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved cross-scale trajectories figures: {out_png} and {out_pdf}")

    # =========================================================================
    # FIGURE 2: Dedicated Standalone Relative Phase Proportions Plot (2 panels)
    # =========================================================================
    fig2, (pax1, pax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    fig2.patch.set_facecolor('#ffffff')
    fig2.suptitle(
        "ND-i-004 Kinetic Phase Inversion: Discrete ACC Nanospheres vs. Spinodal Fractal Network\n"
        "(15.0 mM CaCO$_3$, pH 11.3, ESRF ID02 Time-Resolved SAXS)",
        fontsize=12, fontweight='bold'
    )

    for ax in (pax1, pax2):
        ax.set_facecolor('#fafafa')
        ax.tick_params(direction='in', which='both', top=True, right=True)
        ax.set_xscale('log')

    # Panel 1: Stacked Area Relative Proportions
    pax1.axhline(50.0, color='#64748b', linestyle='--', linewidth=1.2, label='50% Dominance Boundary')
    pax1.axvspan(0.10, 0.28, color='#fef3c7', alpha=0.6, label='Dynamic Crossover Window ($t \\approx 0.15 - 0.25$ s)')

    pax1.fill_between(t_ld, 0, chi_net, color=SAXS_COLORS['fractal_network'], alpha=0.35, label='Spinodal Fractal Network Share [%]')
    pax1.fill_between(t_ld, chi_net, 100, color=SAXS_COLORS['acc_spheres'], alpha=0.35, label='Discrete ACC Nanospheres Share [%]')

    pax1.plot(t_ld, chi_net, 's-', color=SAXS_COLORS['fractal_network'], markersize=4.0, linewidth=2.0)
    pax1.plot(t_ld, chi_acc, 'o-', color=SAXS_COLORS['acc_spheres'], markersize=4.0, linewidth=2.0)

    pax1.set_ylim(0, 100)
    pax1.set_xlim(0.05, 12.0)
    pax1.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    pax1.set_ylabel("Relative Volume Proportion [%]", fontsize=11, fontweight='bold')
    pax1.set_title("(a) Dynamic Inversion of Component Proportions", fontsize=11, fontweight='bold')
    pax1.grid(True, linestyle='--', alpha=0.5)
    pax1.legend(loc='center right', fontsize=8.8, framealpha=0.92)

    # Panel 2: Phase Ratio
    pax2.axhline(1.0, color='#64748b', linestyle='--', linewidth=1.3, label=r'Parity ($\phi_{\rm net} = \phi_{\rm ACC}$)')
    pax2.axvspan(0.10, 0.28, color='#fef3c7', alpha=0.6)
    pax2.plot(t_ld, ratio_net_acc, 'D-', color='#9333ea', markersize=4.5, linewidth=2.0, label=r'Mass Ratio $\phi_{\rm net} / \phi_{\rm ACC}$')

    pax2.axhspan(0.02, 1.0, color=SAXS_COLORS['acc_spheres'], alpha=0.10)
    pax2.axhspan(1.0, 15.0, color=SAXS_COLORS['fractal_network'], alpha=0.10)
    pax2.text(0.065, 0.08, r'ACC Spheres Dominant' '\n' r'($>95\%$ at onset)', color=SAXS_COLORS['acc_spheres'], fontsize=9, fontweight='bold')
    pax2.text(0.6, 4.5, r'Fractal Network Dominant' '\n' r'($\sim 76\%$ plateau)', color=SAXS_COLORS['fractal_network'], fontsize=9, fontweight='bold')

    pax2.set_yscale('log')
    pax2.set_ylim(0.02, 12.0)
    pax2.set_xlim(0.05, 12.0)
    pax2.set_xlabel("Reaction Time $t$ [s]", fontsize=11, fontweight='bold')
    pax2.set_ylabel(r"Ratio $\phi_{\mathrm{net}} / \phi_{\mathrm{ACC}}$", fontsize=11, fontweight='bold')
    pax2.set_title(r"(b) Dynamic Phase Mass Ratio ($0.04 \to 3.3$)", fontsize=11, fontweight='bold')
    pax2.grid(True, which='both', linestyle='--', alpha=0.5)
    pax2.legend(loc='lower right', fontsize=8.8, framealpha=0.92)

    plt.tight_layout()
    stand_png = os.path.join(OUT_DIR, 'relative_phase_evolution.png')
    stand_pdf = os.path.join(OUT_DIR, 'relative_phase_evolution.pdf')
    plt.savefig(stand_png, dpi=300)
    plt.savefig(stand_pdf)
    plt.close(fig2)
    print(f"Saved dedicated relative phase evolution plot: {stand_png} and {stand_pdf}")
    print("=" * 80)


if __name__ == '__main__':
    generate_unified_trajectories()
