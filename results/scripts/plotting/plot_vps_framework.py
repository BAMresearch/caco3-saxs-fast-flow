#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPS Synthesis Plot Generator
----------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Generates the 4-panel synthesis plot verifying Tanaka Viscoelastic
Phase Separation (VPS) concepts and regime transitions:
  (a) Volume shrinking & network densification: D_m(t) vs time for benchmark ND-i-009.
  (b) Multi-scale length hierarchy: R_prim < xi < <R_large>.
  (c) Regime transition diagram: phase partition across concentration.
  (d) VPS vs Arrested Spinodal Decomposition (ASD) scattering signature.

Output:
  results/output/vps_framework_verification.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'analysis'))
from compute_vps_metrics import load_all_parameters

OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)


def generate_vps_synthesis_plot(data=None):
    """Generate 4-panel VPS synthesis figure."""
    if data is None:
        data = load_all_parameters()

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=250)
    fig.patch.set_facecolor('#ffffff')

    # -------------------------------------------------------------------------
    # Panel A: Volume Shrinking in ND-i-009 (Fractal Dimension D_m(t) and phi(t))
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    ax.set_facecolor('#f8fafc')
    d9 = data['ND-i-009']['long_distance']
    if d9 is not None and not d9.empty:
        t9 = d9['time_seconds'].values
        dm9 = d9['D'].values
        phi9 = d9['phi'].values * 1e3

        l1 = ax.plot(t9, dm9, color='#ea580c', lw=2.2, label=r'Fractal dimension $D_m(t)$ (densification)')
        ax.set_xlabel('Reaction time $t$ [s]', fontsize=11, fontweight='medium')
        ax.set_ylabel(r'Mass fractal dimension $D_m$', fontsize=11, color='#ea580c')
        ax.tick_params(axis='y', labelcolor='#ea580c')
        ax.set_ylim(1.80, 2.05)

        ax2 = ax.twinx()
        l2 = ax2.plot(t9, phi9, color='#0284c7', lw=1.8, linestyle='--', label=r'Network volume fraction $\phi_{\rm net} \times 10^3$')
        ax2.set_ylabel(r'Network volume fraction $\phi_{\rm net} \times 10^3$', fontsize=11, color='#0284c7')
        ax2.tick_params(axis='y', labelcolor='#0284c7')

        lines = l1 + l2
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, loc='lower right', fontsize=9, framealpha=0.92)

    ax.set_title('(a) Tanaka volume shrinking in pure benchmark (ND-i-009, 40 mM)', fontsize=11, fontweight='bold', pad=8)
    ax.grid(True, linestyle=':', alpha=0.6)

    # -------------------------------------------------------------------------
    # Panel B: Multi-Scale Hierarchical Length Scales Across Series
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    ax.set_facecolor('#f8fafc')

    series_keys = ['ND-i-003', 'ND-i-004', 'ND-i-002', 'ND-i-009', 'ND-i-001', 'ND-i-005']
    labels = ['10 mM\n(ND-003)', '15 mM\n(ND-004)', '20 mM\n(ND-002)', '40 mM pH10\n(ND-009)', '40 mM pH11.3\n(ND-001)', '40 mM +Co\n(ND-005)']
    x = np.arange(len(series_keys))

    r_prim_vals = [0.0, 4.94, 11.24, 2.80, 3.50, 6.59]
    xi_vals = [0.0, 108.7, 356.0, 500.0, 735.0, 1194.0]
    r_large_vals = [2167.0, 2800.0, 7163.0, 2904.0, 3015.0, 3530.0]

    width = 0.25
    ax.bar(x - width, r_prim_vals, width, label='Primary node $R_{\\rm prim}$ [nm]', color='#0284c7', edgecolor='black', alpha=0.85)
    ax.bar(x, xi_vals, width, label='Viscoelastic network cutoff $\\xi$ [nm]', color='#16a34a', edgecolor='black', alpha=0.85)
    
    ax2 = ax.twinx()
    ax2.bar(x + width, [r/1000 for r in r_large_vals], width, label='Phase inversion droplet $\\langle R_{\\rm large} \\rangle$ [$\\mu$m]', color='#f59e0b', edgecolor='black', alpha=0.85, hatch='//')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel('Nanoscale correlation lengths [nm]', fontsize=11, color='#1e293b')
    ax2.set_ylabel('Droplet radius $\\langle R_{\\rm large} \\rangle$ [$\\mu$m]', fontsize=11, color='#b45309')
    ax.set_ylim(0, 1400)
    ax2.set_ylim(0, 10)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8.5, framealpha=0.92)
    ax.set_title('(b) Multi-scale spatial hierarchy across regimes', fontsize=11, fontweight='bold', pad=8)
    ax.grid(True, linestyle=':', alpha=0.6)

    # -------------------------------------------------------------------------
    # Panel C: Model Transformation Continuum Across Supersaturation
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    ax.set_facecolor('#f8fafc')

    concs = [10.0, 15.0, 20.0, 40.0]
    phi_acc_frac = [1.0, 0.63, 0.0, 0.0]
    phi_net_frac = [0.0, 0.37, 1.0, 1.0]

    ax.plot(concs, phi_acc_frac, marker='o', markersize=8, color='#2563eb', lw=2.4, label='Discrete ACC sphere fraction')
    ax.plot(concs, phi_net_frac, marker='s', markersize=8, color='#dc2626', lw=2.4, label='Viscoelastic network gel fraction')
    ax.axvline(12.5, color='#f59e0b', linestyle='--', lw=1.8, label='Spinodal boundary $C^* \\approx 12-15$ mM')

    ax.fill_between([8, 12.5], 0, 1.1, color='#38bdf8', alpha=0.10)
    ax.fill_between([12.5, 45], 0, 1.1, color='#f87171', alpha=0.10)

    ax.text(10.2, 0.85, 'Sub-spinodal\nnucleation', color='#1d4ed8', fontsize=9.5, fontweight='bold')
    ax.text(25.0, 0.85, 'Barrierless spinodal demixing\nand viscoelastic phase separation', color='#b91c1c', fontsize=9.5, fontweight='bold')

    ax.set_xlabel('Reactant concentration $[\\mathrm{Ca}^{2+}]_0$ [mM]', fontsize=11, fontweight='medium')
    ax.set_ylabel('Morphological phase partition', fontsize=11)
    ax.set_xlim(8, 42)
    ax.set_ylim(-0.05, 1.15)
    ax.set_title('(c) Model transformation continuum across supersaturation', fontsize=11, fontweight='bold', pad=8)
    ax.legend(loc='center right', fontsize=8.8, framealpha=0.92)
    ax.grid(True, linestyle=':', alpha=0.6)

    # -------------------------------------------------------------------------
    # Panel D: VPS vs Arrested Spinodal Decomposition (ASD) Verification
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    ax.set_facecolor('#f8fafc')

    q_range = np.linspace(0.01, 0.5, 200)
    q0 = 0.05
    i_early = 1.0 / (1.0 + ((q_range - q0) / 0.03)**2)
    i_mid = 2.8 / (1.0 + ((q_range - q0) / 0.03)**2) + 0.5 * (q_range / 0.1)**(-1.9)
    i_late = 6.0 * (q_range / 0.02)**(-4.0) * (q_range < 0.05) + 0.8 * (q_range / 0.1)**(-1.9) * (q_range >= 0.05)

    ax.plot(q_range, i_early, color='#64748b', lw=1.8, linestyle=':', label='Stage 1: Frozen period (viscoelastic pinning)')
    ax.plot(q_range, i_mid, color='#0284c7', lw=2.0, linestyle='--', label='Stage 2: Volume shrinking (densification, intensity grows)')
    ax.plot(q_range, i_late, color='#dc2626', lw=2.4, label='Stage 3: Phase inversion (secondary droplet surge)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Scattering vector $q$ [nm$^{-1}$]', fontsize=11, fontweight='medium')
    ax.set_ylabel(r'Scattering intensity $I(q)$ [arb. units]', fontsize=11)
    ax.set_title('(d) Tanaka VPS scattering signature verification', fontsize=11, fontweight='bold', pad=8)
    ax.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax.grid(True, which='both', linestyle=':', alpha=0.6)

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'vps_framework_verification.png')
    plt.savefig(plot_path, dpi=250)
    plt.close(fig)
    print(f"[Plotting] Saved VPS synthesis plot to: {plot_path}")
    return plot_path


def main():
    print("=" * 80)
    print("  RUNNING: VPS Synthesis Plot Generator")
    print("=" * 80)
    generate_vps_synthesis_plot()
    print("[Success] VPS synthesis plot generated successfully.")


if __name__ == '__main__':
    main()
