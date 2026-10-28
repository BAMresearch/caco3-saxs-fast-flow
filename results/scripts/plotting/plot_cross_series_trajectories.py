#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-Series Trajectory Plot Generator
--------------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Renders publication-grade 4-panel synthesis figure comparing time-resolved
structural trajectories across all six calcium carbonate reaction series.

Panels:
  (a) Precipitation extent conversion alpha(t) vs time.
  (b) Meso-scale length scale: ACC radius vs network cutoff xi.
  (c) Structural bifurcation across concentration: spheres vs network.
  (d) Universal phase inversion droplet population across series.

Output:
  results/output/cross_series_trajectories_comparison.png
"""

import os
import sys
import numpy as np
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, RESULTS_DIR)
sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'analysis'))
from saxs_colors import SAXS_COLORS
from extract_cross_series_trajectories import load_series_trajectories

OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)


def generate_synthesis_figure(trajectories=None):
    """Generate 4-panel publication figure showing cross-series kinetic evolution."""
    if trajectories is None:
        trajectories = load_series_trajectories()

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
    ax.plot([15, 20, 40], [4.94, 11.20, 8.35], 'o-', color='#16a34a', linewidth=2.2, markersize=8, label=r'Primary building blocks $R_{\rm prim}$ (pure $\mathrm{CaCO}_3$)')
    ax.plot([20, 40], [7.63, 2.72], 'd:', color='#9333ea', linewidth=1.8, markersize=7, label=r'PNC fluctuation mesh $\xi_{\rm fluc}$ (pure $\mathrm{CaCO}_3$)')
    ax.plot([15, 20, 40], [108.7, 356.0, 735.0], '^--', color='#dc2626', linewidth=2.0, markersize=8, label=r'Network correlation cutoff $\xi$ (pure $\mathrm{CaCO}_3$)')

    # Additive Dopant annotations
    ax.plot(40.0, 7.03, 'v', color='#d97706', markersize=9, label=r'Co$^{2+}$ additive: node $R_{\rm prim}$ ($7.03$ nm)')
    ax.plot(40.0, 8.05, 'D', color='#d97706', markersize=8, label=r'Co$^{2+}$ additive: expanded mesh $\xi_{\rm fluc}$ ($8.05$ nm)')
    ax.plot(40.0, 1194.0, '^', color='#d97706', markersize=9, label=r'Co$^{2+}$ additive: expanded cutoff $\xi$ ($1,194$ nm)')

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
    ax_in.bar(x_in, acc_pct, w_in, color=SAXS_COLORS['acc_spheres'], label='ACC spheres')
    ax_in.bar(x_in, net_pct, w_in, bottom=acc_pct, color=SAXS_COLORS['fractal_network'], label='Fractal gel')
    ax_in.bar(x_in, pnc_pct, w_in, bottom=acc_pct + net_pct, color=SAXS_COLORS['pnc_fluctuations'], label='Tanaka PNC mesh')
    ax_in.bar(x_in, drp_pct, w_in, bottom=acc_pct + net_pct + pnc_pct, color=SAXS_COLORS['droplets_porod'], label='Droplets / Porod')
    
    ax_in.set_xticks(x_in)
    ax_in.set_xticklabels(['10', '15', '20', '40'], fontsize=7.5)
    ax_in.set_xlabel(r'$[\mathrm{Ca}^{2+}]_0$ [mM]', fontsize=7.5)
    ax_in.set_ylim(0, 105)
    ax_in.set_ylabel('Share [%]', fontsize=7.5)
    ax_in.set_title('Phase dominance', fontsize=8.0, fontweight='bold')
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
    print(f"[Plotting] Saved cross-series trajectory comparison plot to: {plot_path}")
    return plot_path


def main():
    print("=" * 80)
    print("  RUNNING: Cross-Series Trajectory Plot Generator")
    print("=" * 80)
    generate_synthesis_figure()
    print("[Success] Cross-series trajectories comparison plot generated successfully.")


if __name__ == '__main__':
    main()
