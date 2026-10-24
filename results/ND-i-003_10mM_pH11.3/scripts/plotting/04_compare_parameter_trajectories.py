#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 04: Unified Multi-Scale Kinetic Trajectories for ND-i-003.
----------------------------------------------------------------
Visualizes the physical mechanisms and kinetics of discrete ACC nanosphere
formation in the low supersaturation regime (10 mM CaCO3, pH 11.3):

1. Panel A: Evolution of Nanosphere Size & Polydispersity across Real Time (t = 0.060 to 10.877 s)
   - Particle mean radius <R>(t) growing from ~47 nm to ~88 nm (diameter 94 - 176 nm)
   - Size dispersion sigma_R(t) ~ 18 to 28 nm (~32% polydispersity)
   - Shaded induction regime (t < 2.45 s)
2. Panel B: Volume Fraction phi(t) vs Stoichiometric Precursor Limit
   - Monotonic growth of ACC volume fraction reaching phi = (4.71 +/- 0.08) x 10^-4
   - Direct comparison with theoretical maximum phi_max = 5.0 x 10^-4 for 10 mM CaCO3
3. Panel C: Cross-Scale Dynamic Range Coupling
   - Forward scattering intensity I(0) at 31 m surging by >4 orders of magnitude
   - Short distance high-q background remaining completely flat (I_bkg ~ 7-10 x 10^-5 cm^-1)
   - Proof of homogeneous nucleation of discrete ACC without spinodal gel network

Output saved to:
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

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)


def generate_unified_trajectories():
    print("=" * 80)
    print("  GENERATING UNIFIED MULTI-SCALE KINETIC TRAJECTORIES (ND-i-003)")
    print("=" * 80)

    out_dir = os.path.join(PROJECT_DIR, 'output', 'unified_trajectories')
    os.makedirs(out_dir, exist_ok=True)

    short_path = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    long_path = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')

    df_sd = pd.read_csv(short_path)
    df_ld = pd.read_csv(long_path)

    t_sd = df_sd['time_seconds'].values
    t_ld = df_ld['time_seconds'].values

    # Induction threshold
    mask_grow = ~df_ld['is_induction']
    t_onset = df_ld.loc[mask_grow, 'time_seconds'].min() if np.any(mask_grow) else 2.45

    fig = plt.figure(figsize=(14, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], hspace=0.28, wspace=0.24)
    fig.suptitle(
        "ND-i-003 Unified Multi-Scale Kinetics (10 mM CaCO3, pH 11.3)\n"
        "Homogeneous Nucleation & Growth of Discrete Polydisperse ACC Nanospheres",
        fontsize=13, fontweight='bold', y=0.98
    )

    # -------------------------------------------------------------
    # Panel A: Hierarchy of Characteristic Length Scales Across Time
    # -------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, :])
    ax1.axvspan(0, t_onset, color='#f1f5f9', alpha=0.9, label=f'Induction Period (t < {t_onset:.2f} s; Not Fitted)')

    # Schulz-Flory Mean Radius & Dispersion
    r_sch = df_ld['Rmean_schulz'].values if 'Rmean_schulz' in df_ld.columns else df_ld['Rmean'].values
    err_r_sch = df_ld['std_err_Rmean_schulz'].values if 'std_err_Rmean_schulz' in df_ld.columns else np.zeros_like(r_sch)
    sig_sch = df_ld['stdDev_schulz'].values if 'stdDev_schulz' in df_ld.columns else df_ld['stdDev'].values

    # Lognormal Mean Radius & Dispersion
    r_log = df_ld['Rmean_lognormal'].values if 'Rmean_lognormal' in df_ld.columns else r_sch
    err_r_log = df_ld['std_err_Rmean_lognormal'].values if 'std_err_Rmean_lognormal' in df_ld.columns else np.zeros_like(r_log)
    sig_log = df_ld['stdDev_lognormal'].values if 'stdDev_lognormal' in df_ld.columns else sig_sch

    r_min_s, r_max_s = r_sch[mask_grow].min(), r_sch[mask_grow].max()
    r_min_l, r_max_l = r_log[mask_grow].min(), r_log[mask_grow].max()
    sig_min_s, sig_max_s = sig_sch[mask_grow].min(), sig_sch[mask_grow].max()

    ax1.plot(t_ld[mask_grow], r_sch[mask_grow], 'o-', color='#dc2626', linewidth=2.2, markersize=5.5,
             label=f'Schulz-Flory $\\langle R \\rangle \\sim {r_min_s:.1f} \\to {r_max_s:.1f}$ nm')
    ax1.fill_between(t_ld[mask_grow],
                     r_sch[mask_grow] - err_r_sch[mask_grow],
                     r_sch[mask_grow] + err_r_sch[mask_grow],
                     color='#dc2626', alpha=0.20)

    ax1.plot(t_ld[mask_grow], r_log[mask_grow], 's--', color='#2563eb', linewidth=1.8, markersize=5.0,
             label=f'Lognormal $\\langle R \\rangle \\sim {r_min_l:.1f} \\to {r_max_l:.1f}$ nm')
    ax1.fill_between(t_ld[mask_grow],
                     r_log[mask_grow] - err_r_log[mask_grow],
                     r_log[mask_grow] + err_r_log[mask_grow],
                     color='#2563eb', alpha=0.15)

    # Size Dispersion sigma_R
    ax1.plot(t_ld[mask_grow], sig_sch[mask_grow], '^:', color='#f97316', linewidth=1.7, markersize=4.5,
             label=f'Size Dispersion $\\sigma_R \\sim {sig_min_s:.1f} \\to {sig_max_s:.1f}$ nm')

    # Short Distance Limit (absence of primary nodes)
    ax1.axhline(1.0, color='#64748b', linestyle=':', linewidth=1.5,
                label=r'Short-Distance SAXS Limit: No detectable nodes / PNCs ($R < 1$ nm)')

    ax1.set_yscale('log')
    ax1.set_ylim(0.5, 300)
    ax1.set_ylabel("Characteristic Radius [nm]", fontsize=11)
    ax1.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax1.set_title(r"(A) Evolution of Discrete ACC Nanosphere Size: Schulz-Flory vs Lognormal Comparison", fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.45)
    ax1.legend(loc='lower right', fontsize=8.8, framealpha=0.92)

    # -------------------------------------------------------------
    # Panel B: Volume Fraction Growth vs Stoichiometric Limit
    # -------------------------------------------------------------
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.axvspan(0, t_onset, color='#f1f5f9', alpha=0.9)

    phi_ld = df_ld['phi'].values
    err_phi = df_ld['std_err_phi'].values if 'std_err_phi' in df_ld.columns else np.zeros_like(phi_ld)

    ax2.plot(t_ld[mask_grow], phi_ld[mask_grow], '^-', color='#16a34a', linewidth=2.2, markersize=5,
             label=r'Fitted ACC Volume Fraction $\phi(t)$')
    ax2.fill_between(t_ld[mask_grow],
                     np.maximum(phi_ld[mask_grow] - err_phi[mask_grow], 1e-6),
                     phi_ld[mask_grow] + err_phi[mask_grow],
                     color='#16a34a', alpha=0.22, label=r'$\phi \pm 1\sigma$ Uncertainty')

    # Stoichiometric maximum for 10 mM CaCO3:
    phi_max = 5.0e-4
    ax2.axhline(phi_max, color='#dc2626', linestyle='--', linewidth=1.8,
                label=r'Stoichiometric Maximum $\phi_{\rm max} \approx 5.0 \times 10^{-4}$ (10 mM)')

    ax2.set_yscale('log')
    ax2.set_ylim(2e-6, 1e-3)
    ax2.set_ylabel(r"ACC Volume Fraction $\phi$", fontsize=11)
    ax2.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax2.set_title(r"(B) ACC Nanosphere Precipitation Kinetics vs Stoichiometry", fontsize=12, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', alpha=0.45)
    ax2.legend(loc='lower right', fontsize=8.5, framealpha=0.92)

    # -------------------------------------------------------------
    # Panel C: Cross-Scale Dynamic Range Coupling (31 m vs 3 m)
    # -------------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axvspan(0, t_onset, color='#f1f5f9', alpha=0.9)

    # Compute forward scattering I(0) for long distance
    i0_ld = []
    from models.scattering_models import calc_sphere_intensity
    for _, r in df_ld.iterrows():
        if not r.get('is_induction', False) and pd.notna(r['Rmean']) and float(r.get('phi', 0.0)) > 0:
            i0_val = calc_sphere_intensity(np.array([0.002]), float(r['Rmean']), float(r['stdDev']), float(r['phi']), bkg=0.0, rho2=1.0)[0]
        else:
            i0_val = 1e-4
        i0_ld.append(i0_val)
    i0_ld = np.array(i0_ld)

    ax3_twin = ax3.twinx()

    l1 = ax3.plot(t_ld[mask_grow], i0_ld[mask_grow], 'o-', color='#2563eb', linewidth=2.0, markersize=4.5,
                  label=r'ACC Forward Scattering $I_0$ [31.0 m]')
    l2 = ax3.plot(t_ld, df_ld['A_porod'] * 1e8, 'D:', color='#16a34a', linewidth=1.8, markersize=4.0,
                  label=r'Droplet Porod Equiv. $A_{\rm porod} \times 10^8$ [31.0 m]')
    l3 = ax3_twin.plot(t_sd, df_sd['bkg'] * 1e4, 's--', color='#9333ea', linewidth=1.5, markersize=4,
                       label=r'Solvent Bkg $I_{\rm bkg} \times 10^4$ [3.0 m]')

    ax3.set_yscale('log')
    ax3.set_ylim(1e-2, 1000)
    ax3_twin.set_ylim(0.0, 2.5)

    ax3.set_ylabel(r"Intensity & Amplitude [$A_p \times 10^8$ / $I_0$]", color='#2563eb', fontsize=10.5)
    ax3_twin.set_ylabel(r"High-$q$ Bkg $I_{\rm bkg} \times 10^4$ [$\mathrm{cm}^{-1}$]", color='#9333ea', fontsize=10.5)
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax3.set_title(r"(C) Cross-Scale Divergence: ACC & Droplets vs Flat Water", fontsize=12, fontweight='bold')
    ax3.grid(True, which='both', linestyle='--', alpha=0.45)

    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='upper right', fontsize=8.2, framealpha=0.92)

    plt.tight_layout()
    unified_plot = os.path.join(out_dir, 'unified_cross_scale_trajectories.png')
    plt.savefig(unified_plot, dpi=160)
    plt.close(fig)
    print(f"Saved unified multi-scale trajectory: {unified_plot}")


def main():
    generate_unified_trajectories()


if __name__ == '__main__':
    main()
