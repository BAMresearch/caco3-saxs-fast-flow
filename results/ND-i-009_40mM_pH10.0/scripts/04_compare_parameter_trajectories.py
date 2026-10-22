#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 04: Unified Multi-Scale Kinetic Trajectories (ND-i-009).
--------------------------------------------------------------
Plots cross-scale kinetic coupling as a function of real reaction time [s]:
1. Hierarchy of structural length scales across 4 decades (0.5 nm to 3.0 um).
2. Network fractal dimension D_f(t) consistency between short and long distance.
3. Coupling of local viscoelastic fluctuations A_fluc(t) to droplet emergence phi_large(t).

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

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)


def generate_unified_multiscale_trajectories():
    print("\n--- Generating Unified Multi-Scale Trajectories vs Real Time ---")
    out_dir = os.path.join(PROJECT_DIR, 'output', 'unified_trajectories')
    os.makedirs(out_dir, exist_ok=True)

    visco_path = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'viscoelastic', 'parameters.csv')
    long_path = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'option_c_forward', 'parameters.csv')

    df_sd = pd.read_csv(visco_path)
    df_ld = pd.read_csv(long_path)

    t_sd = df_sd['time_seconds'].values
    t_ld = df_ld['time_seconds'].values

    fig = plt.figure(figsize=(14, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], hspace=0.28, wspace=0.24)
    fig.suptitle("ND-i-009 Unified Multi-Scale Kinetics: Real-Time Viscoelastic Phase Separation", fontsize=14, fontweight='bold', y=0.98)

    # Panel 1: Hierarchy of Characteristic Length Scales Across 4 Orders of Magnitude
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(t_sd, df_sd['xi_fluc'], color='#ef4444', linewidth=2.0, label=r'Mesh Fluctuation Correlation: $\xi_{\rm fluc}(t) \sim 2.7 \to 2.9$ nm')
    if 'std_err_xi_fluc' in df_sd.columns:
        ax1.fill_between(t_sd, df_sd['xi_fluc'] - df_sd['std_err_xi_fluc'], df_sd['xi_fluc'] + df_sd['std_err_xi_fluc'], color='#ef4444', alpha=0.25)
    
    ax1.plot(t_sd, df_sd['Rmean'], color='#10b981', linewidth=2.2, label=r'Primary Building Block: $R_{\rm prim} \approx 8.3$ nm')
    if 'std_err_Rmean' in df_sd.columns:
        ax1.fill_between(t_sd, df_sd['Rmean'] - df_sd['std_err_Rmean'], df_sd['Rmean'] + df_sd['std_err_Rmean'], color='#10b981', alpha=0.25)
        
    ax1.plot(t_sd[0], df_sd['xi'].iloc[0], 'v', color='#f59e0b', markersize=8, label=r'Pre-Percolation Domain: Frame 0 $\xi = 46.7$ nm')
    ax1.axhline(1e4, color='#f59e0b', linestyle='--', linewidth=1.8, label=r'Percolated Network: $\xi \to \infty$')
    
    # Large spheres active for frames >= 1
    ld_active = df_ld['has_large_spheres'].astype(bool)
    ax1.plot(t_ld[ld_active], df_ld.loc[ld_active, 'mean_radius_large'], color='#2563eb', linewidth=2.2, label=r'Macroscopic Droplets: $R_{\rm large} \sim 2.7 \to 2.9\ \mu$m')
    if 'std_err_mean_radius_large' in df_ld.columns:
        ax1.fill_between(t_ld[ld_active],
                         df_ld.loc[ld_active, 'mean_radius_large'] - df_ld.loc[ld_active, 'std_err_mean_radius_large'],
                         df_ld.loc[ld_active, 'mean_radius_large'] + df_ld.loc[ld_active, 'std_err_mean_radius_large'],
                         color='#2563eb', alpha=0.25)

    ax1.set_yscale('log')
    ax1.set_ylim(0.5, 12000)
    ax1.set_ylabel("Length Scale [nm]", fontsize=11)
    ax1.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax1.set_title(r"(A) Hierarchy of Structural Length Scales Across Real Time ($0.045 - 5.30$ s) $\pm 1\sigma$", fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend(loc='center right', fontsize=9.5, framealpha=0.92)

    # Panel 2: Fractal Network Dimension D_f(t) Consistency
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(t_sd, df_sd['D'], color='#dc2626', linewidth=2.0, label='Short Distance ($3.0$ m, local mesh)')
    if 'std_err_D' in df_sd.columns:
        ax2.fill_between(t_sd, df_sd['D'] - df_sd['std_err_D'], df_sd['D'] + df_sd['std_err_D'], color='#dc2626', alpha=0.25, label=r'SD $\pm 1\sigma$ uncertainty')
        
    ax2.plot(t_ld, df_ld['D'], color='#2563eb', linewidth=2.0, label=r'Long Distance ($31.0$ m, mesoscale mesh)')
    if 'std_err_D' in df_ld.columns:
        ax2.fill_between(t_ld, df_ld['D'] - df_ld['std_err_D'], df_ld['D'] + df_ld['std_err_D'], color='#2563eb', alpha=0.20, label=r'LD $\pm 1\sigma$ uncertainty')
        
    ax2.axhline(1.90, color='#64748b', linestyle=':', label=r'Universal Dimension Reference ($D_f \approx 1.90$)')
    ax2.set_ylabel(r"Mass Fractal Dimension $D_f$", fontsize=11)
    ax2.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax2.set_ylim(1.2, 2.3)
    ax2.set_title(r"(B) Cross-Scale Network Dimension Agreement ($D_f \pm 1\sigma$)", fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=8.5)

    # Panel 3: Kinetic Coupling of Droplet Emergence to VPS Fluctuations
    ax3 = fig.add_subplot(gs[1, 1])
    ax3_twin = ax3.twinx()
    l1 = ax3.plot(t_sd, df_sd['A_fluc'], color='#dc2626', linewidth=2.0, label=r'Local VPS Fluctuation ($A_{\rm fluc}$)')
    if 'std_err_A_fluc' in df_sd.columns:
        ax3.fill_between(t_sd, df_sd['A_fluc'] - df_sd['std_err_A_fluc'], df_sd['A_fluc'] + df_sd['std_err_A_fluc'], color='#dc2626', alpha=0.22)
        
    l2 = ax3_twin.plot(t_ld, df_ld['total_volume_fraction_large'] * 1e3, color='#2563eb', linewidth=2.0, label=r'Droplet Fraction ($\phi_{\rm large} \times 10^3$)')
    if 'std_err_total_volume_fraction_large' in df_ld.columns:
        ax3_twin.fill_between(t_ld,
                              (df_ld['total_volume_fraction_large'] - df_ld['std_err_total_volume_fraction_large']) * 1e3,
                              (df_ld['total_volume_fraction_large'] + df_ld['std_err_total_volume_fraction_large']) * 1e3,
                              color='#2563eb', alpha=0.22)
                              
    ax3.set_ylabel(r"Fluctuation Amplitude $A_{\rm fluc}$ [cm$^{-1}$]", color='#dc2626', fontsize=11)
    ax3_twin.set_ylabel(r"Droplet Volume Fraction $\phi_{\rm large} \times 10^3$", color='#2563eb', fontsize=11)
    ax3.set_xlabel("Reaction Time $t$ [s]", fontsize=11)
    ax3.set_title(r"(C) Kinetic Coupling: Fluctuation vs Droplet Emergence $\pm 1\sigma$", fontsize=12, fontweight='bold')
    ax3.grid(True, linestyle='--', alpha=0.5)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='center left', fontsize=9.0)

    plt.tight_layout()
    unified_plot = os.path.join(out_dir, 'unified_cross_scale_trajectories.png')
    plt.savefig(unified_plot, dpi=160)
    plt.close(fig)
    print(f"Saved unified multi-scale trajectory: {unified_plot}")


def main():
    generate_unified_multiscale_trajectories()


if __name__ == '__main__':
    main()
