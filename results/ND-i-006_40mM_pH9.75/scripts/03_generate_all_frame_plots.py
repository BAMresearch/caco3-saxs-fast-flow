#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 03: Parallel Rendering of 100 Comprehensive Diagnostic Frame Plots for ND-i-006.
--------------------------------------------------------------------------------------
Renders high-resolution diagnostic plots for each of the 50 short-distance (3.0 m)
and 50 long-distance (31.0 m) frames:
- Log-log intensity scattering pattern with multi-scale physical model decomposition:
  * Mesoscale Spinodal Fractal Network (Primary Building Blocks R_prim, Dm, xi)
  * Secondary Large Droplets / Spheres (R_large, phi_large)
  * Macroscopic Porod Surface Tail (A_porod / q^4)
  * Incoherent Flat Background (bkg)
  * NO ACC Nanospheres (experimentally absent)
- Normalized point-by-point residual distribution in units of experimental standard deviation (sigma)

Outputs:
  output/short_distance/fits/frame_00.png to frame_49.png (50 files)
  output/long_distance/fits/frame_00.png to frame_49.png (50 files)
"""

import os
import sys
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    RHO2_CALIBRATED
)


def plot_single_short_frame(args):
    f_idx, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})
    
    valid = (y_frame > 0) & (dy_frame > 0) & (q_data >= 0.022)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = float(row_dict.get('time_seconds', 0.0))
    redchi = float(row_dict.get('chi2_red', 0.0))
    A_porod = float(row_dict.get('A_porod', 0.0))
    phi_net = float(row_dict.get('phi_net', 1e-5))
    R_prim = float(row_dict.get('R_prim', 1.50))
    std_prim = float(row_dict.get('std_prim', 3.03))
    Dm = float(row_dict.get('Dm', 1.80))
    xi = float(row_dict.get('xi', 30000.0))
    bkg = float(row_dict.get('bkg', 2.3e-4))

    q_dense = np.logspace(np.log10(0.02), np.log10(2.6), 250)
    I_net = calc_primary_fractal_intensity(q_dense, R_prim, std_prim, phi_net, Dm, xi, rho2=RHO2_CALIBRATED, integration_steps=60)
    I_porod = calc_porod_intensity(q_dense, A_porod, 4.0) if A_porod > 0 else np.zeros_like(q_dense)
    I_bkg = np.full_like(q_dense, bkg)
    I_tot = I_net + I_porod + I_bkg

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full Model ($\\chi^2_{{\\mathrm{{red}}}} = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_net, color='#16a34a', linestyle='--', linewidth=1.5,
                 label=f'Spinodal Network ($R_{{\\rm prim}} = {R_prim:.1f}$ nm, $D_m = {Dm:.2f}$, $\\phi = {phi_net:.2e}$)')
    if A_porod > 0:
        ax_main.plot(q_dense, I_porod, color='#d97706', linestyle=':', linewidth=1.3,
                     label=f'Porod Surface Tail ($A = {A_porod*1e6:.2f}\\times 10^{{-6}}$)')
    ax_main.axhline(bkg, color='#94a3b8', linestyle=':', linewidth=1.0, label=f'Background ({bkg:.2e} cm$^{{-1}}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.02, 2.65)
    ax_main.set_ylim(max(1e-4, y_v.min() * 0.3), max(1.0, y_v.max() * 3.0))
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_main.set_title(f"ND-i-006 Short Distance Frame {f_idx:02d} (t = {t_s:.3f} s): Mass Fractal Network", fontsize=11, fontweight='bold')
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='upper right', fontsize=8.0, framealpha=0.92)

    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.02, 2.65)
    ax_res.set_ylim(-4.5, 4.5)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_res.set_ylabel("Res. [$\\sigma$]", fontsize=9.5, fontweight='bold')
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    plt.close(fig)


def plot_single_long_frame(args):
    f_idx, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})

    valid = (y_frame > 0) & (dy_frame > 0) & (q_data >= 0.0025) & (q_data <= 0.22)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = float(row_dict.get('time_seconds', 0.0))
    redchi = float(row_dict.get('chi2_red', 0.0))
    Dm = float(row_dict.get('Dm', 1.80))
    xi = float(row_dict.get('xi', 2500.0))
    phi_net = float(row_dict.get('phi_net', 2.5e-5))
    R_prim = float(row_dict.get('R_prim', 1.50))
    std_prim = float(row_dict.get('std_prim', 3.03))
    R_large = float(row_dict.get('mean_radius_large', 3113.0))
    std_large = float(row_dict.get('std_radius_large', 267.0))
    phi_large = float(row_dict.get('total_volume_fraction_large', 1.1e-4))
    A_porod = float(row_dict.get('A_porod', 0.0))
    bkg = float(row_dict.get('bkg', 1e-4))

    q_dense = np.logspace(np.log10(0.002), np.log10(0.25), 250)
    I_frac = calc_primary_fractal_intensity(q_dense, R_prim, std_prim, phi_net, Dm, xi, rho2=RHO2_CALIBRATED, integration_steps=60)
    I_large = calc_large_sphere_intensity(q_dense, R_large, std_large, phi_large, rho2=RHO2_CALIBRATED, integration_steps=60)
    I_porod = calc_porod_intensity(q_dense, A_porod, 4.0) if A_porod > 0 else np.zeros_like(q_dense)
    I_bkg = np.full_like(q_dense, bkg)
    I_tot = I_frac + I_large + I_porod + I_bkg

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full Model ($\\chi^2_{{\\mathrm{{red}}}} = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_large, color='#0284c7', linestyle='-.', linewidth=1.5,
                 label=f'Large Spheres ($\\langle R \\rangle = {R_large:.0f}$ nm, $\\phi = {phi_large:.2e}$)')
    ax_main.plot(q_dense, I_frac, color='#16a34a', linestyle='--', linewidth=1.5,
                 label=f'Fractal Network ($D_m = {Dm:.2f}$, $\\xi = {xi:.0f}$ nm, $\\phi = {phi_net:.2e}$)')
    if A_porod > 0:
        ax_main.plot(q_dense, I_porod, color='#d97706', linestyle=':', linewidth=1.3,
                     label=f'Porod Surface Tail ($A = {A_porod*1e6:.2f}\\times 10^{{-6}}$)')
    ax_main.axhline(bkg, color='#94a3b8', linestyle=':', linewidth=1.0, label=f'Background ({bkg:.2e} cm$^{{-1}}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.002, 0.25)
    ax_main.set_ylim(max(1e-3, y_v.min() * 0.3), max(100.0, y_v.max() * 3.0))
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_main.set_title(f"ND-i-006 Long Distance Frame {f_idx:02d} (t = {t_s:.3f} s): Network & Large Droplets", fontsize=11, fontweight='bold')
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='upper right', fontsize=8.0, framealpha=0.92)

    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.002, 0.25)
    ax_res.set_ylim(-4.5, 4.5)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_res.set_ylabel("Res. [$\\sigma$]", fontsize=9.5, fontweight='bold')
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    plt.close(fig)


def generate_all_diagnostic_plots():
    print("=" * 80)
    print("  PARALLEL DIAGNOSTIC PLOT GENERATION FOR ND-i-006 (100 PLOTS)")
    print("  Model: Spinodal Mass Fractal Network + Large Droplets (No ACC Particles)")
    print("=" * 80)

    # 1. Short Distance Plots
    out_s_dir = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fits')
    os.makedirs(out_s_dir, exist_ok=True)
    df_s_data = pd.read_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_short_distance.csv'))
    q_s = df_s_data['q'].values
    df_s_params = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv'))
    df_s_curves = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fitted_curves.csv'))

    s_tasks = []
    for _, row in df_s_params.iterrows():
        f = int(row['frame_index'])
        y_val = df_s_data[str(f)].values
        dy_val = np.maximum(df_s_data['error_' + str(f)].values, 1e-4)
        fit_c = df_s_curves[str(f)].values
        out_f = os.path.join(out_s_dir, f"frame_{f:02d}.png")
        s_tasks.append((f, q_s, y_val, dy_val, fit_c, row.to_dict(), out_f))

    # 2. Long Distance Plots
    out_l_dir = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fits')
    os.makedirs(out_l_dir, exist_ok=True)
    df_l_data = pd.read_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_long_distance.csv'))
    q_l = df_l_data['q'].values
    df_l_params = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv'))
    df_l_curves = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fitted_curves.csv'))

    l_tasks = []
    for _, row in df_l_params.iterrows():
        f = int(row['frame_index'])
        y_val = df_l_data[str(f)].values
        dy_val = np.maximum(df_l_data['error_' + str(f)].values, 1e-4)
        fit_c = df_l_curves[str(f)].values
        out_f = os.path.join(out_l_dir, f"frame_{f:02d}.png")
        l_tasks.append((f, q_l, y_val, dy_val, fit_c, row.to_dict(), out_f))

    n_workers = min(os.cpu_count() or 4, 8)
    print(f"Rendering 50 Short-Distance plots with {n_workers} workers...")
    with multiprocessing.Pool(n_workers) as pool:
        pool.map(plot_single_short_frame, s_tasks)
    print("Short-Distance plots successfully rendered!")

    print(f"Rendering 50 Long-Distance plots with {n_workers} workers...")
    with multiprocessing.Pool(n_workers) as pool:
        pool.map(plot_single_long_frame, l_tasks)
    print("Long-Distance plots successfully rendered!")
    print("=" * 80)


if __name__ == '__main__':
    generate_all_diagnostic_plots()
