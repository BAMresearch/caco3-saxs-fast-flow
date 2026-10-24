#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 03: Parallel Rendering of 100 Comprehensive Diagnostic Frame Plots for ND-i-004.
--------------------------------------------------------------------------------------
Renders high-resolution diagnostic plots for each of the 50 short-distance (3.0 m)
and 50 long-distance (31.0 m) frames:
- Log-log intensity scattering pattern with multi-scale physical model decomposition:
  * Mesoscale Spinodal Fractal Network (Primary Building Blocks R_prim = 3.5 nm, Dm, xi)
  * ACC Nanospheres (Mean Radius R_acc, Polydispersity std_acc, Volume Fraction phi_acc)
  * Macroscopic Porod Surface Tail (A_porod / q^4)
  * Incoherent Flat Background (bkg)
- Normalized point-by-point residual distribution in units of experimental standard deviation (sigma)
- Legend positioned in the upper right quadrant to avoid obscuring experimental data points

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

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
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
    R_acc = float(row_dict.get('R_acc', 30.0))
    std_acc = float(row_dict.get('std_acc', 8.0))
    phi_acc = float(row_dict.get('phi_acc', 3e-5))
    A_porod = float(row_dict.get('A_porod', 0.0))
    phi_net = float(row_dict.get('phi_net', 1e-5))
    R_prim = float(row_dict.get('R_prim', 3.5))
    std_prim = float(row_dict.get('std_prim', 1.5))
    Dm = float(row_dict.get('Dm', 2.33))
    xi = float(row_dict.get('xi', 200.0))
    bkg = float(row_dict.get('bkg', 2e-4))

    q_dense = np.logspace(np.log10(0.02), np.log10(2.6), 250)
    I_acc = calc_large_sphere_intensity(q_dense, R_acc, std_acc, phi_acc, rho2=RHO2_CALIBRATED, integration_steps=80)
    I_net = calc_primary_fractal_intensity(q_dense, R_prim, std_prim, phi_net, Dm, xi, rho2=RHO2_CALIBRATED, integration_steps=80)
    I_porod = calc_porod_intensity(q_dense, A_porod, 4.0) if A_porod > 0 else np.zeros_like(q_dense)
    I_bkg = np.full_like(q_dense, bkg)
    I_tot = I_acc + I_net + I_porod + I_bkg

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full Model ($\\chi^2_{{\\mathrm{{red}}}} = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_acc, color='#0284c7', linestyle='-.', linewidth=1.5,
                 label=f'ACC Nanospheres ($\\langle R \\rangle = {R_acc:.1f} \\pm {std_acc:.1f}$ nm, $\\phi = {phi_acc:.2e}$)')
    ax_main.plot(q_dense, I_net, color='#16a34a', linestyle='--', linewidth=1.3,
                 label=f'Spinodal Network ($R_{{\\rm prim}} = {R_prim:.1f}$ nm, $\\phi = {phi_net:.2e}$)')
    if A_porod > 0:
        ax_main.plot(q_dense, I_porod, color='#d97706', linestyle=':', linewidth=1.3,
                     label=f'Porod Surface Tail ($A = {A_porod*1e6:.2f}\\times 10^{{-6}}$)')
    ax_main.axhline(bkg, color='#94a3b8', linestyle=':', linewidth=1.0, label=f'Background ({bkg:.2e} cm$^{{-1}}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.02, 2.65)
    ax_main.set_ylim(max(1e-4, y_v.min() * 0.3), max(1.0, y_v.max() * 3.0))
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_main.set_title(f"ND-i-004 Short Distance Frame {f_idx:02d} (t = {t_s:.3f} s): ACC Nanosphere Growth", fontsize=11, fontweight='bold')
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
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f_idx


def plot_single_long_frame(args):
    f_idx, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})
    
    valid = (y_frame > 0) & (dy_frame > 0) & (q_data >= 0.002) & (q_data <= 0.22)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = float(row_dict.get('time_seconds', 0.0))
    redchi = float(row_dict.get('chi2_red', 0.0))
    phi_net = float(row_dict.get('phi_net', 5e-5))
    Dm = float(row_dict.get('Dm', 2.33))
    xi = float(row_dict.get('xi', 200.0))
    A_porod = float(row_dict.get('A_porod', 0.0))
    R_acc = float(row_dict.get('R_acc_transferred', 40.0))
    std_acc = float(row_dict.get('std_acc_transferred', 7.5))
    phi_acc = float(row_dict.get('phi_acc', 3e-5))
    bkg = float(row_dict.get('bkg', 1e-4))

    q_dense = np.logspace(np.log10(0.0018), np.log10(0.25), 250)
    I_prim = calc_primary_fractal_intensity(q_dense, 3.5, 1.5, phi_net, Dm, xi, rho2=RHO2_CALIBRATED, integration_steps=80)
    I_acc = calc_large_sphere_intensity(q_dense, R_acc, std_acc, phi_acc, rho2=RHO2_CALIBRATED, integration_steps=80)
    I_porod = calc_porod_intensity(q_dense, A_porod, 4.0) if A_porod > 0 else np.zeros_like(q_dense)
    I_tot = I_prim + I_acc + I_porod + bkg

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full Model ($\\chi^2_{{\\mathrm{{red}}}} = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_prim, color='#16a34a', linestyle='--', linewidth=1.3,
                 label=f'Network Fractal ($D_m = {Dm:.2f}$, $\\xi = {xi:.0f}$ nm)')
    ax_main.plot(q_dense, I_acc, color='#0284c7', linestyle='-.', linewidth=1.4,
                 label=f'ACC Nanospheres ($\\langle R \\rangle = {R_acc:.1f}$ nm, $\\phi = {phi_acc:.2e}$)')
    if A_porod > 0:
        ax_main.plot(q_dense, I_porod, color='#d97706', linestyle=':', linewidth=1.3,
                     label=f'Porod Surface Tail ($A = {A_porod*1e6:.2f}\\times 10^{{-6}}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.0018, 0.25)
    ax_main.set_ylim(max(1e-3, y_v.min() * 0.3), y_v.max() * 2.5)
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_main.set_title(f"ND-i-004 Long Distance Frame {f_idx:02d} (t = {t_s:.3f} s): Hybrid Multiscale Model", fontsize=11, fontweight='bold')
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='lower left', fontsize=8.0, framealpha=0.92)

    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 0.25)
    ax_res.set_ylim(-4.5, 4.5)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=11, fontweight='bold')
    ax_res.set_ylabel("Res. [$\\sigma$]", fontsize=9.5, fontweight='bold')
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f_idx


def generate_all_diagnostic_frame_plots():
    print("=" * 80)
    print("  RENDERING 100 DIAGNOSTIC FRAME PLOTS FOR ND-i-004 (PARALLEL)")
    print("=" * 80)

    # 1. Setup Directories
    short_fits_dir = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fits')
    long_fits_dir = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fits')
    os.makedirs(short_fits_dir, exist_ok=True)
    os.makedirs(long_fits_dir, exist_ok=True)

    # 2. Load Short Distance Results
    short_data_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_short_distance.csv')
    short_params_csv = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    short_curves_csv = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fitted_curves.csv')

    df_short_data = pd.read_csv(short_data_csv)
    df_short_params = pd.read_csv(short_params_csv)
    df_short_curves = pd.read_csv(short_curves_csv)
    q_short = df_short_data['q'].values

    # 3. Load Long Distance Results
    long_data_csv = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_long_distance.csv')
    long_params_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')
    long_curves_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fitted_curves.csv')

    df_long_data = pd.read_csv(long_data_csv)
    df_long_params = pd.read_csv(long_params_csv)
    df_long_curves = pd.read_csv(long_curves_csv)
    q_long = df_long_data['q'].values

    # 4. Prepare Task Lists
    short_tasks = []
    for _, row in df_short_params.iterrows():
        f_idx = int(row['frame_index'])
        y_frame = df_short_data[str(f_idx)].values
        dy_frame = np.maximum(df_short_data['error_' + str(f_idx)].values, 1e-4)
        fit_curve = df_short_curves[str(f_idx)].values
        out_file = os.path.join(short_fits_dir, f"frame_{f_idx:02d}.png")
        short_tasks.append((f_idx, q_short, y_frame, dy_frame, fit_curve, row.to_dict(), out_file))

    long_tasks = []
    for _, row in df_long_params.iterrows():
        f_idx = int(row['frame_index'])
        y_frame = df_long_data[str(f_idx)].values
        dy_frame = np.maximum(df_long_data['error_' + str(f_idx)].values, 1e-4)
        fit_curve = df_long_curves[str(f_idx)].values
        out_file = os.path.join(long_fits_dir, f"frame_{f_idx:02d}.png")
        long_tasks.append((f_idx, q_long, y_frame, dy_frame, fit_curve, row.to_dict(), out_file))

    n_cpus = max(1, multiprocessing.cpu_count() - 1)
    print(f"Parallel rendering utilizing {n_cpus} CPU workers...")

    # Render Short-Distance Frames
    print(f"\nRendering 50 Short-Distance diagnostic plots -> {short_fits_dir}")
    with multiprocessing.Pool(processes=n_cpus) as pool:
        for done_idx in pool.imap_unordered(plot_single_short_frame, short_tasks):
            pass
    print("All 50 short-distance diagnostic frame plots rendered successfully.")

    # Render Long-Distance Frames
    print(f"\nRendering 50 Long-Distance diagnostic plots -> {long_fits_dir}")
    with multiprocessing.Pool(processes=n_cpus) as pool:
        for done_idx in pool.imap_unordered(plot_single_long_frame, long_tasks):
            pass
    print("All 50 long-distance diagnostic frame plots rendered successfully.")

    print("\n" + "=" * 80)
    print("  ALL 100 DIAGNOSTIC FRAME PLOTS GENERATED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == '__main__':
    generate_all_diagnostic_frame_plots()
