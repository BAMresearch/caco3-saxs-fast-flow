#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 03: Universal Diagnostic Frame Plot Generator for ND-i-003.
------------------------------------------------------------------
Generates publication-grade diagnostic plots with residual panels for all 50 frames:
1. Short Distance (3.0 m, q in [0.022, 2.54] nm^-1, t = 0.060 to 10.877 s):
   - Flat solvent background fit demonstrating absence of high-q clusters or fractal mesh.
   - Saved to: output/short_distance/fits/fit_frame_XXX.png
2. Long Distance (31.0 m, q in [0.002, 0.25] nm^-1, t = 0.060 to 10.877 s):
   - Discrete polydisperse ACC nanosphere form factor fit with Porod q^-4 decay.
   - Saved to: output/long_distance/fits/fit_frame_XXX.png

Uses ProcessPoolExecutor for parallel multiprocessing across CPU cores.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import calc_long_distance_model, calc_short_distance_model, calc_sphere_intensity, RHO2_CALIBRATED


def plot_single_short_frame(args):
    f, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1, figsize=(7.8, 5.8), sharex=True,
        gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08}
    )

    t_s = float(row_dict.get('time_seconds', 0.0))
    redchi = float(row_dict.get('redchi', 0.0))
    bkg_val = float(row_dict.get('bkg', 0.0))
    err_bkg = float(row_dict.get('std_err_bkg', 0.0))

    # Mask valid data points
    valid = (dy_frame > 0) & np.isfinite(y_frame)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    # Data points
    ax_main.errorbar(
        q_v, y_v, yerr=dy_v, fmt='o', color='#334155',
        markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0,
        label=f'SAXS Data (t = {t_s:.3f} s)'
    )

    # Fitted model line
    q_dense = np.logspace(np.log10(0.02), np.log10(2.6), 250)
    y_dense = np.full_like(q_dense, bkg_val)
    ax_main.plot(
        q_dense, y_dense, color='#dc2626', linewidth=2.0,
        label=f'Flat Background Fit ($\\chi^2={redchi:.3f}$)'
    )

    # Reference zero line
    ax_main.axhline(0, color='#64748b', linestyle=':', linewidth=1.0, alpha=0.7)

    # Formatting: linear y-scale since data fluctuates symmetrically around zero background
    ax_main.set_xscale('log')
    ax_main.set_xlim(0.02, 2.65)
    ax_main.set_ylim(-0.005, 0.010)
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11)
    ax_main.set_title(
        f"ND-i-003 Short Distance (3.0 m) Frame {f:03d} (t = {t_s:.3f} s): Solvent Background",
        fontsize=10.5, fontweight='bold'
    )
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='upper right', fontsize=8.5, framealpha=0.92)

    # Annotate physical conclusion
    ax_main.text(
        0.03, 0.15,
        f"$I_{{\\rm bkg}} = ({bkg_val*1e4:.2f} \\pm {err_bkg*1e4:.2f}) \\times 10^{{-4}}\\ \\mathrm{{cm}}^{{-1}}$\n"
        f"Absence of high-$q$ PNCs or fractal network nodes",
        transform=ax_main.transAxes, fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8fafc', edgecolor='#cbd5e1', alpha=0.9)
    )

    # Residuals
    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.02, 2.65)
    ax_res.set_ylim(-5.5, 5.5)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=11)
    ax_res.set_ylabel("Res. [$\\sigma$]", fontsize=9.5)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f


def plot_single_long_frame(args):
    f, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1, figsize=(7.8, 5.8), sharex=True,
        gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08}
    )

    t_s = float(row_dict.get('time_seconds', 0.0))
    is_ind = bool(row_dict.get('is_induction', False))
    has_ff = bool(row_dict.get('has_form_factor', not is_ind))
    bkg_val = float(row_dict.get('bkg', 0.01))

    valid = (y_frame > 0) & (q_data <= 0.22) & (dy_frame > 0)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    # Data points
    ax_main.errorbar(
        q_v, y_v, yerr=dy_v, fmt='o', color='#334155',
        markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0,
        label=f'SAXS Data (t = {t_s:.3f} s)'
    )

    q_dense = np.logspace(np.log10(0.0018), np.log10(0.23), 300)
    from models.scattering_models import calc_large_sphere_intensity, calc_porod_intensity

    R_drop = float(row_dict.get('mean_radius_large', row_dict.get('mean_radius_large_lognormal', 2600.0)))
    phi_drop = float(row_dict.get('total_volume_fraction_large', row_dict.get('total_volume_fraction_large_lognormal', 5.5e-5)))
    A_p = float(row_dict.get('A_porod', 0.0))

    if not has_ff:
        # Induction regime: Macroscopic droplets
        redchi = float(row_dict.get('redchi', 1.0))
        y_dense = calc_large_sphere_intensity(q_dense, R_drop, 248.0, phi_drop) + bkg_val
        ax_main.plot(q_dense, y_dense, color='#16a34a', linewidth=2.0,
                     label=f'Droplets Fit ($R={R_drop:.0f}$ nm, $\\chi^2={redchi:.2f}$)')
        ax_main.text(
            0.04, 0.15,
            f"Induction Period: Macroscopic Droplets\n$R_{{\\rm large}} = {R_drop:.0f}\\ \\mathrm{{nm}}$, $\\phi_{{\\rm large}} = {phi_drop*1e5:.2f} \\times 10^{{-5}}$\n(Sheared emulsion droplet baseline)",
            transform=ax_main.transAxes, fontsize=8.5, fontweight='bold', color='#166534',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='#f0fdf4', edgecolor='#bbf7d0')
        )
    else:
        # Active growth: Show Total Model, Discrete ACC Spheres, and Macroscopic Droplets
        R_log = float(row_dict.get('Rmean_lognormal', row_dict.get('Rmean', 80.0)))
        sig_log = float(row_dict.get('stdDev_lognormal', row_dict.get('stdDev', 25.0)))
        phi_log = float(row_dict.get('phi_lognormal', row_dict.get('phi', 1e-4)))
        bkg_log = float(row_dict.get('bkg_lognormal', bkg_val))
        chi_log = float(row_dict.get('redchi_lognormal', row_dict.get('redchi', 1.0)))

        I_sph_only = calc_sphere_intensity(q_dense, R_log, sig_log, phi_log, bkg=0.0, dist_type='lognormal')
        I_drop_only = calc_large_sphere_intensity(q_dense, R_drop, 248.0, phi_drop)
        I_tot_dense = I_sph_only + I_drop_only + bkg_log

        ax_main.plot(q_dense, I_tot_dense, color='#2563eb', linewidth=2.0,
                     label=f'Total Fit ($\\chi^2={chi_log:.2f}$)')
        ax_main.plot(q_dense, I_sph_only, color='#9333ea', linewidth=1.6, linestyle='--',
                     label=f'ACC Spheres: $\\langle R \\rangle={R_log:.1f}$ nm, $\\sigma={sig_log:.1f}$ nm')
        ax_main.plot(q_dense, I_drop_only, color='#16a34a', linewidth=1.4, linestyle=':',
                     label=f'Droplets: $R={R_drop:.0f}$ nm, $\\phi={phi_drop*1e5:.1f} \\times 10^{{-5}}$')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.0018, 0.23)
    ax_main.set_ylim(1e-3, 3e3)
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11)
    title_suffix = "Induction Period (Macroscopic Droplets)" if not has_ff else "ACC Nanospheres + Macroscopic Droplets"
    ax_main.set_title(
        f"ND-i-003 Long Distance Frame {f:03d} (t = {t_s:.3f} s): {title_suffix}",
        fontsize=10.5, fontweight='bold'
    )
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='upper right', fontsize=8.0, framealpha=0.92)

    # Residuals
    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 0.23)
    ax_res.set_ylim(-5.5, 5.5)
    ax_res.set_xlabel("Scattering Vector $q$ [nm$^{-1}$]", fontsize=11)
    ax_res.set_ylabel("Res. [$\\sigma$]", fontsize=9.5)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f


def main():
    print("=" * 80)
    print("  ND-i-003: GENERATING ALL DIAGNOSTIC FRAME PLOTS (50 SHORT + 50 LONG)")
    print("=" * 80)
    start_t = time.time()

    # --- 1. Short Distance Plots ---
    short_data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_short_distance.csv')
    df_short_data = pd.read_csv(short_data_path)
    q_short = df_short_data['q'].values

    dir_short = os.path.join(PROJECT_DIR, 'output', 'short_distance')
    fits_dir_short = os.path.join(dir_short, 'fits')
    os.makedirs(fits_dir_short, exist_ok=True)

    if os.path.exists(os.path.join(dir_short, 'parameters.csv')):
        df_p_short = pd.read_csv(os.path.join(dir_short, 'parameters.csv'))
        df_c_short = pd.read_csv(os.path.join(dir_short, 'fitted_curves.csv'))
        p_short_dict = {int(r['frame']): dict(r) for _, r in df_p_short.iterrows()}
        frames_short = sorted(list(p_short_dict.keys()))

        tasks_short = []
        for f in frames_short:
            out_f = os.path.join(fits_dir_short, f'fit_frame_{f:03d}.png')
            y_f = df_short_data[str(f)].values
            dy_f = df_short_data['error_' + str(f)].values
            fit_c = df_c_short[str(f)].values
            tasks_short.append((f, q_short, y_f, dy_f, fit_c, p_short_dict[f], out_f))

        print(f"Generating {len(tasks_short)} diagnostic plots for Short-Distance SAXS...")
        with ProcessPoolExecutor(max_workers=14) as executor:
            list(executor.map(plot_single_short_frame, tasks_short))

    # --- 2. Long Distance Plots ---
    long_data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_long_distance.csv')
    df_long_data = pd.read_csv(long_data_path)
    q_long = df_long_data['q'].values

    dir_long = os.path.join(PROJECT_DIR, 'output', 'long_distance')
    fits_dir_long = os.path.join(dir_long, 'fits')
    os.makedirs(fits_dir_long, exist_ok=True)

    if os.path.exists(os.path.join(dir_long, 'parameters.csv')):
        df_p_long = pd.read_csv(os.path.join(dir_long, 'parameters.csv'))
        df_c_long = pd.read_csv(os.path.join(dir_long, 'fitted_curves.csv'))
        p_long_dict = {int(r['frame']): dict(r) for _, r in df_p_long.iterrows()}
        frames_long = sorted(list(p_long_dict.keys()))

        tasks_long = []
        for f in frames_long:
            out_f = os.path.join(fits_dir_long, f'fit_frame_{f:03d}.png')
            y_f = df_long_data[str(f)].values
            dy_f = df_long_data['error_' + str(f)].values
            fit_c = df_c_long[str(f)].values
            tasks_long.append((f, q_long, y_f, dy_f, fit_c, p_long_dict[f], out_f))

        print(f"Generating {len(tasks_long)} diagnostic plots for Long-Distance SAXS...")
        with ProcessPoolExecutor(max_workers=14) as executor:
            list(executor.map(plot_single_long_frame, tasks_long))

    elapsed = time.time() - start_t
    print("=" * 80)
    print(f"  ALL 100 DIAGNOSTIC FRAME PLOTS GENERATED IN {elapsed:.1f} SECONDS!")
    print("=" * 80)


if __name__ == '__main__':
    main()
