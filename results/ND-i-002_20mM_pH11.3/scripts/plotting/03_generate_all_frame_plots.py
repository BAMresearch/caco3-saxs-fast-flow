#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 03: Universal Diagnostic Frame Plot Generator for ND-i-002.
------------------------------------------------------------------
Generates publication-quality diagnostic fit plots for all 50 frames:
1. Short Distance (3.0 m, frames 0 to 49):
   - Primary fractal network + Viscoelastic fluctuations + Residuals
2. Long Distance (31.0 m, frames 0 to 49):
   - Primary fractal network + Secondary droplet population + Residuals

Outputs:
  output/short_distance/fits/frame_XXX.png
  output/long_distance/fits/frame_XXX.png
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

from models.scattering_models import (
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_fb_intensity,
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
    Rmean = float(row_dict.get('Rmean', 10.0))
    stdDev = float(row_dict.get('stdDev', 4.0))
    D_val = float(row_dict.get('D', 1.93))
    xi_net = float(row_dict.get('xi', 400.0))
    phi_val = float(row_dict.get('phi', 2e-5))
    A_fluc = float(row_dict.get('A_fluc', 0.05))
    xi_fluc = float(row_dict.get('xi_fluc', 2.5))
    mu_fluc = float(row_dict.get('mu_fluc', 4.0))
    bkg = float(row_dict.get('bkg', 1e-3))

    q_dense = np.logspace(np.log10(0.02), np.log10(2.6), 250)
    I_prim = calc_primary_fractal_intensity(q_dense, Rmean, stdDev, phi_val, D_val, xi_net, rho2=RHO2_CALIBRATED, integration_steps=100)
    I_fluc = calc_fb_intensity(q_dense, A_fluc, xi_fluc, mu_fluc=mu_fluc)
    I_bkg = np.full_like(q_dense, bkg)
    I_tot = I_prim + I_fluc + I_bkg

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data ($t = {t_s:.3f}$ s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full model ($\\chi^2 = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_prim, color='#16a34a', linestyle='--', linewidth=1.3,
                 label=f'Primary fractal ($\\langle R \\rangle = {Rmean:.1f}$ nm, $D_f = {D_val:.2f}$)')
    ax_main.plot(q_dense, I_fluc, color='#2563eb', linestyle=':', linewidth=1.4,
                 label=f'Viscoelastic fluctuation ($A = {A_fluc:.3f}$, $\\xi_{{\\mathrm{{fluc}}}} = {xi_fluc:.1f}$ nm)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.02, 2.65)
    ax_main.set_ylim(2e-4, 5e2)
    ax_main.set_ylabel(r"Intensity $I(q)$ [$\mathrm{cm}^{-1}$]", fontsize=11)
    ax_main.set_title(f"ND-i-002 short distance frame {f_idx:02d} ($t = {t_s:.3f}$ s): multiscale model", fontsize=11)
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='lower left', fontsize=8.5, framealpha=0.92)

    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.02, 2.65)
    ax_res.set_ylim(-4.5, 4.5)
    ax_res.set_xlabel(r"Scattering vector $q$ [$\mathrm{nm}^{-1}$]", fontsize=11)
    ax_res.set_ylabel(r"Res. [$\sigma$]", fontsize=9.5)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f_idx


def plot_single_long_frame(args):
    f_idx, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})
    
    valid = (y_frame > 0) & (dy_frame > 0) & (q_data >= 0.002) & (q_data <= 0.18)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = float(row_dict.get('time_seconds', 0.0))
    redchi = float(row_dict.get('chi2_red', 0.0))
    Rmean = float(row_dict.get('Rmean_transferred', 10.0))
    stdDev = float(row_dict.get('stdDev_transferred', 4.0))
    D_val = float(row_dict.get('D', 1.93))
    xi_net = float(row_dict.get('xi', 400.0))
    phi_val = float(row_dict.get('phi', 2e-5))
    R_large = float(row_dict.get('mean_radius_large', 4900.0))
    std_large = float(row_dict.get('std_radius_large', 370.0))
    phi_large = float(row_dict.get('total_volume_fraction_large', 1e-3))
    A_fluc = float(row_dict.get('A_fluc_transferred', 0.05))
    xi_fluc = float(row_dict.get('xi_fluc_transferred', 2.5))

    A_porod = float(row_dict.get('A_porod', 0.0))
    R_large_2 = float(row_dict.get('mean_radius_large_2', 0.0))
    std_large_2 = float(row_dict.get('std_radius_large_2', 0.0))
    phi_large_2 = float(row_dict.get('total_volume_fraction_large_2', 0.0))

    q_dense = np.logspace(np.log10(0.0018), np.log10(0.25), 250)
    I_prim = calc_primary_fractal_intensity(q_dense, Rmean, stdDev, phi_val, D_val, xi_net, rho2=RHO2_CALIBRATED, integration_steps=100)
    I_large = calc_large_sphere_intensity(q_dense, R_large, std_large, phi_large, rho2=RHO2_CALIBRATED, integration_steps=100)
    if phi_large_2 > 0 and R_large_2 > 0:
        I_large_2 = calc_large_sphere_intensity(q_dense, R_large_2, std_large_2, phi_large_2, rho2=RHO2_CALIBRATED, integration_steps=100)
        I_large = I_large + I_large_2
        drop_label = f'Bimodal droplets ($\\langle R_1 \\rangle = {R_large:.0f}$ nm, $\\langle R_2 \\rangle = {R_large_2:.0f}$ nm)'
    else:
        drop_label = f'Droplets/Porod ($\\langle R \\rangle = {R_large:.0f}$ nm, $\\sigma = {std_large:.0f}$ nm)'

    I_porod = calc_porod_intensity(q_dense, A_porod, 4.0) if A_porod > 0 else np.zeros_like(q_dense)
    I_fluc = calc_fb_intensity(q_dense, A_fluc, xi_fluc, mu_fluc=4.0)
    I_tot = I_prim + I_large + I_porod + I_fluc

    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data ($t = {t_s:.3f}$ s)')
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Full model ($\\chi^2 = {redchi:.2f}$)')
    ax_main.plot(q_dense, I_large + I_porod, color='#7c3aed', linestyle='-.', linewidth=1.4,
                 label=drop_label)
    ax_main.plot(q_dense, I_prim, color='#16a34a', linestyle='--', linewidth=1.3,
                 label=f'Network fractal ($\\xi = {xi_net:.0f}$ nm, $D_f = {D_val:.2f}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.0018, 0.25)
    ax_main.set_ylim(max(1e-3, y_v.min() * 0.3), y_v.max() * 2.5)
    ax_main.set_ylabel(r"Intensity $I(q)$ [$\mathrm{cm}^{-1}$]", fontsize=11)
    ax_main.set_title(f"ND-i-002 long distance frame {f_idx:02d} ($t = {t_s:.3f}$ s): multiscale model", fontsize=11)
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='lower left', fontsize=8.5, framealpha=0.92)

    residuals = (y_v - fit_v) / dy_v
    ax_res.axhline(0, color='#64748b', linestyle='--', linewidth=0.9)
    ax_res.axhspan(-2, 2, color='#cbd5e1', alpha=0.35)
    ax_res.plot(q_v, residuals, '.', color='#dc2626', markersize=2.5, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 0.25)
    ax_res.set_ylim(-4.5, 4.5)
    ax_res.set_xlabel(r"Scattering vector $q$ [$\mathrm{nm}^{-1}$]", fontsize=11)
    ax_res.set_ylabel(r"Res. [$\sigma$]", fontsize=9.5)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_file, dpi=120)
    plt.close(fig)
    return f_idx


def generate_all_diagnostic_frame_plots():
    print("=" * 80)
    print("  GENERATING UNIVERSAL DIAGNOSTIC FRAME PLOTS FOR ND-i-002 (50 SD + 50 LD)")
    print("=" * 80)
    t0 = time.time()

    sd_fits_dir = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fits')
    ld_fits_dir = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fits')
    os.makedirs(sd_fits_dir, exist_ok=True)
    os.makedirs(ld_fits_dir, exist_ok=True)

    # 1. Short Distance
    sd_data = pd.read_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi002_fitting_results_short_distance.csv'))
    sd_curves = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'short_distance', 'fitted_curves.csv'))
    sd_params = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv'))
    q_sd = sd_data['q'].values

    sd_tasks = []
    for _, row in sd_params.iterrows():
        f = int(row['frame_index'])
        out_f = os.path.join(sd_fits_dir, f"frame_{f:03d}.png")
        sd_tasks.append((f, q_sd, sd_data[str(f)].values, sd_data['error_' + str(f)].values,
                         sd_curves[str(f)].values, row.to_dict(), out_f))

    # 2. Long Distance
    ld_data = pd.read_csv(os.path.join(PROJECT_DIR, 'data', 'dataNDi002_fitting_results_long_distance.csv'))
    ld_curves = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'long_distance', 'fitted_curves.csv'))
    ld_params = pd.read_csv(os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv'))
    q_ld = ld_data['q'].values

    ld_tasks = []
    for _, row in ld_params.iterrows():
        f = int(row['frame_index'])
        out_f = os.path.join(ld_fits_dir, f"frame_{f:03d}.png")
        ld_tasks.append((f, q_ld, ld_data[str(f)].values, ld_data['error_' + str(f)].values,
                         ld_curves[str(f)].values, row.to_dict(), out_f))

    max_workers = max(1, (os.cpu_count() or 4) - 1)
    print(f"Executing parallel diagnostic rendering on {max_workers} CPU workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        print("Rendering 50 Short-Distance diagnostic frame plots...")
        for _ in pool.map(plot_single_short_frame, sd_tasks):
            pass

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        print("Rendering 50 Long-Distance diagnostic frame plots...")
        for _ in pool.map(plot_single_long_frame, ld_tasks):
            pass

    elapsed = time.time() - t0
    print(f"Generated all 100 diagnostic plots in {elapsed:.1f} s.")


if __name__ == '__main__':
    generate_all_diagnostic_frame_plots()
