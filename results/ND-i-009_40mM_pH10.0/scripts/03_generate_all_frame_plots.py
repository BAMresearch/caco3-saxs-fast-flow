#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 03: Universal Diagnostic Frame Plot Generator (ND-i-009).
----------------------------------------------------------------
Generates universal-scale, publication-grade diagnostic plots for all frames:
1. Short Distance:
   - Viscoelastic VPS fluctuations (Option 2: evolutionary xi_fluc & flowing Df, xi -> inf) [351 plots]
2. Long Distance:
   - Option C Forward Model (Constant Polydispersity sigma/R = 0.065) [350 plots]

Outputs saved to:
  output/short_distance/viscoelastic/fits/
  output/long_distance/option_c_forward/fits/

Uses parallel multiprocessing across all CPU cores.
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

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import decompose_model_components, RHO2_CALIBRATED


def plot_single_short_viscoelastic(args):
    f, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})
    valid = (y_frame > 0) & (dy_frame > 0)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = row_dict.get('time_seconds', 0.045 + 0.015 * f)
    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')

    redchi = row_dict.get('redchi', 0.0)
    Rmean = row_dict.get('Rmean', 8.5)
    D_val = row_dict.get('D', 1.90)
    xi_fluc = row_dict.get('xi_fluc', 2.5)
    A_fluc = row_dict.get('A_fluc', 0.1)
    xi_net = row_dict.get('xi', 1e7)

    param_dict = dict(row_dict)
    param_dict['rho2'] = RHO2_CALIBRATED
    param_dict['bkg'] = 0.010
    param_dict['total_volume_fraction_large'] = 0.0
    param_dict['mu_fluc'] = 4.0

    q_dense = np.logspace(np.log10(0.02), np.log10(2.6), 250)
    I_tot, I_prim, _, I_high, _ = decompose_model_components(param_dict, q_dense)

    xi_label = f"$\\xi={xi_net:.1f}$ nm" if xi_net < 1e5 else r"$\xi \to \infty$"
    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Total Fit ($\\chi^2={redchi:.2f}$)')
    ax_main.plot(q_dense, I_prim, color='#16a34a', linestyle='--', linewidth=1.2, label=f'Fractal Mesh ($R={Rmean:.1f}$ nm, $D={D_val:.2f}$, {xi_label})')
    if A_fluc > 0:
        ax_main.plot(q_dense, I_high, color='#2563eb', linestyle=':', linewidth=1.3, label=f'VPS Fluc. ($\\xi={xi_fluc:.1f}$ nm, $A={A_fluc:.2f}$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.02, 2.65)
    ax_main.set_ylim(5e-5, 1e4)
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11)
    ax_main.set_title(f"ND-i-009 Short Distance (3 m) Frame {f:03d} (t = {t_s:.3f} s): Viscoelastic VPS", fontsize=11, fontweight='bold')
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='lower left', fontsize=8.5, framealpha=0.92)

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



def plot_single_long_frame_forward(args):
    f, q_data, y_frame, dy_frame, fit_curve, row_dict, out_file = args
    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.08})
    valid = (y_frame > 0) & (q_data <= 0.22) & (dy_frame > 0)
    q_v = q_data[valid]
    y_v = y_frame[valid]
    dy_v = dy_frame[valid]
    fit_v = fit_curve[valid]

    t_s = row_dict.get('time_seconds', 0.045 + 0.015 * f)
    ax_main.errorbar(q_v, y_v, yerr=dy_v, fmt='o', color='#334155', markersize=2.2, alpha=0.35, elinewidth=0.5, capsize=0, label=f'Data (t = {t_s:.3f} s)')

    redchi = row_dict.get('redchi', 0.0)
    D_val = row_dict.get('D', 1.95)
    has_large = row_dict.get('has_large_spheres', True)
    phi_large = row_dict.get('total_volume_fraction_large', 0.0) if has_large else 0.0
    R_large = row_dict.get('mean_radius_large', 2800.0) if has_large else 0.0

    param_dict = dict(row_dict)
    param_dict['rho2'] = RHO2_CALIBRATED
    param_dict['bkg'] = 0.0
    param_dict['Rmean'] = 8.44
    param_dict['stdDev'] = 3.33
    param_dict['xi'] = 1e7
    param_dict['A_OZ'] = 0.0
    if not has_large:
        param_dict['total_volume_fraction_large'] = 0.0

    q_dense = np.logspace(np.log10(0.0018), np.log10(0.23), 250)
    I_tot, I_prim, I_large, _, _ = decompose_model_components(param_dict, q_dense)

    ax_main.plot(q_dense, I_tot, color='#dc2626', linewidth=2.0, label=f'Total Fit ($\\chi^2={redchi:.2f}$)')
    ax_main.plot(q_dense, I_prim, color='#16a34a', linestyle='--', linewidth=1.2, label=f'Mass Fractal ($D={D_val:.2f}, \\xi\\to\\infty$)')
    if has_large and phi_large > 0 and np.max(I_large) > 0:
        ax_main.plot(q_dense, I_large, color='#2563eb', linestyle='-.', linewidth=1.3, label=f'Droplets ($R={R_large:.0f}$ nm, $\\sigma/R=0.065$)')

    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlim(0.0018, 0.23)
    ax_main.set_ylim(5e-2, 2e6)
    ax_main.set_ylabel("Intensity $I(q)$ [cm$^{-1}$]", fontsize=11)
    ax_main.set_title(f"ND-i-009 Long Distance (31 m) Frame {f:03d} (t = {t_s:.3f} s): Option C Forward", fontsize=11, fontweight='bold')
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_main.legend(loc='lower left', fontsize=8.5, framealpha=0.92)

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
    print("  GENERATING ALL UNIVERSAL DIAGNOSTIC FIT PLOTS (SHORT & LONG DISTANCE)")
    print("=" * 80)
    start_t = time.time()

    # --- 1. Short Distance Plots (Viscoelastic) ---
    short_data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi009_fitting_results_short_distance.csv')
    df_short_data = pd.read_csv(short_data_path)
    q_short = df_short_data['q'].values

    dir_visco = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'viscoelastic')
    fits_dir_visco = os.path.join(dir_visco, 'fits')
    os.makedirs(fits_dir_visco, exist_ok=True)

    if os.path.exists(os.path.join(dir_visco, 'parameters.csv')):
        df_p_visco = pd.read_csv(os.path.join(dir_visco, 'parameters.csv'))
        df_c_visco = pd.read_csv(os.path.join(dir_visco, 'fitted_curves.csv'))
        p_visco_dict = {int(r['frame']): dict(r) for _, r in df_p_visco.iterrows()}
        frames_sd = sorted(list(p_visco_dict.keys()))

        tasks_visco = []
        for f in frames_sd:
            out_f = os.path.join(fits_dir_visco, f'fit_frame_{f:03d}.png')
            y_f = df_short_data[str(f)].values
            dy_f = df_short_data['error_' + str(f)].values
            fit_c = df_c_visco[str(f)].values
            tasks_visco.append((f, q_short, y_f, dy_f, fit_c, p_visco_dict[f], out_f))

        print(f"Generating {len(tasks_visco)} plots for Short-Distance Viscoelastic VPS...")
        with ProcessPoolExecutor(max_workers=14) as executor:
            list(executor.map(plot_single_short_viscoelastic, tasks_visco))

    # --- 2. Long Distance Plots (Option C Forward) ---
    long_data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi009_fitting_results_long_distance.csv')
    df_long_data = pd.read_csv(long_data_path)
    q_long = df_long_data['q'].values

    dir_ld = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'option_c_forward')
    fits_dir_ld = os.path.join(dir_ld, 'fits')
    os.makedirs(fits_dir_ld, exist_ok=True)

    if os.path.exists(os.path.join(dir_ld, 'parameters.csv')):
        df_p_ld = pd.read_csv(os.path.join(dir_ld, 'parameters.csv'))
        df_c_ld = pd.read_csv(os.path.join(dir_ld, 'fitted_curves.csv'))
        p_ld_dict = {int(r['frame']): dict(r) for _, r in df_p_ld.iterrows()}
        frames_ld = sorted(list(p_ld_dict.keys()))

        tasks_ld = []
        for f in frames_ld:
            out_f = os.path.join(fits_dir_ld, f'fit_frame_{f:03d}.png')
            # Only regenerate if missing
            if not os.path.exists(out_f):
                y_f = df_long_data[str(f)].values
                dy_f = df_long_data['error_' + str(f)].values
                fit_c = df_c_ld[str(f)].values
                tasks_ld.append((f, q_long, y_f, dy_f, fit_c, p_ld_dict[f], out_f))

        if tasks_ld:
            print(f"Generating {len(tasks_ld)} plots for Long-Distance Option C Forward...")
            with ProcessPoolExecutor(max_workers=14) as executor:
                list(executor.map(plot_single_long_frame_forward, tasks_ld))
        else:
            print("Long-distance frame plots already up-to-date.")

    elapsed = time.time() - start_t
    print("=" * 80)
    print(f"  ALL DIAGNOSTIC FRAME PLOTS GENERATED IN {elapsed:.1f} SECONDS!")
    print("=" * 80)


if __name__ == '__main__':
    main()

