#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Long Distance Time-Resolved Fitting for ND-i-006.
-----------------------------------------------------------
Fits the 50 long-distance SAXS frames (31.0 m camera, q in [0.002, 0.22] nm^-1, t = 0.060 s to 10.877 s).

Scientific Principle for 40.0 mM CaCO3 (pH 9.75):
- Mesoscale spinodal / mass-fractal percolation (Dm ~ 1.8 - 2.4, xi ~ 200 - 30000 nm).
- Secondary macroscopic droplets / large spheres: mean_radius_large (~3000 nm), std_radius_large, total_volume_fraction_large.
- Macroscopic domain coarsening: Asymptotic Porod scattering (I ~ A_porod / q^4) for giant domains where q*R >> 1.
- Incoherent background (bkg).
- NO ACC NANOSPHERES (experimentally absent in this condition).
- Sequential continuity regularization prevents unphysical parameter jumps.
- Parameter standard errors (+/- 1 sigma) computed via the covariance matrix.

Outputs:
  output/long_distance/parameters.csv
  output/long_distance/fitted_curves.csv
  output/long_distance/trajectories.png
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_long_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'long_distance')


def fit_long_distance_series():
    print("=" * 80)
    print("  LONG-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-006, 50 FRAMES)")
    print("  Model: Spinodal Mass Fractal Network + Large Droplets / Spheres (No ACC Particles)")
    print("=" * 80)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    df_short = pd.read_csv(SHORT_CSV)

    d_base = 1.80
    xi_base = 30000.0
    r_prim_base = 1.50
    std_prim_base = 3.03
    r_large_base = 3113.0
    std_large_base = 267.0
    phi_large_base = 1.12e-4
    bkg_base = 2.3e-4

    if os.path.exists(BASELINE_JSON):
        try:
            with open(BASELINE_JSON, 'r') as f:
                bp = json.load(f)
            d_base = bp.get('Dm', bp.get('D', 1.80))
            xi_base = bp.get('xi', 30000.0)
            r_prim_base = bp.get('R_prim', 1.50)
            std_prim_base = bp.get('std_prim', 3.03)
            r_large_base = bp.get('mean_radius_large', 3113.0)
            std_large_base = bp.get('std_radius_large', 267.0)
            phi_large_base = bp.get('total_volume_fraction_large', 1.12e-4)
            bkg_base = bp.get('bkg', 2.3e-4)
        except Exception:
            pass

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Baseline anchor: Dm = {d_base:.3f}, R_large = {r_large_base:.1f} nm, phi_large = {phi_large_base:.2e}")

    fitted_records = []
    fitted_curves = {'q': q}
    prev_params = None

    for pos in reversed(range(n_frames)):
        row = df_meta.iloc[pos]
        f_idx = int(row['frame_index'])
        t_sec = float(row['time_seconds'])
        f_file = str(row['long_file'])

        y = df_data[str(f_idx)].values
        dy = np.maximum(df_data['error_' + str(f_idx)].values, 1e-4)

        mask = (q >= 0.0025) & (q <= 0.22) & (y > 0)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = dy[mask]

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)

        # Transferred building block parameters
        p.add('R_prim', value=r_prim_base, vary=False)
        p.add('std_prim', value=std_prim_base, vary=False)
        p.add('std_radius_large', value=std_large_base, vary=False)

        if prev_params is None:
            p.add('phi_net', value=2.5e-5, min=1e-7, max=0.01, vary=True)
            p.add('Dm', value=d_base, min=1.1, max=2.85, vary=True)
            p.add('xi', value=min(xi_base, 5000.0), min=20.0, max=50000.0, vary=True)
            p.add('mean_radius_large', value=r_large_base, min=1000.0, max=8000.0, vary=True)
            p.add('total_volume_fraction_large', value=phi_large_base, min=1e-8, max=0.005, vary=True)
            p.add('A_porod', value=1.0e-6, min=1e-10, max=0.01, vary=True)
            p.add('bkg', value=bkg_base, min=1e-6, max=0.01, vary=True)
        else:
            dm_val = prev_params['Dm'].value
            xi_val = prev_params['xi'].value
            phi_net_val = prev_params['phi_net'].value
            r_large_val = prev_params['mean_radius_large'].value
            phi_large_val = prev_params['total_volume_fraction_large'].value
            a_porod_val = prev_params['A_porod'].value
            bkg_val = prev_params['bkg'].value

            p.add('phi_net', value=phi_net_val, min=1e-7, max=0.01, vary=True)
            p.add('Dm', value=dm_val, min=max(1.1, dm_val - 0.25), max=min(2.85, dm_val + 0.25), vary=True)
            p.add('xi', value=xi_val, min=max(20.0, xi_val * 0.4), max=min(50000.0, xi_val * 2.5), vary=True)
            p.add('mean_radius_large', value=r_large_val, min=max(1000.0, r_large_val - 800.0), max=min(8000.0, r_large_val + 800.0), vary=True)
            p.add('total_volume_fraction_large', value=phi_large_val, min=1e-8, max=0.005, vary=True)
            p.add('A_porod', value=a_porod_val, min=1e-10, max=0.01, vary=True)
            p.add('bkg', value=bkg_val, min=1e-6, max=0.01, vary=True)

        def model_ld(pars, q_eval):
            I_frac = calc_primary_fractal_intensity(
                q_eval, pars['R_prim'].value, pars['std_prim'].value,
                pars['phi_net'].value, pars['Dm'].value, pars['xi'].value,
                integration_steps=60
            )
            I_large = calc_large_sphere_intensity(
                q_eval, pars['mean_radius_large'].value, pars['std_radius_large'].value,
                pars['total_volume_fraction_large'].value,
                integration_steps=60
            )
            I_porod = calc_porod_intensity(q_eval, pars['A_porod'].value, 4.0)
            return I_frac + I_large + I_porod + pars['bkg'].value

        def res_long(pars):
            c = model_ld(pars, q_fit)
            residuals = (y_fit - c) / dy_fit
            if prev_params is not None:
                dDm = (pars['Dm'].value - prev_params['Dm'].value) / 0.15
                dxi = (np.log(pars['xi'].value) - np.log(prev_params['xi'].value)) / 0.35
                dRlarge = (pars['mean_radius_large'].value - prev_params['mean_radius_large'].value) / 300.0
                return np.concatenate([residuals, [dDm, dxi, dRlarge]])
            return residuals

        out_fit = minimize(res_long, p, method='least_squares', max_nfev=200)

        p_cov = out_fit.params.copy()
        for k in ['phi_net', 'Dm', 'xi', 'mean_radius_large', 'total_volume_fraction_large', 'A_porod', 'bkg']:
            p_cov[k].vary = True
        try:
            out_cov = minimize(res_long, p_cov, method='leastsq', calc_covar=True)
            final_p = out_cov.params
        except Exception:
            final_p = out_fit.params

        fit_full = model_ld(final_p, q)
        fit_fit = model_ld(final_p, q_fit)
        chi2_red = float(np.sum(((y_fit - fit_fit) / dy_fit) ** 2) / max(1, len(y_fit) - 6))

        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'filename': f_file,
            'chi2_red': chi2_red,
            'phi_net': final_p['phi_net'].value,
            'std_err_phi_net': float(final_p['phi_net'].stderr) if final_p['phi_net'].stderr is not None else 0.0,
            'Dm': final_p['Dm'].value,
            'std_err_Dm': float(final_p['Dm'].stderr) if final_p['Dm'].stderr is not None else 0.0,
            'xi': final_p['xi'].value,
            'std_err_xi': float(final_p['xi'].stderr) if final_p['xi'].stderr is not None else 0.0,
            'mean_radius_large': final_p['mean_radius_large'].value,
            'std_err_mean_radius_large': float(final_p['mean_radius_large'].stderr) if final_p['mean_radius_large'].stderr is not None else 0.0,
            'std_radius_large': final_p['std_radius_large'].value,
            'total_volume_fraction_large': final_p['total_volume_fraction_large'].value,
            'std_err_total_volume_fraction_large': float(final_p['total_volume_fraction_large'].stderr) if final_p['total_volume_fraction_large'].stderr is not None else 0.0,
            'A_porod': final_p['A_porod'].value,
            'std_err_A_porod': float(final_p['A_porod'].stderr) if final_p['A_porod'].stderr is not None else 0.0,
            'R_prim': final_p['R_prim'].value,
            'std_prim': final_p['std_prim'].value,
            'bkg': final_p['bkg'].value,
            'std_err_bkg': float(final_p['bkg'].stderr) if final_p['bkg'].stderr is not None else 0.0,
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_full
        prev_params = final_p

        if f_idx % 10 == 0 or f_idx == n_frames - 1 or f_idx == 0:
            print(f"  Frame {f_idx:02d} (t = {t_sec:6.3f} s): chi2 = {chi2_red:.3f} | Dm = {final_p['Dm'].value:.3f} ± {rec['std_err_Dm']:.3f}, xi = {final_p['xi'].value:.0f} nm, R_large = {final_p['mean_radius_large'].value:.0f} nm")

    fitted_records.sort(key=lambda x: x['frame_index'])
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)

    df_curves = pd.DataFrame(fitted_curves)
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)

    # Overview Trajectory Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("ND-i-006 Long-Distance Kinetic Evolution: Fractal Network & Large Droplets (40.0 mM CaCO3, pH 9.75)", fontsize=13, fontweight='bold')
    t_arr = df_params['time_seconds']

    # Panel 1: Dm
    axes[0, 0].plot(t_arr, df_params['Dm'], 'o-', color='#d97706', label=r'$D_m$')
    axes[0, 0].fill_between(t_arr, df_params['Dm'] - df_params['std_err_Dm'], df_params['Dm'] + df_params['std_err_Dm'], color='#fbbf24', alpha=0.3)
    axes[0, 0].set_ylabel(r"Mass Fractal Dimension $D_m$", fontweight='bold')
    axes[0, 0].set_title(r"Fractal Dimension $D_m(t)$", fontweight='bold')
    axes[0, 0].grid(True, linestyle='--', alpha=0.5)
    axes[0, 0].legend()

    # Panel 2: Correlation Length xi & Large Sphere Radius
    axes[0, 1].plot(t_arr, df_params['xi'], 's-', color='#2563eb', label=r'$\xi$ [nm]')
    axes[0, 1].plot(t_arr, df_params['mean_radius_large'], 'd-', color='#0284c7', label=r'$\langle R_{\mathrm{large}} \rangle$ [nm]')
    axes[0, 1].set_ylabel("Length Scale [nm]", fontweight='bold')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_title(r"Network Cutoff $\xi(t)$ & Large Sphere Radius", fontweight='bold')
    axes[0, 1].grid(True, which='both', linestyle='--', alpha=0.5)
    axes[0, 1].legend()

    # Panel 3: Volume Fractions
    axes[1, 0].plot(t_arr, df_params['phi_net'] * 1e4, '^-', color='#16a34a', label=r'$\phi_{\mathrm{net}} \times 10^4$')
    axes[1, 0].plot(t_arr, df_params['total_volume_fraction_large'] * 1e4, 'o-', color='#0284c7', label=r'$\phi_{\mathrm{large}} \times 10^4$')
    axes[1, 0].set_ylabel(r"$\phi \times 10^4$", fontweight='bold')
    axes[1, 0].set_xlabel("Time [s]", fontweight='bold')
    axes[1, 0].set_title("Volume Fraction Kinetics", fontweight='bold')
    axes[1, 0].grid(True, linestyle='--', alpha=0.5)
    axes[1, 0].legend()

    # Panel 4: Chi2
    axes[1, 1].plot(t_arr, df_params['chi2_red'], 'v-', color='#dc2626', label=r'$\chi^2_{\mathrm{red}}$')
    axes[1, 1].axhline(1.0, color='gray', linestyle=':')
    axes[1, 1].set_ylabel(r"Reduced $\chi^2$", fontweight='bold')
    axes[1, 1].set_xlabel("Time [s]", fontweight='bold')
    axes[1, 1].set_title("Goodness-of-Fit Metric", fontweight='bold')
    axes[1, 1].grid(True, linestyle='--', alpha=0.5)
    axes[1, 1].legend()

    plt.tight_layout()
    traj_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(traj_path, dpi=300)
    plt.close()

    elapsed = time.time() - start_time
    print(f"Long distance fitting completed in {elapsed:.1f} s. Saved to {OUT_DIR}")


if __name__ == '__main__':
    fit_long_distance_series()
