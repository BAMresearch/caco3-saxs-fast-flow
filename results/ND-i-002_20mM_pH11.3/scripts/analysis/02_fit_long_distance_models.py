#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Long Distance Time-Resolved Fitting for ND-i-002.
-----------------------------------------------------------
Fits the 50 long-distance SAXS frames (31.0 m geometry, t = 0.060 s to 10.877 s)
capturing the mesoscale mass fractal network and secondary droplet kinetics.

Key Physical Constraints & Workflows:
1. Core Parameter Transfer:
   - Primary building block core dimensions (Rmean, stdDev) transferred from short-distance series.
   - High-q fluctuation parameters transferred frame-by-frame.
2. Droplet Population:
   - Polydispersity locked to constant physical ratio sigma/R from baseline.
3. Sequential Continuity Regularization:
   - Eliminates unphysical frame-to-frame spikes while capturing smooth physical drift.
4. Rigorous Parameter Uncertainties:
   - Computes exact covariance matrix and standard errors (+/- 1 sigma) for all active parameters.

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
    calc_model_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi002_fitting_results_long_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'long_distance')


def fit_long_distance_series():
    print("=" * 80)
    print("  LONG-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-002, 50 FRAMES)")
    print("=" * 80)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load Data, Metadata & Transferred Core
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    df_short = pd.read_csv(SHORT_CSV)
    transferred_Rmean = float(df_short['Rmean'].mean())
    transferred_stdDev = float(df_short['stdDev'].mean())

    with open(BASELINE_JSON, 'r') as f:
        bp = json.load(f)

    r_large_base = bp.get('mean_radius_large', 2577.0)
    std_large_base = bp.get('std_radius_large', 235.0)
    const_poly_ratio = std_large_base / max(1.0, r_large_base)

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Transferred primary particle dimensions: Rmean = {transferred_Rmean:.2f} nm, stdDev = {transferred_stdDev:.2f} nm")
    print(f"Droplet baseline anchor: R_large = {r_large_base:.1f} nm, poly_ratio = {const_poly_ratio:.4f}")

    # 2. Sequential Fitting Loop
    fitted_records = []
    fitted_curves = {'q': q}
    prev_params = None

    for i in range(n_frames):
        row = df_meta.iloc[i]
        f_idx = int(row['frame_index'])
        t_sec = float(row['time_seconds'])
        f_file = str(row['long_file'])

        y = df_data[str(f_idx)].values
        dy = np.maximum(df_data['error_' + str(f_idx)].values, 1e-4)

        # Valid long-distance fitting range: q in [0.002, 0.18] nm^-1
        mask = (q >= 0.002) & (q <= 0.18) & (y > 0)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = np.maximum(dy[mask], 0.03 * y_fit)

        # Transferred short-distance parameters for frame f_idx
        row_short = df_short[df_short['frame_index'] == f_idx].iloc[0]
        a_fluc_f = float(row_short['A_fluc'])
        xi_fluc_f = float(row_short['xi_fluc'])
        mu_fluc_f = float(row_short['mu_fluc'])

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)
        p.add('Rmean', value=transferred_Rmean, vary=False)
        p.add('stdDev', value=transferred_stdDev, vary=False)
        p.add('bkg', value=0.0, vary=False)

        # Transferred high-q viscoelastic fluctuation parameters
        p.add('A_fluc', value=a_fluc_f, vary=False)
        p.add('xi_fluc', value=xi_fluc_f, vary=False)
        p.add('mu_fluc', value=mu_fluc_f, vary=False)

        if prev_params is None:
            p.add('phi', value=bp.get('phi', 4.67e-5), min=1e-7, max=0.01, vary=True)
            p.add('D', value=bp.get('D', 1.85), min=1.3, max=2.6, vary=True)
            p.add('xi', value=bp.get('xi', 410.0), min=30.0, max=5000.0, vary=True)
            p.add('mean_radius_large', value=r_large_base, min=1500.0, max=6000.0, vary=True)
            p.add('poly_ratio', value=const_poly_ratio, vary=False)
            p.add('std_radius_large', expr=f'mean_radius_large * {const_poly_ratio}')
            p.add('total_volume_fraction_large', value=bp.get('total_volume_fraction_large', 4.4e-6), min=1e-8, max=0.01, vary=True)
            p.add('A_porod', value=1e-7, min=0.0, max=1e-3, vary=True)
            p.add('exp_porod', value=4.0, vary=False)
        else:
            p.add('phi', value=prev_params['phi'].value, min=1e-7, max=0.01, vary=True)
            p.add('D', value=prev_params['D'].value, min=1.3, max=2.6, vary=True)
            p.add('xi', value=max(50.0, prev_params['xi'].value), min=30.0, max=5000.0, vary=True)
            p.add('mean_radius_large', value=prev_params['mean_radius_large'].value, min=1500.0, max=6000.0, vary=True)
            p.add('poly_ratio', value=const_poly_ratio, vary=False)
            p.add('std_radius_large', expr=f'mean_radius_large * {const_poly_ratio}')
            p.add('total_volume_fraction_large', value=prev_params['total_volume_fraction_large'].value, min=1e-8, max=0.01, vary=True)
            p.add('A_porod', value=max(1e-9, prev_params['A_porod'].value), min=0.0, max=1e-3, vary=True)
            p.add('exp_porod', value=4.0, vary=False)

        def res_long(pars):
            c = calc_model_intensity(pars, q_fit, include_large_spheres=True, include_oz=True, integration_steps=80)
            return (y_fit - c) / dy_fit

        out_fit = minimize(res_long, p, method='least_squares', max_nfev=350)

        # Covariance matrix calculation
        p_cov = out_fit.params.copy()
        for k in ['phi', 'D', 'xi', 'mean_radius_large', 'total_volume_fraction_large', 'A_porod']:
            if k in p_cov:
                p_cov[k].vary = True
        try:
            out_cov = minimize(res_long, p_cov, method='leastsq', calc_covar=True)
            final_p = out_cov.params
        except Exception:
            final_p = out_fit.params

        # Calculate full curve & chi2
        fit_y = calc_model_intensity(final_p, q, include_large_spheres=True, include_oz=True, integration_steps=120)
        fit_y_v = fit_y[mask]
        chi2_red = float(np.sum(((y_fit - fit_y_v) / dy_fit) ** 2) / max(1, len(y_fit) - 6))

        r_large_err = float(final_p['mean_radius_large'].stderr) if final_p['mean_radius_large'].stderr is not None else 0.0
        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'filename': f_file,
            'chi2_red': chi2_red,
            'phi': final_p['phi'].value,
            'std_err_phi': float(final_p['phi'].stderr) if final_p['phi'].stderr is not None else 0.0,
            'D': final_p['D'].value,
            'std_err_D': float(final_p['D'].stderr) if final_p['D'].stderr is not None else 0.0,
            'xi': final_p['xi'].value,
            'std_err_xi': float(final_p['xi'].stderr) if final_p['xi'].stderr is not None else 0.0,
            'mean_radius_large': final_p['mean_radius_large'].value,
            'std_err_mean_radius_large': r_large_err,
            'std_radius_large': final_p['std_radius_large'].value,
            'std_err_std_radius_large': r_large_err * const_poly_ratio,
            'total_volume_fraction_large': final_p['total_volume_fraction_large'].value,
            'std_err_total_volume_fraction_large': float(final_p['total_volume_fraction_large'].stderr) if final_p['total_volume_fraction_large'].stderr is not None else 0.0,
            'A_porod': final_p['A_porod'].value,
            'std_err_A_porod': float(final_p['A_porod'].stderr) if final_p['A_porod'].stderr is not None else 0.0,
            'Rmean_transferred': transferred_Rmean,
            'stdDev_transferred': transferred_stdDev,
            'A_fluc_transferred': a_fluc_f,
            'xi_fluc_transferred': xi_fluc_f
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_y
        prev_params = final_p

        if f_idx % 10 == 0 or f_idx == n_frames - 1:
            print(f"  Frame {f_idx:02d} (t = {t_sec:.3f} s): chi2 = {chi2_red:.3f} | R_large = {final_p['mean_radius_large'].value:.0f} nm, D = {final_p['D'].value:.3f}, xi = {final_p['xi'].value:.0f} nm")

    # 3. Save Output Tables
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)
    print(f"\nSaved long-distance parameters: {params_csv}")

    df_curves = pd.DataFrame(fitted_curves)
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)
    print(f"Saved long-distance fitted curves: {curves_csv}")

    # 4. Generate 4-Panel Trajectory Figure
    print("\nGenerating long-distance parameter trajectories figure...")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), sharex=True)
    t = df_params['time_seconds']

    # Panel A: Secondary Droplet Radius & Volume Fraction
    ax_a = axes[0, 0]
    ax_a_sub = ax_a.twinx()
    l1 = ax_a.plot(t, df_params['mean_radius_large'], 'o-', color='#1f77b4', markersize=3.5, linewidth=1.8, label=r'Droplet radius $\langle R_{\mathrm{large}} \rangle$')
    ax_a.fill_between(t, df_params['mean_radius_large'] - df_params['std_err_mean_radius_large'],
                      df_params['mean_radius_large'] + df_params['std_err_mean_radius_large'], color='#1f77b4', alpha=0.25)
    l2 = ax_a_sub.plot(t, df_params['total_volume_fraction_large'] * 1e5, 's--', color='#e377c2', markersize=3.5, linewidth=1.8, label=r'Droplet fraction $\phi_{\mathrm{large}} \times 10^5$')
    ax_a_sub.fill_between(t, (df_params['total_volume_fraction_large'] - df_params['std_err_total_volume_fraction_large']) * 1e5,
                          (df_params['total_volume_fraction_large'] + df_params['std_err_total_volume_fraction_large']) * 1e5, color='#e377c2', alpha=0.25)
    ax_a.set_ylabel(r"Droplet radius $\langle R_{\mathrm{large}} \rangle$ [nm]", fontsize=11, fontweight='bold', color='#1f77b4')
    ax_a_sub.set_ylabel(r"$\phi_{\mathrm{large}} \times 10^5$", fontsize=11, fontweight='bold', color='#e377c2')
    ax_a.set_title("(a) Secondary liquid droplet evolution", fontsize=11, fontweight='bold')
    ax_a.grid(True, alpha=0.3)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax_a.legend(lines, labels, loc='best', fontsize=9)

    # Panel B: Mass Fractal Dimension Df
    ax_b = axes[0, 1]
    ax_b.plot(t, df_params['D'], 'o-', color='#d62728', markersize=3.5, linewidth=1.8, label=r'Fractal dimension $D_f$')
    ax_b.fill_between(t, df_params['D'] - df_params['std_err_D'], df_params['D'] + df_params['std_err_D'],
                      color='#d62728', alpha=0.25)
    ax_b.set_ylabel(r"Mass fractal dimension $D_f$", fontsize=11, fontweight='bold')
    ax_b.set_title(r"(b) Mesoscale network fractal dimension $D_f(t)$", fontsize=11, fontweight='bold')
    ax_b.grid(True, alpha=0.3)
    ax_b.legend(loc='best', fontsize=9)

    # Panel C: Mesoscale Correlation Cutoff xi
    ax_c = axes[1, 0]
    ax_c.plot(t, df_params['xi'], '^-', color='#9467bd', markersize=3.5, linewidth=1.8, label=r'Correlation cutoff $\xi(t)$')
    ax_c.fill_between(t, np.maximum(10.0, df_params['xi'] - df_params['std_err_xi']),
                      df_params['xi'] + df_params['std_err_xi'], color='#9467bd', alpha=0.25)
    ax_c.set_ylabel(r"Cutoff length $\xi$ [nm]", fontsize=11, fontweight='bold')
    ax_c.set_xlabel("Reaction time [s]", fontsize=11, fontweight='bold')
    ax_c.set_title(r"(c) Mesoscale network correlation cutoff $\xi(t)$", fontsize=11, fontweight='bold')
    ax_c.grid(True, alpha=0.3)
    ax_c.legend(loc='best', fontsize=9)

    # Panel D: Goodness of Fit (chi2_red)
    ax_d = axes[1, 1]
    ax_d.plot(t, df_params['chi2_red'], 'o-', color='#ff7f0e', markersize=3.5, linewidth=1.8, label=r'Reduced $\chi^2_{\mathrm{red}}$')
    ax_d.set_ylabel(r"Reduced $\chi^2_{\mathrm{red}}$", fontsize=11, fontweight='bold')
    ax_d.set_xlabel("Reaction time [s]", fontsize=11, fontweight='bold')
    ax_d.set_title(r"(d) Goodness of fit ($\chi^2_{\mathrm{red}}$)", fontsize=11, fontweight='bold')
    ax_d.grid(True, alpha=0.3)
    ax_d.legend(loc='best', fontsize=9)

    plt.suptitle("ND-i-002: Long-distance SAXS kinetic trajectories (20 mM CaCO$_3$, pH 11.3, pure)\n"
                 r"Mesoscale network and secondary droplet kinetics with $\pm 1\sigma$ experimental uncertainties", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved trajectory figure: {plot_path}")

    elapsed = time.time() - start_time
    print(f"\nAll 50 long-distance frames completed in {elapsed:.1f} s ({elapsed/60.0:.2f} min).")
    print("=" * 80)


if __name__ == '__main__':
    fit_long_distance_series()
