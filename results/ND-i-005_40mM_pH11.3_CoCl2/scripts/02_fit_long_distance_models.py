#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Long Distance Time-Resolved Fitting for ND-i-005.
-----------------------------------------------------------
Fits the 50 long-distance SAXS frames (31.0 m geometry, t = 0.060 s to 10.877 s)
capturing the mesoscale mass fractal network and secondary droplet kinetics.

Key Physical Constraints & Workflows:
1. Core Parameter Transfer:
   - Primary building block core dimensions (Rmean, stdDev) transferred from short-distance series.
   - High-q PNC parameters transferred frame-by-frame to accurately model the overlap window.
2. Droplet Population:
   - Polydispersity locked to constant physical ratio sigma/R = 0.143 (from baseline).
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

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_model_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_long_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'long_distance')

CONST_POLY_RATIO = 0.143  # Relative droplet polydispersity sigma_large / R_large


def fit_long_distance_series():
    print("=" * 80)
    print("  LONG-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-005, 50 FRAMES)")
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

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Transferred primary particle dimensions: Rmean = {transferred_Rmean:.2f} nm, stdDev = {transferred_stdDev:.2f} nm")

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
        dy_fit = dy[mask]

        # Transferred short-distance parameters for frame f_idx
        row_short = df_short[df_short['frame_index'] == f_idx].iloc[0]
        a_fluc_f = float(row_short['A_fluc'])
        xi_fluc_f = float(row_short['xi_fluc'])

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)
        p.add('Rmean', value=transferred_Rmean, vary=False)
        p.add('stdDev', value=transferred_stdDev, vary=False)
        p.add('bkg', value=0.0, vary=False)
        p.add('poly_ratio', value=CONST_POLY_RATIO, vary=False)
        p.add('A_fluc', value=a_fluc_f, vary=False)
        p.add('xi_fluc', value=xi_fluc_f, vary=False)
        p.add('mu_fluc', value=4.0, vary=False)

        if prev_params is None:
            p.add('mean_radius_large', value=bp['mean_radius_large'], min=1500.0, max=5500.0, vary=True)
            p.add('std_radius_large', expr=f'mean_radius_large * {CONST_POLY_RATIO}')
            p.add('total_volume_fraction_large', value=bp['total_volume_fraction_large'], min=1e-8, max=1e-3, vary=True)
            p.add('phi', value=bp['phi'], min=1e-6, max=0.01, vary=True)
            p.add('D', value=bp['D'] + 0.15, min=1.5, max=2.6, vary=True)
            p.add('xi', value=bp['xi'], min=100.0, max=15000.0, vary=True)
        else:
            ref = prev_params
            p.add('mean_radius_large', value=ref['mean_radius_large'].value,
                  min=max(1500.0, ref['mean_radius_large'].value - 120.0),
                  max=min(5500.0, ref['mean_radius_large'].value + 120.0), vary=True)
            p.add('std_radius_large', expr=f'mean_radius_large * {CONST_POLY_RATIO}')
            p.add('total_volume_fraction_large', value=ref['total_volume_fraction_large'].value,
                  min=max(1e-8, ref['total_volume_fraction_large'].value * 0.75),
                  max=min(1e-3, ref['total_volume_fraction_large'].value * 1.25), vary=True)
            p.add('phi', value=ref['phi'].value,
                  min=max(1e-6, ref['phi'].value * 0.8),
                  max=min(0.01, ref['phi'].value * 1.2), vary=True)
            p.add('D', value=ref['D'].value,
                  min=max(1.5, ref['D'].value - 0.035),
                  max=min(2.6, ref['D'].value + 0.035), vary=True)
            p.add('xi', value=ref['xi'].value,
                  min=max(100.0, ref['xi'].value * 0.75),
                  max=min(15000.0, ref['xi'].value * 1.25), vary=True)

        def res_func(params):
            calc = calc_model_intensity(params, q_fit, include_large_spheres=True, include_oz=False, integration_steps=80)
            res = (y_fit - calc) / dy_fit
            if prev_params is not None:
                dD = (params['D'].value - prev_params['D'].value) / 0.02
                dR = (params['mean_radius_large'].value - prev_params['mean_radius_large'].value) / 60.0
                dphi = (params['phi'].value - prev_params['phi'].value) / max(1e-6, prev_params['phi'].value * 0.15)
                res = np.append(res, [1.5 * dD, 1.0 * dR, 0.5 * dphi])
            return res

        out = minimize(res_func, p, method='least_squares', max_nfev=350)
        prev_params = out.params

        # Compute exact covariance matrix and standard errors
        p_cov = Parameters()
        p_cov.add('rho2', value=RHO2_CALIBRATED, vary=False)
        p_cov.add('Rmean', value=transferred_Rmean, vary=False)
        p_cov.add('stdDev', value=transferred_stdDev, vary=False)
        p_cov.add('bkg', value=0.0, vary=False)
        p_cov.add('poly_ratio', value=CONST_POLY_RATIO, vary=False)
        p_cov.add('A_fluc', value=a_fluc_f, vary=False)
        p_cov.add('xi_fluc', value=xi_fluc_f, vary=False)
        p_cov.add('mu_fluc', value=4.0, vary=False)

        p_cov.add('mean_radius_large', value=out.params['mean_radius_large'].value, min=1e-8, max=1e5, vary=True)
        p_cov.add('std_radius_large', expr=f'mean_radius_large * {CONST_POLY_RATIO}')
        p_cov.add('total_volume_fraction_large', value=out.params['total_volume_fraction_large'].value, min=1e-8, max=1e5, vary=True)
        p_cov.add('phi', value=out.params['phi'].value, min=1e-8, max=1e5, vary=True)
        p_cov.add('D', value=out.params['D'].value, min=1e-8, max=1e5, vary=True)
        p_cov.add('xi', value=out.params['xi'].value, min=1e-8, max=1e5, vary=True)

        def res_chi2_only(params):
            calc = calc_model_intensity(params, q_fit, include_large_spheres=True, include_oz=False, integration_steps=80)
            return (y_fit - calc) / dy_fit

        out_cov = minimize(res_chi2_only, p_cov, method='leastsq', calc_covar=True)
        final_p = out_cov.params
        redchi = float(out_cov.redchi)

        # Full range calculated intensity
        full_calc = calc_model_intensity(final_p, q, include_large_spheres=True, include_oz=False, integration_steps=80)
        fitted_curves[str(f_idx)] = full_calc

        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'long_file': f_file,
            'chi2_red': redchi,
            'aic': float(out_cov.aic),
            'bic': float(out_cov.bic),
            'Rmean': transferred_Rmean,
            'std_err_Rmean': 0.0,
            'stdDev': transferred_stdDev,
            'std_err_stdDev': 0.0,
            'A_fluc': a_fluc_f,
            'xi_fluc': xi_fluc_f
        }
        for k in ['mean_radius_large', 'total_volume_fraction_large', 'phi', 'D', 'xi']:
            val = float(final_p[k].value)
            err = float(final_p[k].stderr) if final_p[k].stderr is not None else 0.0
            rec[k] = val
            rec[f'std_err_{k}'] = err

        rec['std_radius_large'] = float(final_p['mean_radius_large'].value * CONST_POLY_RATIO)
        rec['std_err_std_radius_large'] = float(rec['std_err_mean_radius_large'] * CONST_POLY_RATIO)

        fitted_records.append(rec)
        if (i + 1) % 10 == 0 or i == 0 or i == n_frames - 1:
            print(f"  Frame {f_idx:2d} (t = {t_sec:6.3f}s): chi2_red = {redchi:.3f}, "
                  f"R_large = {rec['mean_radius_large']:.0f} +/- {rec['std_err_mean_radius_large']:.0f} nm, "
                  f"D = {rec['D']:.3f} +/- {rec['std_err_D']:.3f}, "
                  f"xi = {rec['xi']:.0f} +/- {rec['std_err_xi']:.0f} nm")

    # 3. Save Output Tables
    df_params = pd.DataFrame(fitted_records)
    csv_path = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(csv_path, index=False)
    print(f"\nSaved long-distance parameters CSV: {csv_path}")

    df_curves = pd.DataFrame(fitted_curves)
    curves_path = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_path, index=False)
    print(f"Saved fitted curves CSV: {curves_path}")

    # 4. Generate Publication-Quality Trajectory Figure with Shaded Uncertainties
    print("\nGenerating long-distance trajectory figure with experimental uncertainties...")
    t = df_params['time_seconds'].values

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("ND-i-005 Long Distance (31.0 m): Time-Resolved Kinetic Trajectories\n(40 mM CaCO$_3$, pH 11.3 with Cobalt Additive)", fontsize=13, fontweight='bold')

    # Panel A: Large Droplet Radius
    ax_a = axes[0, 0]
    ax_a.plot(t, df_params['mean_radius_large'], 'o-', color='#1f77b4', markersize=3.5, linewidth=1.8, label=r'Droplet Radius $\langle R_{\mathrm{large}} \rangle$')
    ax_a.fill_between(t, df_params['mean_radius_large'] - df_params['std_err_mean_radius_large'],
                      df_params['mean_radius_large'] + df_params['std_err_mean_radius_large'],
                      color='#1f77b4', alpha=0.25, label=r'$\pm 1\sigma$ Droplet Radius Error')
    ax_a.plot(t, df_params['xi'], '^--', color='#e377c2', markersize=3.5, linewidth=1.5, label=r'Network Cutoff $\xi(t)$')
    ax_a.fill_between(t, np.maximum(10.0, df_params['xi'] - df_params['std_err_xi']),
                      df_params['xi'] + df_params['std_err_xi'],
                      color='#e377c2', alpha=0.18, label=r'$\pm 1\sigma$ Correlation Cutoff')
    ax_a.set_ylabel("Length Scale [nm]", fontsize=11, fontweight='bold')
    ax_a.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_a.set_title(r"(A) Droplet Population & Network Correlation Length", fontsize=11, fontweight='bold')
    ax_a.grid(True, linestyle='--', alpha=0.5)
    ax_a.legend(loc='best', frameon=True, fontsize=9)
    ax_a.set_xscale('log')

    # Panel B: Mass Fractal Dimension Df(t)
    ax_b = axes[0, 1]
    ax_b.plot(t, df_params['D'], 'o-', color='#2ca02c', markersize=3.5, linewidth=1.8, label=r'Fractal Dimension $D_f(t)$')
    ax_b.fill_between(t, df_params['D'] - df_params['std_err_D'], df_params['D'] + df_params['std_err_D'],
                      color='#2ca02c', alpha=0.25, label=r'$\pm 1\sigma$ Confidence Band')
    ax_b.set_ylabel("Mass Fractal Dimension $D_f$", fontsize=11, fontweight='bold')
    ax_b.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_b.set_title("(B) Mesoscale Fractal Network Architecture $D_f(t)$", fontsize=11, fontweight='bold')
    ax_b.grid(True, linestyle='--', alpha=0.5)
    ax_b.legend(loc='best', frameon=True, fontsize=9)
    ax_b.set_ylim(1.9, 2.5)
    ax_b.set_xscale('log')

    # Panel C: Volume Fractions
    ax_c = axes[1, 0]
    ax_c.plot(t, df_params['phi'] * 1e4, 'o-', color='#9467bd', markersize=3.5, linewidth=1.8, label=r'Primary Fractal $\phi_{\mathrm{prim}} \times 10^4$')
    ax_c.fill_between(t, (df_params['phi'] - df_params['std_err_phi']) * 1e4, (df_params['phi'] + df_params['std_err_phi']) * 1e4,
                      color='#9467bd', alpha=0.25)
    ax_c.plot(t, df_params['total_volume_fraction_large'] * 1e5, 's-', color='#d62728', markersize=3.5, linewidth=1.8, label=r'Droplet Fraction $\phi_{\mathrm{large}} \times 10^5$')
    ax_c.fill_between(t, (df_params['total_volume_fraction_large'] - df_params['std_err_total_volume_fraction_large']) * 1e5,
                      (df_params['total_volume_fraction_large'] + df_params['std_err_total_volume_fraction_large']) * 1e5,
                      color='#d62728', alpha=0.25)
    ax_c.set_ylabel("Volume Fraction", fontsize=11, fontweight='bold')
    ax_c.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_c.set_title("(C) Phase Volume Fraction Evolution", fontsize=11, fontweight='bold')
    ax_c.grid(True, linestyle='--', alpha=0.5)
    ax_c.legend(loc='best', frameon=True, fontsize=9)
    ax_c.set_xscale('log')

    # Panel D: Goodness of Fit
    ax_d = axes[1, 1]
    ax_d.plot(t, df_params['chi2_red'], 'o-', color='#333333', markersize=3.5, linewidth=1.5, label=r'$\chi^2_{\mathrm{red}}(t)$')
    ax_d.axhline(1.0, color='red', linestyle='--', linewidth=1.0, alpha=0.7, label=r'Ideal $\chi^2_{\mathrm{red}} = 1.0$')
    ax_d.set_ylabel(r"Reduced $\chi^2$", fontsize=11, fontweight='bold')
    ax_d.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_d.set_title(r"(D) Goodness-of-Fit Metric $\chi^2_{\mathrm{red}}$", fontsize=11, fontweight='bold')
    ax_d.grid(True, linestyle='--', alpha=0.5)
    ax_d.legend(loc='best', frameon=True, fontsize=9)
    ax_d.set_ylim(0.0, 1.5)
    ax_d.set_xscale('log')

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved trajectory graphic: {plot_path}")

    elapsed = time.time() - start_time
    print(f"\nLong-distance series fitting completed in {elapsed:.1f} s.")
    return df_params


if __name__ == '__main__':
    fit_long_distance_series()
