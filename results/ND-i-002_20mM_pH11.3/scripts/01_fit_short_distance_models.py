#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 01: Short Distance Time-Resolved Fitting for ND-i-002.
------------------------------------------------------------
Fits the 50 short-distance SAXS frames (3.0 m geometry, t = 0.060 s to 10.877 s)
using a two-step physical procedure:
1. Fit the polydisperse spherical building blocks (Rmean, stdDev, phi, D) in the first place
   on the lower q-range (q <= 0.40 nm^-1) where primary particles dominate.
2. Explain the remaining intensity in the fit with the viscoelastic concentration fluctuation
   component (A_fluc, xi_fluc, mu_fluc = 4.0).
3. Jointly refine all parameters across the full range with sequential continuity regularization.
4. Compute rigorous parameter uncertainties (+/- 1 sigma) via the covariance matrix.

Outputs:
  output/short_distance/parameters.csv
  output/short_distance/fitted_curves.csv
  output/short_distance/trajectories.png
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
    decompose_model_components,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi002_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'short_distance')


def fit_short_distance_series():
    print("=" * 80)
    print("  SHORT-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-002, 50 FRAMES)")
    print("  Model: Polydisperse Primary Spheres + Teixeira Fractal + Viscoelastic Fluctuations")
    print("=" * 80)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load Data & Metadata
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    with open(BASELINE_JSON, 'r') as f:
        bp = json.load(f)

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Transferring baseline correlation length xi = {bp['xi']:.1f} nm")

    # 2. Sequential Fitting Loop
    fitted_records = []
    fitted_curves = {'q': q}
    prev_params = None

    for i in range(n_frames):
        row = df_meta.iloc[i]
        f_idx = int(row['frame_index'])
        t_sec = float(row['time_seconds'])
        f_file = str(row['short_file'])

        y = df_data[str(f_idx)].values
        dy = np.maximum(df_data['error_' + str(f_idx)].values, 1e-4)

        # Valid fitting mask (q >= 0.022 nm^-1, positive intensity)
        mask = (q >= 0.022) & (y > 0)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = dy[mask]

        # STEP 1: Fit polydisperse spherical building blocks on q <= 0.40 nm^-1
        s1_mask = (q_fit <= 0.40)
        q_s1 = q_fit[s1_mask]
        y_s1 = y_fit[s1_mask]
        dy_s1 = dy_fit[s1_mask]

        p1 = Parameters()
        p1.add('rho2', value=RHO2_CALIBRATED, vary=False)
        p1.add('xi', value=bp['xi'], vary=False)
        p1.add('total_volume_fraction_large', value=0.0, vary=False)
        p1.add('A_fluc', value=0.0, vary=False)
        p1.add('xi_fluc', value=2.5, vary=False)
        p1.add('mu_fluc', value=4.0, vary=False)

        if prev_params is None:
            p1.add('phi', value=bp.get('phi', 0.00018), min=1e-6, max=0.05, vary=True)
            p1.add('Rmean', value=bp.get('Rmean', 10.0), min=1.5, max=30.0, vary=True)
            p1.add('stdDev', value=bp.get('stdDev', 4.0), min=0.3, max=15.0, vary=True)
            p1.add('D', value=bp.get('D', 1.93), min=1.3, max=2.8, vary=True)
            p1.add('bkg', value=bp.get('bkg', 0.008), min=1e-5, max=0.2, vary=True)
        else:
            p1.add('phi', value=prev_params['phi'].value, min=prev_params['phi'].value * 0.7, max=prev_params['phi'].value * 1.3, vary=True)
            p1.add('Rmean', value=prev_params['Rmean'].value, min=max(1.5, prev_params['Rmean'].value - 1.0), max=prev_params['Rmean'].value + 1.0, vary=True)
            p1.add('stdDev', value=prev_params['stdDev'].value, min=max(0.3, prev_params['stdDev'].value - 0.8), max=prev_params['stdDev'].value + 0.8, vary=True)
            p1.add('D', value=prev_params['D'].value, min=max(1.3, prev_params['D'].value - 0.05), max=min(2.8, prev_params['D'].value + 0.05), vary=True)
            p1.add('bkg', value=prev_params['bkg'].value, min=max(1e-5, prev_params['bkg'].value * 0.7), max=prev_params['bkg'].value * 1.3, vary=True)

        def res_s1(pars):
            c = calc_model_intensity(pars, q_s1, include_large_spheres=False, include_oz=False, integration_steps=80)
            return (y_s1 - c) / dy_s1

        out_s1 = minimize(res_s1, p1, method='least_squares')

        # STEP 2: Explain remaining intensity with viscoelastic fluctuation component
        calc_s1_full = calc_model_intensity(out_s1.params, q_fit, include_large_spheres=False, include_oz=False, integration_steps=80)
        excess = y_fit - calc_s1_full
        s2_mask = (q_fit > 0.35)

        p2 = out_s1.params.copy()
        for k in ['phi', 'Rmean', 'stdDev', 'D', 'bkg']:
            p2[k].vary = False

        init_a_fluc = float(np.clip(np.median(excess[s2_mask]) if np.any(s2_mask) else 0.05, 0.001, 1.5))
        if prev_params is not None and 'A_fluc' in prev_params:
            init_a_fluc = prev_params['A_fluc'].value

        p2['A_fluc'].value = init_a_fluc
        p2['A_fluc'].vary = True
        p2['A_fluc'].min = 1e-4
        p2['A_fluc'].max = 3.0

        p2['xi_fluc'].value = 2.5 if prev_params is None else prev_params['xi_fluc'].value
        p2['xi_fluc'].vary = True
        p2['xi_fluc'].min = 0.5
        p2['xi_fluc'].max = 15.0

        def res_s2(pars):
            c = calc_model_intensity(pars, q_fit, include_large_spheres=False, include_oz=True, integration_steps=80)
            return (y_fit - c) / dy_fit

        out_s2 = minimize(res_s2, p2, method='least_squares')

        # STEP 3: Joint Refinement with Sequential Regularization
        p3 = out_s2.params.copy()
        p3['phi'].vary = True
        p3['Rmean'].vary = True
        p3['stdDev'].vary = True
        p3['D'].vary = True
        p3['bkg'].vary = True
        p3['A_fluc'].vary = True
        p3['xi_fluc'].vary = True

        def res_joint(pars):
            c = calc_model_intensity(pars, q_fit, include_large_spheres=False, include_oz=True, integration_steps=90)
            residuals = (y_fit - c) / dy_fit

            reg = 0.0
            if prev_params is not None:
                dD = (pars['D'].value - prev_params['D'].value) / 0.008
                dR = (pars['Rmean'].value - prev_params['Rmean'].value) / 0.15
                dxi_f = (pars['xi_fluc'].value - prev_params['xi_fluc'].value) / 0.20
                reg = np.array([dD, dR, dxi_f])
                return np.concatenate([residuals, reg])
            return residuals

        out_joint = minimize(res_joint, p3, method='least_squares')

        # STEP 4: Exact Covariance Matrix Computation (+/- 1 sigma)
        p_cov = out_joint.params.copy()
        for k in ['phi', 'Rmean', 'stdDev', 'D', 'bkg', 'A_fluc', 'xi_fluc']:
            p_cov[k].vary = True
        out_cov = minimize(res_joint, p_cov, method='leastsq', calc_covar=True)
        final_p = out_cov.params

        # Calculate final model curve & chi2
        fit_y = calc_model_intensity(final_p, q, include_large_spheres=False, include_oz=True, integration_steps=120)
        fit_y_v = fit_y[mask]
        chi2_red = float(np.sum(((y_fit - fit_y_v) / dy_fit) ** 2) / max(1, len(y_fit) - 7))

        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'filename': f_file,
            'chi2_red': chi2_red,
            'phi': final_p['phi'].value,
            'std_err_phi': float(final_p['phi'].stderr) if final_p['phi'].stderr is not None else 0.0,
            'Rmean': final_p['Rmean'].value,
            'std_err_Rmean': float(final_p['Rmean'].stderr) if final_p['Rmean'].stderr is not None else 0.0,
            'stdDev': final_p['stdDev'].value,
            'std_err_stdDev': float(final_p['stdDev'].stderr) if final_p['stdDev'].stderr is not None else 0.0,
            'D': final_p['D'].value,
            'std_err_D': float(final_p['D'].stderr) if final_p['D'].stderr is not None else 0.0,
            'xi': final_p['xi'].value,
            'A_fluc': final_p['A_fluc'].value,
            'std_err_A_fluc': float(final_p['A_fluc'].stderr) if final_p['A_fluc'].stderr is not None else 0.0,
            'xi_fluc': final_p['xi_fluc'].value,
            'std_err_xi_fluc': float(final_p['xi_fluc'].stderr) if final_p['xi_fluc'].stderr is not None else 0.0,
            'mu_fluc': final_p['mu_fluc'].value,
            'bkg': final_p['bkg'].value,
            'std_err_bkg': float(final_p['bkg'].stderr) if final_p['bkg'].stderr is not None else 0.0
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_y
        prev_params = final_p

        if f_idx % 10 == 0 or f_idx == n_frames - 1:
            print(f"  Frame {f_idx:02d} (t = {t_sec:.3f} s): chi2 = {chi2_red:.3f} | Rmean = {final_p['Rmean'].value:.2f} nm, D = {final_p['D'].value:.3f}, xi_fluc = {final_p['xi_fluc'].value:.2f} nm")

    # 3. Save Output Tables
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)
    print(f"\nSaved short-distance parameters: {params_csv}")

    df_curves = pd.DataFrame(fitted_curves)
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)
    print(f"Saved short-distance fitted curves: {curves_csv}")

    # 4. Generate 4-Panel Trajectory Figure
    print("\nGenerating short-distance parameter trajectories figure...")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), sharex=True)
    t = df_params['time_seconds']

    # Panel A: Primary Building Block Core & Polydispersity
    ax_a = axes[0, 0]
    ax_a.plot(t, df_params['Rmean'], 'o-', color='#1f77b4', markersize=3.5, linewidth=1.8, label=r'Primary Radius $\langle R_{\mathrm{prim}} \rangle$')
    ax_a.fill_between(t, df_params['Rmean'] - df_params['std_err_Rmean'], df_params['Rmean'] + df_params['std_err_Rmean'],
                      color='#1f77b4', alpha=0.25)
    ax_a.plot(t, df_params['stdDev'], 's--', color='#2ca02c', markersize=3.5, linewidth=1.8, label=r'Dispersion $\sigma_{\mathrm{prim}}$')
    ax_a.fill_between(t, df_params['stdDev'] - df_params['std_err_stdDev'], df_params['stdDev'] + df_params['std_err_stdDev'],
                      color='#2ca02c', alpha=0.25)
    ax_a.set_ylabel("Length Scale [nm]", fontsize=11, fontweight='bold')
    ax_a.set_title("(a) Primary Building Block Core Dimensions", fontsize=11, fontweight='bold')
    ax_a.grid(True, alpha=0.3)
    ax_a.legend(loc='best', fontsize=9)

    # Panel B: Mass Fractal Dimension Df and Correlation xi_fluc
    ax_b = axes[0, 1]
    ax_b_sub = ax_b.twinx()
    l1 = ax_b.plot(t, df_params['D'], 'o-', color='#d62728', markersize=3.5, linewidth=1.8, label=r'Fractal Dimension $D_f$')
    ax_b.fill_between(t, df_params['D'] - df_params['std_err_D'], df_params['D'] + df_params['std_err_D'],
                      color='#d62728', alpha=0.25)
    l2 = ax_b_sub.plot(t, df_params['xi_fluc'], 's-', color='#9467bd', markersize=3.5, linewidth=1.8, label=r'Fluctuation Correlation $\xi_{\mathrm{fluc}}$')
    ax_b_sub.fill_between(t, df_params['xi_fluc'] - df_params['std_err_xi_fluc'], df_params['xi_fluc'] + df_params['std_err_xi_fluc'],
                          color='#9467bd', alpha=0.25)
    ax_b.set_ylabel(r"Mass Fractal Dimension $D_f$", fontsize=11, fontweight='bold', color='#d62728')
    ax_b_sub.set_ylabel(r"Correlation Length $\xi_{\mathrm{fluc}}$ [nm]", fontsize=11, fontweight='bold', color='#9467bd')
    ax_b.set_title(r"(b) Network Fractal Dimension $D_f$ & Viscoelastic $\xi_{\mathrm{fluc}}$", fontsize=11, fontweight='bold')
    ax_b.grid(True, alpha=0.3)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax_b.legend(lines, labels, loc='best', fontsize=9)

    # Panel C: Viscoelastic Amplitude A_fluc and Volume Fraction phi
    ax_c = axes[1, 0]
    ax_c_sub = ax_c.twinx()
    l1 = ax_c.plot(t, df_params['A_fluc'], 'o-', color='#e377c2', markersize=3.5, linewidth=1.8, label=r'Fluctuation Amplitude $A_{\mathrm{fluc}}$')
    ax_c.fill_between(t, df_params['A_fluc'] - df_params['std_err_A_fluc'], df_params['A_fluc'] + df_params['std_err_A_fluc'],
                      color='#e377c2', alpha=0.25)
    l2 = ax_c_sub.plot(t, df_params['phi'] * 1e4, '^-', color='#8c564b', markersize=3.5, linewidth=1.8, label=r'Primary Vol. Frac. $\phi \times 10^4$')
    ax_c_sub.fill_between(t, (df_params['phi'] - df_params['std_err_phi']) * 1e4, (df_params['phi'] + df_params['std_err_phi']) * 1e4,
                          color='#8c564b', alpha=0.25)
    ax_c.set_ylabel(r"$A_{\mathrm{fluc}}$ [cm$^{-1}$]", fontsize=11, fontweight='bold', color='#e377c2')
    ax_c_sub.set_ylabel(r"Primary $\phi \times 10^4$", fontsize=11, fontweight='bold', color='#8c564b')
    ax_c.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_c.set_title("(c) Viscoelastic Amplitude & Building Block Volume Fraction", fontsize=11, fontweight='bold')
    ax_c.grid(True, alpha=0.3)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax_c.legend(lines, labels, loc='best', fontsize=9)

    # Panel D: Goodness of Fit (chi2_red)
    ax_d = axes[1, 1]
    ax_d.plot(t, df_params['chi2_red'], 'o-', color='#ff7f0e', markersize=3.5, linewidth=1.8, label=r'Reduced $\chi^2_{\mathrm{red}}$')
    ax_d.set_ylabel(r"Reduced $\chi^2_{\mathrm{red}}$", fontsize=11, fontweight='bold')
    ax_d.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_d.set_title(r"(d) Goodness of Fit ($\chi^2_{\mathrm{red}}$)", fontsize=11, fontweight='bold')
    ax_d.grid(True, alpha=0.3)
    ax_d.legend(loc='best', fontsize=9)

    plt.suptitle("ND-i-002 Short-Distance SAXS Kinetic Trajectories (20 mM CaCO$_3$, pH 11.3, Pure)\n"
                 r"Point-to-point evolutionary model with $\pm 1\sigma$ experimental uncertainties", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved trajectory figure: {plot_path}")

    elapsed = time.time() - start_time
    print(f"\nAll 50 short-distance frames completed in {elapsed:.1f} s ({elapsed/60.0:.2f} min).")
    print("=" * 80)


if __name__ == '__main__':
    fit_short_distance_series()
