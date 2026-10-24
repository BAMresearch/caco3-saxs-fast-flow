#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 01: Short Distance Time-Resolved Fitting for ND-i-005.
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

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_model_intensity,
    decompose_model_components,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi005_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'short_distance')


def fit_short_distance_series():
    print("=" * 80)
    print("  SHORT-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-005, 50 FRAMES)")
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

        # STEP 1: Fit polydisperse spherical building blocks in the first place on q <= 0.40 nm^-1
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
            p1.add('phi', value=bp['phi'], min=1e-6, max=0.01, vary=True)
            p1.add('Rmean', value=bp['Rmean'], min=1.5, max=15.0, vary=True)
            p1.add('stdDev', value=bp['stdDev'], min=0.3, max=8.0, vary=True)
            p1.add('D', value=bp['D'], min=1.4, max=2.6, vary=True)
            p1.add('bkg', value=bp['bkg'], min=1e-5, max=0.05, vary=True)
        else:
            ref = prev_params
            p1.add('phi', value=float(ref['phi'].value), min=1e-6, max=0.01, vary=True)
            p1.add('Rmean', value=float(ref['Rmean'].value), min=1.5, max=15.0, vary=True)
            p1.add('stdDev', value=float(ref['stdDev'].value), min=0.3, max=8.0, vary=True)
            p1.add('D', value=float(ref['D'].value), min=1.4, max=2.6, vary=True)
            p1.add('bkg', value=float(ref['bkg'].value), min=1e-5, max=0.05, vary=True)

        def res_s1(params):
            calc = calc_model_intensity(params, q_s1, include_large_spheres=False, include_oz=True, integration_steps=80)
            res = (y_s1 - calc) / dy_s1
            if prev_params is not None:
                dD = (params['D'].value - prev_params['D'].value) / 0.02
                dR = (params['Rmean'].value - prev_params['Rmean'].value) / 0.5
                res = np.append(res, [1.5 * dD, 0.8 * dR])
            return res

        out_s1 = minimize(res_s1, p1, method='least_squares', max_nfev=250)

        # STEP 2: Explain remaining intensity with viscoelastic fluctuation component
        p2 = out_s1.params.copy()
        p2['phi'].vary = False
        p2['Rmean'].vary = False
        p2['stdDev'].vary = False
        p2['D'].vary = False
        p2['bkg'].vary = True

        init_A = float(prev_params['A_fluc'].value) if prev_params is not None else 0.5
        init_xi_fluc = float(prev_params['xi_fluc'].value) if prev_params is not None else 8.0
        p2['A_fluc'].set(value=init_A, min=1e-6, max=5.0, vary=True)
        p2['xi_fluc'].set(value=init_xi_fluc, min=0.5, max=25.0, vary=True)
        p2['mu_fluc'].set(value=4.0, vary=False)

        def res_s2(params):
            calc = calc_model_intensity(params, q_fit, include_large_spheres=False, include_oz=True, integration_steps=80)
            return (y_fit - calc) / dy_fit

        out_s2 = minimize(res_s2, p2, method='least_squares', max_nfev=250)

        # STEP 3: Joint refinement of all parameters across full q-range with sequential continuity
        p3 = out_s2.params.copy()
        p3['phi'].set(min=1e-6, max=0.01, vary=True)
        p3['Rmean'].set(min=1.5, max=15.0, vary=True)
        p3['stdDev'].set(min=0.3, max=8.0, vary=True)
        p3['D'].set(min=1.4, max=2.6, vary=True)
        p3['bkg'].set(min=1e-5, max=0.05, vary=True)
        p3['A_fluc'].set(min=1e-6, max=5.0, vary=True)
        p3['xi_fluc'].set(min=0.5, max=25.0, vary=True)
        p3['mu_fluc'].set(value=4.0, vary=False)

        def res_joint(params):
            calc = calc_model_intensity(params, q_fit, include_large_spheres=False, include_oz=True, integration_steps=80)
            res = (y_fit - calc) / dy_fit
            if prev_params is not None:
                dD = (params['D'].value - prev_params['D'].value) / 0.02
                dR = (params['Rmean'].value - prev_params['Rmean'].value) / 0.5
                dA = (params['A_fluc'].value - prev_params['A_fluc'].value) / max(0.03, prev_params['A_fluc'].value * 0.2)
                dxi = (params['xi_fluc'].value - prev_params['xi_fluc'].value) / 0.8
                res = np.append(res, [1.5 * dD, 0.8 * dR, 0.8 * dA, 0.8 * dxi])
            return res

        out_joint = minimize(res_joint, p3, method='least_squares', max_nfev=350)
        prev_params = out_joint.params

        # STEP 4: Exact Covariance Matrix and Standard Errors
        p_cov = Parameters()
        p_cov.add('rho2', value=RHO2_CALIBRATED, vary=False)
        p_cov.add('xi', value=bp['xi'], vary=False)
        p_cov.add('total_volume_fraction_large', value=0.0, vary=False)
        p_cov.add('mu_fluc', value=4.0, vary=False)
        for k in ['phi', 'Rmean', 'stdDev', 'D', 'bkg', 'A_fluc', 'xi_fluc']:
            p_cov.add(k, value=out_joint.params[k].value, min=1e-8, max=1e5, vary=True)

        def res_chi2_only(params):
            calc = calc_model_intensity(params, q_fit, include_large_spheres=False, include_oz=True, integration_steps=80)
            return (y_fit - calc) / dy_fit

        out_cov = minimize(res_chi2_only, p_cov, method='leastsq', calc_covar=True)
        final_p = out_cov.params
        redchi = float(out_cov.redchi)

        # Full-range fitted curve
        full_calc = calc_model_intensity(final_p, q, include_large_spheres=False, include_oz=True, integration_steps=80)
        fitted_curves[str(f_idx)] = full_calc

        # Collect record
        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'short_file': f_file,
            'chi2_red': redchi,
            'aic': float(out_cov.aic),
            'bic': float(out_cov.bic),
            'xi': float(bp['xi']),
            'std_err_xi': 0.0,
            'mu_fluc': 4.0,
            'std_err_mu_fluc': 0.0
        }
        for k in ['phi', 'Rmean', 'stdDev', 'D', 'bkg', 'A_fluc', 'xi_fluc']:
            val = float(final_p[k].value)
            err = float(final_p[k].stderr) if final_p[k].stderr is not None else 0.0
            rec[k] = val
            rec[f'std_err_{k}'] = err

        fitted_records.append(rec)
        if (i + 1) % 10 == 0 or i == 0 or i == n_frames - 1:
            print(f"  Frame {f_idx:2d} (t = {t_sec:6.3f}s): chi2_red = {redchi:.3f}, "
                  f"Rmean = {rec['Rmean']:.2f} +/- {rec['std_err_Rmean']:.2f} nm, "
                  f"D = {rec['D']:.3f} +/- {rec['std_err_D']:.3f}, "
                  f"A_fluc = {rec['A_fluc']:.4f} +/- {rec['std_err_A_fluc']:.4f}, "
                  f"xi_fluc = {rec['xi_fluc']:.2f} +/- {rec['std_err_xi_fluc']:.2f} nm")

    # 3. Save Results
    df_params = pd.DataFrame(fitted_records)
    csv_path = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(csv_path, index=False)
    print(f"\nSaved short-distance parameters CSV: {csv_path}")

    df_curves = pd.DataFrame(fitted_curves)
    curves_path = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_path, index=False)
    print(f"Saved fitted curves CSV: {curves_path}")

    # 4. Generate Publication-Quality Trajectory Figure with Shaded Uncertainties
    print("\nGenerating short-distance trajectory figure with experimental uncertainties...")
    t = df_params['time_seconds'].values

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("ND-i-005 Short Distance (3.0 m): Time-Resolved Kinetic Trajectories\n(40 mM CaCO$_3$, pH 11.3 with Cobalt Additive : Viscoelastic Fluctuation Model)", fontsize=13, fontweight='bold')

    # Panel A: Length Scales (Rmean, stdDev, xi_fluc)
    ax_a = axes[0, 0]
    ax_a.plot(t, df_params['Rmean'], 'o-', color='#1f77b4', markersize=3.5, linewidth=1.8, label=r'Primary Core $\langle R_{\mathrm{prim}} \rangle$')
    ax_a.fill_between(t, df_params['Rmean'] - df_params['std_err_Rmean'], df_params['Rmean'] + df_params['std_err_Rmean'],
                      color='#1f77b4', alpha=0.25, label=r'$\pm 1\sigma$ Core Error')
    ax_a.plot(t, df_params['stdDev'], '--', color='#17becf', linewidth=1.5, label=r'Core Polydispersity $\sigma_{\mathrm{prim}}$')
    ax_a.plot(t, df_params['xi_fluc'], 's-', color='#d62728', markersize=3.5, linewidth=1.8, label=r'Viscoelastic Fluctuation Mesh $\xi_{\mathrm{fluc}}$')
    ax_a.fill_between(t, df_params['xi_fluc'] - df_params['std_err_xi_fluc'], df_params['xi_fluc'] + df_params['std_err_xi_fluc'],
                      color='#d62728', alpha=0.25, label=r'$\pm 1\sigma$ Fluctuation Error')
    ax_a.set_ylabel("Length Scale [nm]", fontsize=11, fontweight='bold')
    ax_a.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_a.set_title("(a) Primary Core & Fluctuation Length Scales", fontsize=11, fontweight='bold')
    ax_a.grid(True, alpha=0.3)
    ax_a.legend(loc='best', fontsize=9)

    # Panel B: Mass Fractal Dimension D
    ax_b = axes[0, 1]
    ax_b.plot(t, df_params['D'], 'o-', color='#2ca02c', markersize=3.5, linewidth=1.8, label=r'Fractal Dimension $D_f$')
    ax_b.fill_between(t, df_params['D'] - df_params['std_err_D'], df_params['D'] + df_params['std_err_D'],
                      color='#2ca02c', alpha=0.25, label=r'$\pm 1\sigma$ Fractal Dimension Error')
    ax_b.set_ylabel(r"Mass Fractal Dimension $D_f$", fontsize=11, fontweight='bold')
    ax_b.set_xlabel("Reaction Time [s]", fontsize=11, fontweight='bold')
    ax_b.set_title(r"(b) Primary Network Mass Fractal Dimension $D_f$", fontsize=11, fontweight='bold')
    ax_b.grid(True, alpha=0.3)
    ax_b.legend(loc='best', fontsize=9)

    # Panel C: Viscoelastic Amplitude A_fluc and Volume Fraction phi
    ax_c = axes[1, 0]
    ax_c_sub = ax_c.twinx()
    l1 = ax_c.plot(t, df_params['A_fluc'], 'o-', color='#e377c2', markersize=3.5, linewidth=1.8, label=r'Fluctuation Amplitude $A_{\mathrm{fluc}}$')
    ax_c.fill_between(t, df_params['A_fluc'] - df_params['std_err_A_fluc'], df_params['A_fluc'] + df_params['std_err_A_fluc'],
                      color='#e377c2', alpha=0.25)
    l2 = ax_c_sub.plot(t, df_params['phi'] * 1e4, '^-', color='#9467bd', markersize=3.5, linewidth=1.8, label=r'Primary Vol. Frac. $\phi \times 10^4$')
    ax_c_sub.fill_between(t, (df_params['phi'] - df_params['std_err_phi']) * 1e4, (df_params['phi'] + df_params['std_err_phi']) * 1e4,
                          color='#9467bd', alpha=0.25)
    ax_c.set_ylabel(r"$A_{\mathrm{fluc}}$ [cm$^{-1}$]", fontsize=11, fontweight='bold', color='#e377c2')
    ax_c_sub.set_ylabel(r"Primary $\phi \times 10^4$", fontsize=11, fontweight='bold', color='#9467bd')
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
