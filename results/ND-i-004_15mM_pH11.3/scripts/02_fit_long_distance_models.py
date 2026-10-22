#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Long Distance Time-Resolved Fitting for ND-i-004.
-----------------------------------------------------------
Fits the 50 long-distance SAXS frames (31.0 m camera, q in [0.002, 0.22] nm^-1, t = 0.060 s to 10.877 s).

Scientific Principle for 15.0 mM CaCO3 (pH 11.3, Hybrid Regime):
- At 15 mM, spinodal decomposition / mass-fractal percolation and discrete ACC nanosphere growth co-occur.
- Mesoscale network: Teixeira mass fractal S(q, Dm, xi) captures the intermediate-to-low q power-law slope (Dm ~ 1.75 - 2.2).
- Macroscopic domain coarsening: Asymptotic Porod scattering (I ~ A_porod / q^4) for giant domains where q*R >> 1.
- Cross-camera transfer: ACC nanosphere core (R_acc, std_acc, phi_acc) transferred directly from Short Distance,
  providing seamless multi-scale consistency.
- Sequential continuity regularization prevents unphysical parameter jumps.
- Parameter uncertainties (+/- 1 sigma) computed via the covariance matrix.

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
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_long_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'long_distance')


def fit_long_distance_series():
    print("=" * 80, flush=True)
    print("  LONG-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-004, 50 FRAMES)", flush=True)
    print("  Model: Spinodal Mass Fractal Network + Coarsened Porod Domain + Transferred ACC Nanospheres", flush=True)
    print("=" * 80, flush=True)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load Data, Metadata & Transferred Core
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    df_short = pd.read_csv(SHORT_CSV)

    with open(BASELINE_JSON, 'r') as f:
        bp = json.load(f)

    d_base = bp.get('D', 1.79)

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)", flush=True)
    print(f"Baseline anchor: D = {d_base:.3f}", flush=True)

    # 2. Sequential Backward Fitting Loop (49 down to 0)
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

        # Valid long-distance fitting range: q in [0.0025, 0.22] nm^-1
        mask = (q >= 0.0025) & (q <= 0.22) & (y > 0)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = dy[mask]

        # Transferred short-distance parameters for frame f_idx
        row_short = df_short[df_short['frame_index'] == f_idx].iloc[0]
        r_acc_f = float(row_short['R_acc'])
        std_acc_f = float(row_short['std_acc'])
        phi_acc_f = float(row_short['phi_acc'])
        a_porod_f = float(row_short.get('A_porod', 1e-7))

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)

        # Primary nodes in the network
        p.add('R_prim', value=3.5, vary=False)
        p.add('std_prim', value=1.5, vary=False)

        # Transferred ACC nanospheres (fixed from short distance)
        p.add('R_acc', value=r_acc_f, vary=False)
        p.add('std_acc', value=std_acc_f, vary=False)
        p.add('phi_acc', value=phi_acc_f, min=max(1e-8, phi_acc_f * 0.5), max=max(1e-4, phi_acc_f * 2.0), vary=True)

        # Network Fractal Parameters
        if prev_params is None:
            # Anchor at Frame 49: mature percolated network
            p.add('phi_net', value=3e-4, min=1e-6, max=0.01, vary=True)
            p.add('Dm', value=2.25, min=1.5, max=2.6, vary=True)
            p.add('xi', value=200.0, min=20.0, max=2000.0, vary=True)
            p.add('A_porod', value=max(1e-8, a_porod_f * 0.01), min=1e-9, max=0.01, vary=True)
            p.add('bkg', value=1e-4, min=1e-6, max=0.01, vary=True)
        else:
            dm_val = prev_params['Dm'].value
            xi_val = prev_params['xi'].value
            phi_net_val = prev_params['phi_net'].value

            p.add('phi_net', value=phi_net_val, min=1e-6, max=max(1e-3, phi_net_val * 2.0), vary=True)
            p.add('Dm', value=dm_val, min=max(1.4, dm_val - 0.2), max=min(2.6, dm_val + 0.1), vary=True)
            p.add('xi', value=xi_val, min=max(20.0, xi_val * 0.4), max=min(2000.0, xi_val * 1.5), vary=True)
            p.add('A_porod', value=prev_params['A_porod'].value, min=1e-9, max=0.01, vary=True)
            p.add('bkg', value=prev_params['bkg'].value, min=1e-6, max=0.01, vary=True)

        def model_ld(pars, q_eval):
            I_frac = calc_primary_fractal_intensity(q_eval, pars['R_prim'].value, pars['std_prim'].value,
                                                    pars['phi_net'].value, pars['Dm'].value, pars['xi'].value, integration_steps=60)
            I_acc = calc_large_sphere_intensity(q_eval, pars['R_acc'].value, pars['std_acc'].value,
                                                 pars['phi_acc'].value, integration_steps=60)
            I_porod = calc_porod_intensity(q_eval, pars['A_porod'].value, 4.0)
            return I_frac + I_acc + I_porod + pars['bkg'].value

        def res_long(pars):
            c = model_ld(pars, q_fit)
            residuals = (y_fit - c) / dy_fit

            if prev_params is not None:
                dDm = (pars['Dm'].value - prev_params['Dm'].value) / 0.08
                dxi = (np.log(pars['xi'].value) - np.log(prev_params['xi'].value)) / 0.25
                return np.concatenate([residuals, [dDm, dxi]])
            return residuals

        out_fit = minimize(res_long, p, method='least_squares', max_nfev=200)

        # Covariance matrix calculation
        p_cov = out_fit.params.copy()
        for k in ['phi_net', 'Dm', 'xi', 'A_porod', 'phi_acc', 'bkg']:
            p_cov[k].vary = True
        try:
            out_cov = minimize(res_long, p_cov, method='leastsq', calc_covar=True)
            final_p = out_cov.params
        except Exception:
            final_p = out_fit.params

        # Calculate full curve & chi2
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
            'A_porod': final_p['A_porod'].value,
            'std_err_A_porod': float(final_p['A_porod'].stderr) if final_p['A_porod'].stderr is not None else 0.0,
            'phi_acc': final_p['phi_acc'].value,
            'std_err_phi_acc': float(final_p['phi_acc'].stderr) if final_p['phi_acc'].stderr is not None else 0.0,
            'R_acc_transferred': r_acc_f,
            'std_acc_transferred': std_acc_f,
            'bkg': final_p['bkg'].value,
            'std_err_bkg': float(final_p['bkg'].stderr) if final_p['bkg'].stderr is not None else 0.0,
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_full
        prev_params = final_p

        if f_idx % 5 == 0 or f_idx == n_frames - 1 or f_idx == 0:
            print(f"  Frame {f_idx:02d} (t = {t_sec:6.3f} s): chi2 = {chi2_red:.3f} | Dm = {final_p['Dm'].value:.3f} ± {rec['std_err_Dm']:.3f}, xi = {final_p['xi'].value:.0f} nm, A_porod = {final_p['A_porod'].value:.2e}", flush=True)

    # 3. Sort chronologically from frame 0 to 49
    fitted_records.sort(key=lambda x: x['frame_index'])
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)
    print(f"\nSaved long-distance parameters: {params_csv}", flush=True)

    # Reorder columns of fitted_curves so 'q' is first, followed by '0', '1', ..., '49'
    curve_cols = ['q'] + [str(i) for i in range(n_frames)]
    df_curves = pd.DataFrame(fitted_curves)[curve_cols]
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)
    print(f"Saved long-distance fitted curves: {curves_csv}", flush=True)

    # 4. Generate Trajectory Figure
    print("Generating long-distance trajectory plots...", flush=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    t = df_params['time_seconds'].values

    # Panel A: Fractal Dimension Dm
    ax = axes[0, 0]
    ax.errorbar(t, df_params['Dm'], yerr=df_params['std_err_Dm'], fmt='o-', color='#d97706', label=r'Fractal Dimension $D_m$')
    ax.axhline(1.78, color='#64748b', linestyle=':', label='Diffusion-Limited / DLCA ($D_m \\approx 1.78$)')
    ax.set_ylabel(r'Dimension $D_m$', fontsize=11, fontweight='bold')
    ax.set_title(r'Network Mass Fractal Dimension $D_m$', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', fontsize=8.5)

    # Panel B: Correlation Length xi
    ax = axes[0, 1]
    ax.errorbar(t, df_params['xi'], yerr=df_params['std_err_xi'], fmt='s-', color='#2563eb', label=r'Correlation Length $\xi$ [nm]')
    ax.set_ylabel(r'Correlation Cutoff $\xi$ [nm]', fontsize=11, fontweight='bold')
    ax.set_title(r'Percolation Network Correlation Length $\xi$', fontsize=11, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', fontsize=8.5)

    # Panel C: Macroscopic Domain Porod Amplitude
    ax = axes[1, 0]
    ax.errorbar(t, df_params['A_porod'] * 1e6, yerr=df_params['std_err_A_porod'] * 1e6, fmt='d-', color='#7c3aed', label=r'$A_{\rm porod} \times 10^6$ [cm$^{-1}\cdot$nm$^{-4}$]')
    ax.set_xlabel('Reaction Time $t$ [s]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Porod Amplitude $A_{\rm porod} \times 10^6$', fontsize=11, fontweight='bold')
    ax.set_title('Macroscopic Precipitate Surface Scattering (Porod)', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper left', fontsize=8.5)

    # Panel D: Goodness of Fit
    ax = axes[1, 1]
    ax.plot(t, df_params['chi2_red'], 'o-', color='#ea580c', label=r'Reduced $\chi^2_{\rm red}$')
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('Reaction Time $t$ [s]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Reduced $\chi^2$', fontsize=11, fontweight='bold')
    ax.set_title('Goodness of Fit Quality', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=8.5)

    plt.suptitle("ND-i-004 Long-Distance SAXS Kinetic Trajectories (15.0 mM CaCO$_3$, pH 11.3, Hybrid)\n"
                 "Time-Resolved Mesoscale Spinodal Network & Coarsening Evolution", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved long-distance trajectories: {plot_path}", flush=True)
    print(f"Long-distance time-resolved fitting completed in {time.time() - start_time:.1f} s.", flush=True)
    print("=" * 80, flush=True)


if __name__ == '__main__':
    fit_long_distance_series()
