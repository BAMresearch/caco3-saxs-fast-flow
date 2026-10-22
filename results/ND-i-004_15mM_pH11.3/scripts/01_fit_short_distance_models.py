#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 01: Short Distance Time-Resolved Fitting for ND-i-004.
------------------------------------------------------------
Fits the 50 short-distance SAXS frames (3.0 m camera, q in [0.022, 2.54] nm^-1, t = 0.060 s to 10.877 s).

Scientific Principle for 15.0 mM CaCO3 (pH 11.3, Hybrid Regime):
- Mesoscale Spinodal Fractal Network with Primary Building Blocks (R_prim = 3.5 nm, std_prim = 1.5 nm):
  Scatters across the intermediate-to-high q regime via the primary building block form factor P(q).
- Discrete ACC Nanosphere Growth:
  Accurately tracks nanosphere evolution from nascent nuclei (R ~ 15 nm) to mature size (R ~ 44 nm)
  with narrow physical polydispersity (sigma ~ 7.5 nm, ~17%), capturing sharp form-factor Bessel oscillations.
- Macroscopic Porod Surface Scattering:
  Asymptotic Porod scattering (I ~ A_porod / q^4) for coarsened phase-separated domains (q*R >> 1).
- Backward Evolutionary Propagation (Frame 49 -> Frame 0):
  Seeded from the mature, high-SNR Frame 49 and propagated backward to eliminate unphysical parameter trapping.
- Formal +/- 1 sigma parameter standard errors computed via the covariance matrix.

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
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi004_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'short_distance')


def fit_short_distance_series():
    print("=" * 80)
    print("  SHORT-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-004, 50 FRAMES)")
    print("  Model: Spinodal Network (Primary Building Blocks) + ACC Nanospheres + Porod Surface Tail")
    print("  Strategy: Backward Evolutionary Propagation (Frame 49 -> Frame 0)")
    print("=" * 80)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load Data, Metadata & Baseline Anchor
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    # Load baseline anchor values if available
    d_base = 2.45
    xi_base = 250.0
    r_prim_base = 3.5
    std_prim_base = 1.5
    if os.path.exists(BASELINE_JSON):
        try:
            with open(BASELINE_JSON, 'r') as f:
                bp = json.load(f)
            d_base = bp.get('Dm', bp.get('D', 2.45))
            xi_base = bp.get('xi', 250.0)
            r_prim_base = bp.get('R_prim', 3.5)
            std_prim_base = bp.get('std_prim', 1.5)
        except Exception:
            pass

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Primary building block core: R_prim = {r_prim_base:.2f} nm, std_prim = {std_prim_base:.2f} nm")

    # 2. Sequential Backward Fitting Loop (49 down to 0)
    fitted_records = []
    fitted_curves = {'q': q}
    prev_params = None

    for pos in reversed(range(n_frames)):
        row = df_meta.iloc[pos]
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

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)

        # Primary Building Blocks of the Mass Fractal Network
        p.add('R_prim', value=r_prim_base, vary=False)
        p.add('std_prim', value=std_prim_base, vary=False)
        p.add('Dm', value=d_base, vary=False)
        p.add('xi', value=xi_base, vary=False)

        if prev_params is None:
            # Seed at Frame 49: mature nanosphere population with distinct Bessel oscillations
            p.add('R_acc', value=44.25, min=38.0, max=50.0, vary=True)
            p.add('std_acc', value=7.17, min=4.0, max=11.0, vary=True)
            p.add('phi_acc', value=2.7e-5, min=1e-7, max=0.001, vary=True)
            p.add('A_porod', value=7.8e-6, min=1e-9, max=0.01, vary=True)
            p.add('phi_net', value=2.0e-5, min=1e-7, max=0.005, vary=True)
            p.add('bkg', value=1.5e-4, min=1e-6, max=0.01, vary=True)
        else:
            # Backward inheritance: subsequent frame provides optimal starting condition
            r_val = prev_params['R_acc'].value
            std_val = prev_params['std_acc'].value
            phi_val = prev_params['phi_acc'].value
            por_val = prev_params['A_porod'].value
            net_val = prev_params['phi_net'].value
            bkg_val = prev_params['bkg'].value

            p.add('R_acc', value=r_val,
                  min=max(10.0, r_val - 4.0),
                  max=min(55.0, r_val + 1.0), vary=True)
            p.add('std_acc', value=std_val,
                  min=max(4.0, std_val - 1.5),
                  max=min(12.0, std_val + 1.5), vary=True)
            p.add('phi_acc', value=phi_val,
                  min=1e-7, max=max(1e-4, phi_val * 2.0), vary=True)
            p.add('A_porod', value=por_val,
                  min=1e-9, max=max(1e-4, por_val * 2.5), vary=True)
            p.add('phi_net', value=net_val,
                  min=1e-7, max=0.005, vary=True)
            p.add('bkg', value=bkg_val,
                  min=1e-6, max=0.01, vary=True)

        def model_sd(pars, q_eval):
            I_net = calc_primary_fractal_intensity(q_eval, pars['R_prim'].value, pars['std_prim'].value,
                                                   pars['phi_net'].value, pars['Dm'].value, pars['xi'].value, integration_steps=60)
            I_acc = calc_large_sphere_intensity(q_eval, pars['R_acc'].value, pars['std_acc'].value,
                                                 pars['phi_acc'].value, integration_steps=60)
            I_porod = calc_porod_intensity(q_eval, pars['A_porod'].value, 4.0)
            return I_net + I_acc + I_porod + pars['bkg'].value

        def res_sd(pars):
            c = model_sd(pars, q_fit)
            residuals = (y_fit - c) / dy_fit
            if prev_params is not None:
                # Gentle backward continuity regularization
                dR = (pars['R_acc'].value - prev_params['R_acc'].value) / 1.5
                dPhi = (pars['phi_acc'].value - prev_params['phi_acc'].value) / (prev_params['phi_acc'].value * 0.35 + 1e-7)
                dPorod = (pars['A_porod'].value - prev_params['A_porod'].value) / (prev_params['A_porod'].value * 0.40 + 1e-8)
                reg = np.array([dR, dPhi, dPorod])
                return np.concatenate([residuals, reg])
            return residuals

        out_fit = minimize(res_sd, p, method='least_squares', max_nfev=200)

        # Compute exact covariance matrix
        p_cov = out_fit.params.copy()
        for k in ['R_acc', 'std_acc', 'phi_acc', 'A_porod', 'phi_net', 'bkg']:
            p_cov[k].vary = True
        try:
            out_cov = minimize(res_sd, p_cov, method='leastsq', calc_covar=True)
            final_p = out_cov.params
        except Exception:
            final_p = out_fit.params

        # Calculate full curve across full q vector
        fit_full = model_sd(final_p, q)
        fit_fit = model_sd(final_p, q_fit)
        chi2_red = float(np.sum(((y_fit - fit_fit) / dy_fit) ** 2) / max(1, len(y_fit) - 6))

        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'filename': f_file,
            'chi2_red': chi2_red,
            'R_acc': final_p['R_acc'].value,
            'std_err_R_acc': float(final_p['R_acc'].stderr) if final_p['R_acc'].stderr is not None else 0.0,
            'std_acc': final_p['std_acc'].value,
            'std_err_std_acc': float(final_p['std_acc'].stderr) if final_p['std_acc'].stderr is not None else 0.0,
            'phi_acc': final_p['phi_acc'].value,
            'std_err_phi_acc': float(final_p['phi_acc'].stderr) if final_p['phi_acc'].stderr is not None else 0.0,
            'A_porod': final_p['A_porod'].value,
            'std_err_A_porod': float(final_p['A_porod'].stderr) if final_p['A_porod'].stderr is not None else 0.0,
            'phi_net': final_p['phi_net'].value,
            'std_err_phi_net': float(final_p['phi_net'].stderr) if final_p['phi_net'].stderr is not None else 0.0,
            'R_prim': final_p['R_prim'].value,
            'std_prim': final_p['std_prim'].value,
            'Dm': final_p['Dm'].value,
            'xi': final_p['xi'].value,
            'bkg': final_p['bkg'].value,
            'std_err_bkg': float(final_p['bkg'].stderr) if final_p['bkg'].stderr is not None else 0.0,
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_full
        prev_params = final_p

        if f_idx % 10 == 0 or f_idx == 49 or f_idx == 0:
            print(f"  Frame {f_idx:02d} (t = {t_sec:6.3f} s): chi2 = {chi2_red:.3f} | R_acc = {final_p['R_acc'].value:.2f} ± {rec['std_err_R_acc']:.2f} nm, std = {final_p['std_acc'].value:.2f} nm, phi_net = {final_p['phi_net'].value:.2e}")

    # 3. Sort chronologically from frame 0 to 49
    fitted_records.sort(key=lambda x: x['frame_index'])
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)
    print(f"\nSaved short-distance parameters: {params_csv}")

    # Reorder columns of fitted_curves so 'q' is first, followed by '0', '1', ..., '49'
    curve_cols = ['q'] + [str(i) for i in range(n_frames)]
    df_curves = pd.DataFrame(fitted_curves)[curve_cols]
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)
    print(f"Saved short-distance fitted curves: {curves_csv}")

    # 4. Generate Trajectory Figure
    print("Generating short-distance trajectory plots...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    t = df_params['time_seconds'].values

    # Panel A: R_acc and std_acc
    ax = axes[0, 0]
    ax.errorbar(t, df_params['R_acc'], yerr=df_params['std_err_R_acc'], fmt='o-', color='#0284c7', label=r'Mean ACC Radius $\langle R_{\rm ACC} \rangle$')
    ax.fill_between(t, df_params['R_acc'] - df_params['std_acc'], df_params['R_acc'] + df_params['std_acc'], color='#38bdf8', alpha=0.25, label=r'Polydispersity $\pm \sigma_{\rm ACC}$')
    ax.set_ylabel(r'Radius $R$ [nm]', fontsize=11, fontweight='bold')
    ax.set_title('ACC Nanosphere Size & Polydispersity', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=9)

    # Panel B: Volume Fractions
    ax = axes[0, 1]
    ax.errorbar(t, df_params['phi_acc'] * 1e4, yerr=df_params['std_err_phi_acc'] * 1e4, fmt='s-', color='#16a34a', label=r'ACC Volume Fraction $\phi_{\rm ACC} \times 10^4$')
    ax_twin = ax.twinx()
    ax_twin.plot(t, df_params['phi_net'] * 1e4, 'd--', color='#9333ea', label=r'Fractal Network $\phi_{\rm net} \times 10^4$')
    ax.set_ylabel(r'$\phi_{\rm ACC} \times 10^4$', fontsize=11, fontweight='bold', color='#16a34a')
    ax_twin.set_ylabel(r'$\phi_{\rm net} \times 10^4$', fontsize=11, fontweight='bold', color='#9333ea')
    ax.set_title('Volume Fractions (ACC & Spinodal Network)', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax_twin.get_legend_handles_labels()
    ax.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', fontsize=8.5)

    # Panel C: Porod Surface Kinetics
    ax = axes[1, 0]
    ax.plot(t, df_params['A_porod'] * 1e6, 'd-', color='#d97706', label=r'Porod Amplitude $A_{\rm porod} \times 10^6$ [cm$^{-1}$ nm$^4$]')
    ax.set_xlabel('Reaction Time $t$ [s]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$A_{\rm porod} \times 10^6$', fontsize=11, fontweight='bold')
    ax.set_title('Macroscopic Porod Surface Kinetics', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=9)

    # Panel D: Goodness of Fit chi2_red
    ax = axes[1, 1]
    ax.plot(t, df_params['chi2_red'], 'o-', color='#ea580c', label=r'Reduced $\chi^2_{\rm red}$')
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('Reaction Time $t$ [s]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Reduced $\chi^2$', fontsize=11, fontweight='bold')
    ax.set_title('Goodness of Fit Quality', fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=9)

    plt.suptitle("ND-i-004 Short-Distance SAXS Kinetic Trajectories (15.0 mM CaCO$_3$, pH 11.3, Hybrid)\n"
                 "Spinodal Fractal Network (Building Blocks) + ACC Nanosphere Growth & Oscillations", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved short-distance trajectories: {plot_path}")
    print(f"Short-distance time-resolved fitting completed in {time.time() - start_time:.1f} s.")
    print("=" * 80)


if __name__ == '__main__':
    fit_short_distance_series()
