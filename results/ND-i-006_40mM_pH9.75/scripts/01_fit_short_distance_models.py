#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 01: Short Distance Time-Resolved Fitting for ND-i-006.
------------------------------------------------------------
Fits the 50 short-distance SAXS frames (3.0 m camera, q in [0.022, 2.54] nm^-1, t = 0.060 s to 10.877 s).

Scientific Principle for 40.0 mM CaCO3 (pH 9.75):
- Mesoscale Spinodal Fractal Network with Primary Building Blocks (R_prim, std_prim):
  Scatters across the intermediate-to-high q regime via the primary building block form factor P(q) and S(q).
- Macroscopic Porod Surface Scattering:
  Asymptotic Porod scattering (I ~ A_porod / q^4) for coarsened phase-separated large droplets/spheres (q*R >> 1).
- Incoherent flat background (bkg).
- NO ACC NANOSPHERES (experimentally absent in this condition).
- Backward Evolutionary Propagation (Frame 49 -> Frame 0):
  Seeded from the mature, high-SNR Frame 49 and propagated backward to eliminate unphysical parameter trapping.
- Parameter standard errors (+/- 1 sigma) computed via the covariance matrix.

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
    calc_porod_intensity,
    RHO2_CALIBRATED
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'short_distance')


def fit_short_distance_series():
    print("=" * 80)
    print("  SHORT-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-006, 50 FRAMES)")
    print("  Model: Spinodal Network (Primary Building Blocks) + Porod Surface Tail (No ACC Particles)")
    print("  Strategy: Backward Evolutionary Propagation (Frame 49 -> Frame 0)")
    print("=" * 80)
    start_time = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load Data, Metadata & Baseline Anchor
    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    d_base = 1.80
    xi_base = 30000.0
    r_prim_base = 1.50
    std_prim_base = 3.03
    bkg_base = 2.3e-4
    a_porod_base = 1.0e-6

    if os.path.exists(BASELINE_JSON):
        try:
            with open(BASELINE_JSON, 'r') as f:
                bp = json.load(f)
            d_base = bp.get('Dm', bp.get('D', 1.80))
            xi_base = bp.get('xi', 30000.0)
            r_prim_base = bp.get('R_prim', 1.50)
            std_prim_base = bp.get('std_prim', 3.03)
            bkg_base = bp.get('bkg', 2.3e-4)
            a_porod_base = bp.get('A_porod', 1.0e-6)
        except Exception:
            pass

    print(f"Loaded {n_frames} frames from metadata (t: {df_meta['time_seconds'].min():.4f}s to {df_meta['time_seconds'].max():.4f}s)")
    print(f"Primary building block core: R_prim = {r_prim_base:.2f} nm, std_prim = {std_prim_base:.2f} nm, Dm = {d_base:.3f}")

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

        mask = (q >= 0.022) & (y > 0)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = dy[mask]

        p = Parameters()
        p.add('rho2', value=RHO2_CALIBRATED, vary=False)

        p.add('R_prim', value=r_prim_base, vary=False)
        p.add('std_prim', value=std_prim_base, vary=False)
        p.add('xi', value=xi_base, vary=False)

        if prev_params is None:
            # Seed at Frame 49
            p.add('phi_net', value=2.5e-5, min=1e-7, max=0.005, vary=True)
            p.add('Dm', value=d_base, min=1.1, max=2.85, vary=True)
            p.add('A_porod', value=a_porod_base, min=1e-10, max=0.01, vary=True)
            p.add('bkg', value=bkg_base, min=1e-6, max=0.01, vary=True)
        else:
            phi_val = prev_params['phi_net'].value
            dm_val = prev_params['Dm'].value
            por_val = prev_params['A_porod'].value
            bkg_val = prev_params['bkg'].value

            p.add('phi_net', value=phi_val, min=1e-7, max=0.005, vary=True)
            p.add('Dm', value=dm_val, min=max(1.1, dm_val - 0.25), max=min(2.85, dm_val + 0.25), vary=True)
            p.add('A_porod', value=por_val, min=1e-10, max=max(1e-4, por_val * 2.5), vary=True)
            p.add('bkg', value=bkg_val, min=1e-6, max=0.01, vary=True)

        def model_sd(pars, q_eval):
            I_net = calc_primary_fractal_intensity(
                q_eval, pars['R_prim'].value, pars['std_prim'].value,
                pars['phi_net'].value, pars['Dm'].value, pars['xi'].value,
                integration_steps=60
            )
            I_porod = calc_porod_intensity(q_eval, pars['A_porod'].value, 4.0)
            return I_net + I_porod + pars['bkg'].value

        def res_sd(pars):
            c = model_sd(pars, q_fit)
            residuals = (y_fit - c) / dy_fit
            if prev_params is not None:
                dDm = (pars['Dm'].value - prev_params['Dm'].value) / 0.15
                dPhi = (pars['phi_net'].value - prev_params['phi_net'].value) / (prev_params['phi_net'].value * 0.40 + 1e-7)
                reg = np.array([dDm, dPhi])
                return np.concatenate([residuals, reg])
            return residuals

        out_fit = minimize(res_sd, p, method='least_squares', max_nfev=200)

        p_cov = out_fit.params.copy()
        for k in ['phi_net', 'Dm', 'A_porod', 'bkg']:
            p_cov[k].vary = True
        try:
            out_cov = minimize(res_sd, p_cov, method='leastsq', calc_covar=True)
            final_p = out_cov.params
        except Exception:
            final_p = out_fit.params

        fit_full = model_sd(final_p, q)
        fit_fit = model_sd(final_p, q_fit)
        chi2_red = float(np.sum(((y_fit - fit_fit) / dy_fit) ** 2) / max(1, len(y_fit) - 4))

        rec = {
            'frame_index': f_idx,
            'time_seconds': t_sec,
            'filename': f_file,
            'chi2_red': chi2_red,
            'phi_net': final_p['phi_net'].value,
            'std_err_phi_net': float(final_p['phi_net'].stderr) if final_p['phi_net'].stderr is not None else 0.0,
            'Dm': final_p['Dm'].value,
            'std_err_Dm': float(final_p['Dm'].stderr) if final_p['Dm'].stderr is not None else 0.0,
            'R_prim': final_p['R_prim'].value,
            'std_prim': final_p['std_prim'].value,
            'xi': final_p['xi'].value,
            'A_porod': final_p['A_porod'].value,
            'std_err_A_porod': float(final_p['A_porod'].stderr) if final_p['A_porod'].stderr is not None else 0.0,
            'bkg': final_p['bkg'].value,
            'std_err_bkg': float(final_p['bkg'].stderr) if final_p['bkg'].stderr is not None else 0.0,
        }
        fitted_records.append(rec)
        fitted_curves[str(f_idx)] = fit_full
        prev_params = final_p

        if f_idx % 10 == 0 or f_idx == n_frames - 1 or f_idx == 0:
            print(f"  Frame {f_idx:02d} (t = {t_sec:6.3f} s): chi2 = {chi2_red:.3f} | Dm = {final_p['Dm'].value:.3f} ± {rec['std_err_Dm']:.3f}, phi_net = {final_p['phi_net'].value:.2e}")

    # 3. Save parameters and curves sorted chronologically
    fitted_records.sort(key=lambda x: x['frame_index'])
    df_params = pd.DataFrame(fitted_records)
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    df_params.to_csv(params_csv, index=False)

    df_curves = pd.DataFrame(fitted_curves)
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_curves.to_csv(curves_csv, index=False)

    # 4. Generate Trajectory Overview Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("ND-i-006 Short-Distance Kinetic Evolution: Mass Fractal Network (40.0 mM CaCO3, pH 9.75)", fontsize=13, fontweight='bold')

    t_arr = df_params['time_seconds']

    # Panel 1: Dm
    axes[0, 0].plot(t_arr, df_params['Dm'], 'o-', color='#d97706', label=r'$D_m$')
    axes[0, 0].fill_between(t_arr, df_params['Dm'] - df_params['std_err_Dm'], df_params['Dm'] + df_params['std_err_Dm'], color='#fbbf24', alpha=0.3, label=r'$\pm 1\sigma$')
    axes[0, 0].set_ylabel(r"Mass Fractal Dimension $D_m$", fontweight='bold')
    axes[0, 0].set_title(r"Fractal Dimension Evolution $D_m(t)$", fontweight='bold')
    axes[0, 0].grid(True, linestyle='--', alpha=0.5)
    axes[0, 0].legend()

    # Panel 2: Volume Fraction
    axes[0, 1].plot(t_arr, df_params['phi_net'] * 1e4, '^-', color='#16a34a', label=r'$\phi_{\mathrm{net}} \times 10^4$')
    axes[0, 1].fill_between(t_arr, (df_params['phi_net'] - df_params['std_err_phi_net']) * 1e4, (df_params['phi_net'] + df_params['std_err_phi_net']) * 1e4, color='#16a34a', alpha=0.25)
    axes[0, 1].set_ylabel(r"$\phi_{\mathrm{net}} \times 10^4$", fontweight='bold')
    axes[0, 1].set_title("Network Volume Fraction Kinetics", fontweight='bold')
    axes[0, 1].grid(True, linestyle='--', alpha=0.5)
    axes[0, 1].legend()

    # Panel 3: Porod Tail
    axes[1, 0].plot(t_arr, df_params['A_porod'] * 1e6, 'd-', color='#d97706', label=r'$A_{\mathrm{porod}} \times 10^6$')
    axes[1, 0].set_ylabel(r"$A_{\mathrm{porod}} \times 10^6$", fontweight='bold')
    axes[1, 0].set_xlabel("Time [s]", fontweight='bold')
    axes[1, 0].set_title("Macroscopic Porod Surface Tail", fontweight='bold')
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
    print(f"Short distance fitting completed in {elapsed:.1f} s. Saved to {OUT_DIR}")


if __name__ == '__main__':
    fit_short_distance_series()
