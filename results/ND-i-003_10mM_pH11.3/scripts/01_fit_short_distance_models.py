#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 01: Short Distance Background Modelling for ND-i-003.
----------------------------------------------------------
Fits the 50 short-distance SAXS frames (3.0 m camera, q in [0.022, 2.54] nm^-1, t = 0.060 to 10.877 s).

Scientific Principle for Low Supersaturation (10 mM, pH 11.3):
- No spinodal decomposition or mass fractal gelation occurs.
- Primary building blocks (r0 ~ 8 nm) and high-q concentration fluctuations (xi ~ 2.8 nm)
  observed at 40 mM are completely absent in this low-supersaturation regime.
- The subtracted scattering profile across all 50 frames is statistically identical
  to the flat solvent (water) background (mean I ~ 1e-4 cm^-1, chi^2_red ~ 0.3 - 0.5).

Outputs:
  output/short_distance/parameters.csv
  output/short_distance/fitted_curves.csv
  output/short_distance/trajectories.png
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
from lmfit import minimize, Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import calc_short_distance_model
from models.fitting_engine import create_short_distance_parameters

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_short_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'short_distance')


def fit_short_distance_series():
    print("=" * 80)
    print("  SHORT-DISTANCE SAXS FITTING: FLAT SOLVENT BASELINE VERIFICATION (ND-i-003)")
    print("  Model: Flat Solvent Background (Verifying Absence of Primary Nodes/Fractal Mesh)")
    print("=" * 80)
    start_t = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    fitted_records = []
    fitted_curves = {'q': q}

    for i in range(n_frames):
        row = df_meta.iloc[i]
        f_idx = int(row['frame_index'])
        t_sec = float(row['time_seconds'])
        f_file = str(row['short_file'])

        y = df_data[str(f_idx)].values
        dy = np.maximum(df_data['error_' + str(f_idx)].values, 1e-4)

        mask = (q >= 0.025) & (q <= 2.50)
        q_fit = q[mask]
        y_fit = y[mask]
        dy_fit = dy[mask]

        params = create_short_distance_parameters(bkg_init=float(np.median(y_fit)))

        def residual(p, q_pts, y_pts, dy_pts):
            ymod = calc_short_distance_model(p, q_pts)
            return (y_pts - ymod) / dy_pts

        out = minimize(residual, params, args=(q_fit, y_fit, dy_fit), method='least_squares')

        bkg_val = float(out.params['bkg'].value)
        bkg_err = float(out.params['bkg'].stderr or 0.0)

        # Calculate high-q metrics
        mean_I = float(np.mean(y_fit))
        std_I = float(np.std(y_fit))
        pos_fraction = float((y_fit > 0).mean())

        record = {
            'frame': f_idx,
            'filename': f_file,
            'time_seconds': t_sec,
            'bkg': bkg_val,
            'std_err_bkg': bkg_err,
            'mean_intensity': mean_I,
            'std_intensity': std_I,
            'pos_fraction': pos_fraction,
            'chisqr': float(out.chisqr),
            'redchi': float(out.redchi)
        }
        fitted_records.append(record)

        # Store curve
        fitted_curves[str(f_idx)] = calc_short_distance_model(out.params, q)

        if f_idx % 10 == 0 or f_idx == n_frames - 1:
            print(f"  Frame {f_idx:2d} (t = {t_sec:6.3f}s): bkg = {bkg_val:10.3e} +/- {bkg_err:9.2e}, "
                  f"mean_I = {mean_I:10.3e}, redchi = {out.redchi:.3f}, pos_frac = {pos_fraction:.2f}")

    df_params = pd.DataFrame(fitted_records)
    df_curves = pd.DataFrame(fitted_curves)

    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    df_params.to_csv(params_csv, index=False)
    df_curves.to_csv(curves_csv, index=False)

    # Plot short-distance summary trajectories
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    fig.suptitle('ND-i-003: Short Distance SAXS (3.0 m Camera, 10 mM CaCO3, pH 11.3)\nVerification of Flat Solvent Background & Absence of Fractal Mesh', fontsize=12, fontweight='bold')

    t = df_params['time_seconds'].values

    # 1. Background intensity
    ax = axes[0, 0]
    ax.errorbar(t, df_params['bkg'], yerr=df_params['std_err_bkg'], fmt='o-', color='#1f77b4', ecolor='#aec7e8', capsize=3, markersize=4, label='Fitted Solvent Background')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.6)
    ax.set_ylabel(r'Background $I_{\rm bkg}$ (cm$^{-1}$)', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9)

    # 2. Mean Intensity across scattering curve
    ax = axes[0, 1]
    ax.plot(t, df_params['mean_intensity'], 's-', color='#2ca02c', markersize=4, label='Mean Subtracted Intensity')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.6)
    ax.set_ylabel(r'Mean Subtracted $I(q)$ (cm$^{-1}$)', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9)

    # 3. Fraction of positive points (pure noise = ~0.5)
    ax = axes[1, 0]
    ax.plot(t, df_params['pos_fraction'] * 100.0, '^-', color='#ff7f0e', markersize=4, label='% Positive Data Points')
    ax.axhline(50.0, color='red', linestyle='--', alpha=0.7, label='Pure Gaussian Noise (50%)')
    ax.set_xlabel(r'Reaction Time $t$ (s)', fontsize=11)
    ax.set_ylabel(r'Positive Points (%)', fontsize=11)
    ax.set_ylim(40, 75)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9)

    # 4. Reduced Chi2
    ax = axes[1, 1]
    ax.plot(t, df_params['redchi'], 'd-', color='#d62728', markersize=4, label=r'$\chi^2_{\rm red}$ (Flat Model)')
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5, label='Ideal $\\chi^2_{\\rm red} = 1$')
    ax.set_xlabel(r'Reaction Time $t$ (s)', fontsize=11)
    ax.set_ylabel(r'Reduced $\chi^2$', fontsize=11)
    ax.set_ylim(0, 1.5)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'trajectories.png'), dpi=300)
    plt.close(fig)

    elapsed = time.time() - start_t
    print(f"\nShort-distance fitting completed in {elapsed:.2f} s. Outputs saved to: {OUT_DIR}")
    return df_params


if __name__ == '__main__':
    fit_short_distance_series()
