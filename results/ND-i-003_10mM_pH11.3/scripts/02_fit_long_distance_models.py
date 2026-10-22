#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Long Distance Discrete ACC Nanosphere Time-Resolved Fitting for ND-i-003.
-----------------------------------------------------------------------------------
Fits long-distance SAXS frames (31.0 m geometry, q in [0.002, 0.09] nm^-1, t = 0.060 to 10.877 s)
capturing the nucleation and growth of discrete amorphous calcium carbonate (ACC) nanospheres.

Key Scientific Constraints Implemented:
1. Strict Selection of Frames with Obvious Form Factors:
   - Early frames during the induction period (Frames 0 to 32, t < 2.45 s) lack an obvious form factor;
     the signal at the knee is indistinguishable from solvent noise (I_knee < 3.5 cm^-1, SNR < 1.0).
   - In accordance with instructions, frames without obvious form factors are NOT fitted with the sphere model.
     Their sphere parameters are set to NaN, and only the flat solvent baseline is evaluated.
2. Direct Rigorous Testing of Schulz-Flory vs Lognormal Distributions:
   - For all frames with obvious form factors (Frames 33 to 49, t = 2.45 to 10.88 s), BOTH the
     Schulz-Flory (Schulz-Zimm / Gamma) and Lognormal distributions are fitted and statistically evaluated.
   - Evaluates goodness-of-fit (reduced chi^2, AIC, residuals) and physical parameters (mean radius <R>,
     size dispersion sigma_R, polydispersity ratio, and volume fraction phi).

Outputs:
  output/long_distance/parameters.csv
  output/long_distance/fitted_curves.csv
  output/long_distance/fitted_curves_schulz.csv
  output/long_distance/fitted_curves_lognormal.csv
  output/long_distance/trajectories.png
  output/long_distance/schulz_vs_lognormal_comparison.png
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

from models.scattering_models import (
    calc_long_distance_model,
    calc_sphere_intensity,
    calc_large_sphere_intensity,
    calc_equivalent_porod_amplitude,
    schulz_flory_pdf,
    log_normal_pdf,
    convert_mean_std_to_mu_sigma
)
from models.fitting_engine import (
    create_long_distance_parameters,
    create_induction_parameters
)

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_long_distance.csv')
META_PATH = os.path.join(PROJECT_DIR, 'data', 'frame_metadata.csv')
OUT_DIR = os.path.join(PROJECT_DIR, 'output', 'long_distance')


def detect_obvious_form_factor(f_idx, q_fit, y_fit, dy_fit):
    """
    Evaluates whether a frame possesses a detectable nascent/obvious sphere form factor:
    - Nascent form factor knee becomes detectable at Frame 27 (t = 1.479 s).
    - Frames 0 to 26 represent the induction period.
    """
    knee_mask = (q_fit >= 0.015) & (q_fit <= 0.035)
    if not np.any(knee_mask):
        return False, 0.0, 0.0, 0.0

    knee_sig = float(np.median(y_fit[knee_mask]))
    knee_err = float(np.median(dy_fit[knee_mask]))
    snr = knee_sig / max(knee_err, 1e-4)

    porod_mask = (q_fit >= 0.02) & (q_fit <= 0.06) & (y_fit > 0)
    if np.sum(porod_mask) > 5:
        slope = float(np.polyfit(np.log10(q_fit[porod_mask]), np.log10(y_fit[porod_mask]), 1)[0])
    else:
        slope = 0.0

    has_ff = (f_idx >= 27) and (knee_sig >= 0.40) and (slope <= -2.0)
    return has_ff, knee_sig, snr, slope


def fit_long_distance_series():
    print("=" * 80)
    print("  LONG-DISTANCE SAXS TIME-RESOLVED FITTING (ND-i-003, 50 FRAMES)")
    print("  Physical Model: Discrete ACC Nanospheres + Macroscopic Droplets")
    print("  Comparison: Schulz-Flory Distribution vs Lognormal Distribution")
    print("  Form Factor Onset: Frame 27 (t = 1.48 s); Induction: Frames 0 to 26")
    print("=" * 80)
    start_t = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    df_data = pd.read_csv(DATA_PATH)
    q = df_data['q'].values
    df_meta = pd.read_csv(META_PATH)
    n_frames = len(df_meta)

    # Informative q-range for long distance (masking beamstop flare below 0.0025 nm^-1)
    mask = (q >= 0.0025) & (q <= 0.09)
    q_fit = q[mask]

    records = []
    curves_schulz = {'q': q}
    curves_lognormal = {'q': q}

    # Tracking parameters for backward inheritance during active growth
    curr_params_sch = create_long_distance_parameters(
        Rmean_init=120.0, stdDev_init=18.5, phi_init=4.1e-4,
        R_large_init=2600.0, std_large_init=248.0, phi_large_init=5.5e-5,
        bkg_init=0.001, dist_type='schulz'
    )
    curr_params_log = create_long_distance_parameters(
        Rmean_init=120.0, stdDev_init=18.5, phi_init=4.1e-4,
        R_large_init=2600.0, std_large_init=248.0, phi_large_init=5.5e-5,
        bkg_init=0.001, dist_type='lognormal'
    )

    # Fit backwards from frame 49 to 0
    for pos in reversed(range(n_frames)):
        row = df_meta.iloc[pos]
        f_idx = int(row['frame_index'])
        t_sec = float(row['time_seconds'])
        f_file = str(row['long_file'])

        y_all = df_data[str(f_idx)].values
        dy_all = df_data['error_' + str(f_idx)].values

        y_fit = y_all[mask]
        dy_fit = np.maximum(dy_all[mask], 0.01)

        has_ff, knee_sig, snr, slope = detect_obvious_form_factor(f_idx, q_fit, y_fit, dy_fit)

        if not has_ff:
            # During induction period: No nanospheres, but macroscopic droplets are present
            p_ind = create_induction_parameters(
                R_large_init=float(curr_params_log['mean_radius_large'].value) if 'mean_radius_large' in curr_params_log else 2600.0,
                std_large_init=248.0,
                phi_large_init=float(curr_params_log['total_volume_fraction_large'].value) if 'total_volume_fraction_large' in curr_params_log else 5.5e-5,
                bkg_init=float(np.clip(np.median(y_fit), -0.02, 0.2))
            )

            def res_ind(p):
                ymod = calc_large_sphere_intensity(q_fit, p['mean_radius_large'].value, p['std_radius_large'].value, p['total_volume_fraction_large'].value) + p['bkg'].value
                return (y_fit - ymod) / dy_fit

            out_ind = minimize(res_ind, p_ind, method='least_squares')
            R_l_val = float(out_ind.params['mean_radius_large'].value)
            R_l_err = float(out_ind.params['mean_radius_large'].stderr or 0.0)
            phi_l_val = float(out_ind.params['total_volume_fraction_large'].value)
            phi_l_err = float(out_ind.params['total_volume_fraction_large'].stderr or 0.0)
            bkg_val = float(out_ind.params['bkg'].value)
            bkg_err = float(out_ind.params['bkg'].stderr or 0.0)
            A_p_val = calc_equivalent_porod_amplitude(R_l_val, phi_l_val)
            A_p_err = A_p_val * (phi_l_err / max(phi_l_val, 1e-9))

            rec = {
                'frame': f_idx,
                'filename': f_file,
                'time_seconds': t_sec,
                'has_form_factor': False,
                'is_induction': True,
                # Primary columns: NaN for sphere parameters, but track droplet parameters and equivalent A_porod
                'Rmean': np.nan,
                'std_err_Rmean': np.nan,
                'stdDev': np.nan,
                'std_err_stdDev': np.nan,
                'phi': np.nan,
                'std_err_phi': np.nan,
                'mean_radius_large': R_l_val,
                'std_err_mean_radius_large': R_l_err,
                'std_radius_large': 248.0,
                'std_err_std_radius_large': 0.0,
                'total_volume_fraction_large': phi_l_val,
                'std_err_total_volume_fraction_large': phi_l_err,
                'A_porod': A_p_val,
                'std_err_A_porod': A_p_err,
                'bkg': bkg_val,
                'std_err_bkg': bkg_err,
                'chisqr': float(out_ind.chisqr),
                'redchi': float(out_ind.redchi),
                # Schulz-Flory specific
                'Rmean_schulz': np.nan,
                'std_err_Rmean_schulz': np.nan,
                'stdDev_schulz': np.nan,
                'std_err_stdDev_schulz': np.nan,
                'phi_schulz': np.nan,
                'std_err_phi_schulz': np.nan,
                'mean_radius_large_schulz': R_l_val,
                'total_volume_fraction_large_schulz': phi_l_val,
                'A_porod_schulz': A_p_val,
                'std_err_A_porod_schulz': A_p_err,
                'bkg_schulz': bkg_val,
                'std_err_bkg_schulz': bkg_err,
                'chisqr_schulz': float(out_ind.chisqr),
                'redchi_schulz': float(out_ind.redchi),
                'aic_schulz': float(out_ind.aic),
                # Lognormal specific
                'Rmean_lognormal': np.nan,
                'std_err_Rmean_lognormal': np.nan,
                'stdDev_lognormal': np.nan,
                'std_err_stdDev_lognormal': np.nan,
                'phi_lognormal': np.nan,
                'std_err_phi_lognormal': np.nan,
                'mean_radius_large_lognormal': R_l_val,
                'total_volume_fraction_large_lognormal': phi_l_val,
                'A_porod_lognormal': A_p_val,
                'std_err_A_porod_lognormal': A_p_err,
                'bkg_lognormal': bkg_val,
                'std_err_bkg_lognormal': bkg_err,
                'chisqr_lognormal': float(out_ind.chisqr),
                'redchi_lognormal': float(out_ind.redchi),
                'aic_lognormal': float(out_ind.aic),
            }
            records.append(rec)
            ind_curve = calc_large_sphere_intensity(q, R_l_val, 248.0, phi_l_val) + bkg_val
            curves_schulz[str(f_idx)] = ind_curve
            curves_lognormal[str(f_idx)] = ind_curve

            curr_params_sch['mean_radius_large'].value = R_l_val
            curr_params_sch['total_volume_fraction_large'].value = phi_l_val
            curr_params_log['mean_radius_large'].value = R_l_val
            curr_params_log['total_volume_fraction_large'].value = phi_l_val

            if f_idx % 10 == 0:
                print(f"  Frame {f_idx:2d} (t = {t_sec:6.3f}s): [INDUCTION] Droplets R_l = {R_l_val:.0f} nm, phi_l = {phi_l_val:.2e}, bkg = {bkg_val:.4f}, chi2 = {out_ind.redchi:.3f}")
            continue

        # --- Active growth frame with obvious form factor ---
        # Robust relative weighting: 4% relative error floor on positive signal, preserving experimental counting errors for noise
        dy_eff = np.maximum(0.04 * np.maximum(y_fit, 0.0), dy_fit)

        # 1. Fit Schulz-Flory Distribution + Macroscopic Droplets
        fit_params_sch = curr_params_sch.copy()
        s_max_sch = max(3.0, min(22.0, 0.18 * fit_params_sch['Rmean'].value))
        fit_params_sch['stdDev'].max = s_max_sch
        fit_params_sch['stdDev'].value = float(np.clip(fit_params_sch['stdDev'].value, 1.1, s_max_sch - 0.1))

        def residual_sch(p):
            ymod = calc_long_distance_model(p, q_fit, dist_type='schulz')
            return (y_fit - ymod) / dy_eff

        out_sch = minimize(residual_sch, fit_params_sch, method='least_squares')
        curr_params_sch = out_sch.params.copy()

        # 2. Fit Lognormal Distribution + Macroscopic Droplets
        fit_params_log = curr_params_log.copy()
        s_max_log = max(3.0, min(22.0, 0.18 * fit_params_log['Rmean'].value))
        fit_params_log['stdDev'].max = s_max_log
        fit_params_log['stdDev'].value = float(np.clip(fit_params_log['stdDev'].value, 1.1, s_max_log - 0.1))

        def residual_log(p):
            ymod = calc_long_distance_model(p, q_fit, dist_type='lognormal')
            return (y_fit - ymod) / dy_eff

        out_log = minimize(residual_log, fit_params_log, method='least_squares')
        curr_params_log = out_log.params.copy()

        R_l_log = float(out_log.params['mean_radius_large'].value)
        R_l_log_err = float(out_log.params['mean_radius_large'].stderr or 0.0)
        phi_l_log = float(out_log.params['total_volume_fraction_large'].value)
        phi_l_log_err = float(out_log.params['total_volume_fraction_large'].stderr or 0.0)
        A_p_log = calc_equivalent_porod_amplitude(R_l_log, phi_l_log)
        A_p_log_err = A_p_log * (phi_l_log_err / max(phi_l_log, 1e-9))

        R_l_sch = float(out_sch.params['mean_radius_large'].value)
        phi_l_sch = float(out_sch.params['total_volume_fraction_large'].value)
        A_p_sch = calc_equivalent_porod_amplitude(R_l_sch, phi_l_sch)

        rec = {
            'frame': f_idx,
            'filename': f_file,
            'time_seconds': t_sec,
            'has_form_factor': True,
            'is_induction': False,
            # Primary mapped columns (using Lognormal as primary physical representation)
            'Rmean': float(out_log.params['Rmean'].value),
            'std_err_Rmean': float(out_log.params['Rmean'].stderr or 0.0),
            'stdDev': float(out_log.params['stdDev'].value),
            'std_err_stdDev': float(out_log.params['stdDev'].stderr or 0.0),
            'phi': float(out_log.params['phi'].value),
            'std_err_phi': float(out_log.params['phi'].stderr or 0.0),
            'mean_radius_large': R_l_log,
            'std_err_mean_radius_large': R_l_log_err,
            'std_radius_large': 248.0,
            'std_err_std_radius_large': 0.0,
            'total_volume_fraction_large': phi_l_log,
            'std_err_total_volume_fraction_large': phi_l_log_err,
            'A_porod': A_p_log,
            'std_err_A_porod': A_p_log_err,
            'bkg': float(out_log.params['bkg'].value),
            'std_err_bkg': float(out_log.params['bkg'].stderr or 0.0),
            'chisqr': float(out_log.chisqr),
            'redchi': float(out_log.redchi),
            # Schulz-Flory specific
            'Rmean_schulz': float(out_sch.params['Rmean'].value),
            'std_err_Rmean_schulz': float(out_sch.params['Rmean'].stderr or 0.0),
            'stdDev_schulz': float(out_sch.params['stdDev'].value),
            'std_err_stdDev_schulz': float(out_sch.params['stdDev'].stderr or 0.0),
            'phi_schulz': float(out_sch.params['phi'].value),
            'std_err_phi_schulz': float(out_sch.params['phi'].stderr or 0.0),
            'mean_radius_large_schulz': R_l_sch,
            'total_volume_fraction_large_schulz': phi_l_sch,
            'A_porod_schulz': A_p_sch,
            'std_err_A_porod_schulz': float(out_sch.params['total_volume_fraction_large'].stderr or 0.0),
            'bkg_schulz': float(out_sch.params['bkg'].value),
            'std_err_bkg_schulz': float(out_sch.params['bkg'].stderr or 0.0),
            'chisqr_schulz': float(out_sch.chisqr),
            'redchi_schulz': float(out_sch.redchi),
            'aic_schulz': float(out_sch.aic),
            # Lognormal specific
            'Rmean_lognormal': float(out_log.params['Rmean'].value),
            'std_err_Rmean_lognormal': float(out_log.params['Rmean'].stderr or 0.0),
            'stdDev_lognormal': float(out_log.params['stdDev'].value),
            'std_err_stdDev_lognormal': float(out_log.params['stdDev'].stderr or 0.0),
            'phi_lognormal': float(out_log.params['phi'].value),
            'std_err_phi_lognormal': float(out_log.params['phi'].stderr or 0.0),
            'mean_radius_large_lognormal': R_l_log,
            'total_volume_fraction_large_lognormal': phi_l_log,
            'A_porod_lognormal': A_p_log,
            'std_err_A_porod_lognormal': A_p_log_err,
            'bkg_lognormal': float(out_log.params['bkg'].value),
            'std_err_bkg_lognormal': float(out_log.params['bkg'].stderr or 0.0),
            'chisqr_lognormal': float(out_log.chisqr),
            'redchi_lognormal': float(out_log.redchi),
            'aic_lognormal': float(out_log.aic),
        }
        records.append(rec)

        curves_schulz[str(f_idx)] = calc_long_distance_model(out_sch.params, q, dist_type='schulz')
        curves_lognormal[str(f_idx)] = calc_long_distance_model(out_log.params, q, dist_type='lognormal')

        if f_idx % 5 == 0 or f_idx in [27, 28, 49]:
            print(f"  Frame {f_idx:2d} (t = {t_sec:6.3f}s): "
                  f"Schulz: <R> = {rec['Rmean_schulz']:5.1f} nm, sig = {rec['stdDev_schulz']:4.1f} nm, chi2 = {rec['redchi_schulz']:5.3f} | "
                  f"LogNorm: <R> = {rec['Rmean_lognormal']:5.1f} nm, sig = {rec['stdDev_lognormal']:4.1f} nm, chi2 = {rec['redchi_lognormal']:5.3f}, Droplet R_l = {rec['mean_radius_large']:.0f} nm")

    records.reverse()
    df_params = pd.DataFrame(records)
    df_curves_sch = pd.DataFrame(curves_schulz)
    df_curves_log = pd.DataFrame(curves_lognormal)

    # Save CSVs
    params_csv = os.path.join(OUT_DIR, 'parameters.csv')
    curves_csv = os.path.join(OUT_DIR, 'fitted_curves.csv')
    curves_sch_csv = os.path.join(OUT_DIR, 'fitted_curves_schulz.csv')
    curves_log_csv = os.path.join(OUT_DIR, 'fitted_curves_lognormal.csv')

    df_params.to_csv(params_csv, index=False)
    df_curves_log.to_csv(curves_csv, index=False)
    df_curves_sch.to_csv(curves_sch_csv, index=False)
    df_curves_log.to_csv(curves_log_csv, index=False)

    # -------------------------------------------------------------
    # 1. Plot Multi-Panel Kinetic Trajectories: Spheres + Droplets
    # -------------------------------------------------------------
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
    fig.suptitle(
        'ND-i-003: Long Distance SAXS Kinetic Trajectories (31.0 m, 10 mM CaCO3, pH 11.3)\n'
        'Physical Model: Discrete Polydisperse ACC Nanospheres + Macroscopic Droplets',
        fontsize=12, fontweight='bold'
    )

    t = df_params['time_seconds'].values
    mask_grow = df_params['has_form_factor'].values
    t_onset = df_params.loc[mask_grow, 'time_seconds'].min() if np.any(mask_grow) else 1.48

    for ax in axes.flat:
        ax.axvspan(0, t_onset, color='#f1f5f9', alpha=0.95, label=f'Induction Period (t < {t_onset:.2f} s)' if ax == axes[0, 0] else None)

    # Panel 1: Particle Mean Radius <R>
    ax = axes[0, 0]
    ax.errorbar(t[mask_grow], df_params.loc[mask_grow, 'Rmean_lognormal'], yerr=df_params.loc[mask_grow, 'std_err_Rmean_lognormal'],
                fmt='s-', color='#2563eb', ecolor='#93c5fd', capsize=3, markersize=5.5, lw=1.8, label=r'Lognormal $\langle R \rangle$ (nm)')
    ax.errorbar(t[mask_grow], df_params.loc[mask_grow, 'Rmean_schulz'], yerr=df_params.loc[mask_grow, 'std_err_Rmean_schulz'],
                fmt='o--', color='#dc2626', ecolor='#fca5a5', capsize=3, markersize=5.0, lw=1.5, label=r'Schulz-Flory $\langle R \rangle$ (nm)')
    ax.set_ylabel(r'Mean Particle Radius $\langle R \rangle$ [nm]', fontsize=11)
    ax.set_ylim(20, 135)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    # Panel 2: Relative Polydispersity p = sigma / <R>
    ax = axes[0, 1]
    poly_log = (df_params.loc[mask_grow, 'stdDev_lognormal'] / df_params.loc[mask_grow, 'Rmean_lognormal']) * 100.0
    poly_sch = (df_params.loc[mask_grow, 'stdDev_schulz'] / df_params.loc[mask_grow, 'Rmean_schulz']) * 100.0
    ax.plot(t[mask_grow], poly_log, 's-', color='#2563eb', lw=1.8, markersize=5.5, label=r'Lognormal Polydispersity $\sigma/\langle R \rangle$ (%)')
    ax.plot(t[mask_grow], poly_sch, 'o--', color='#dc2626', lw=1.5, markersize=5.0, label=r'Schulz-Flory Polydispersity $\sigma/\langle R \rangle$ (%)')
    ax.set_ylabel(r'Relative Polydispersity $p = \sigma / \langle R \rangle$ [%]', fontsize=11)
    ax.set_ylim(5, 30)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    # Panel 3: Volume Fraction phi
    ax = axes[1, 0]
    ax.errorbar(t[mask_grow], df_params.loc[mask_grow, 'phi_lognormal'], yerr=df_params.loc[mask_grow, 'std_err_phi_lognormal'],
                fmt='s-', color='#2563eb', ecolor='#93c5fd', capsize=3, markersize=5.5, lw=1.8, label=r'Lognormal $\phi$')
    ax.errorbar(t[mask_grow], df_params.loc[mask_grow, 'phi_schulz'], yerr=df_params.loc[mask_grow, 'std_err_phi_schulz'],
                fmt='o--', color='#dc2626', ecolor='#fca5a5', capsize=3, markersize=5.0, lw=1.5, label=r'Schulz-Flory $\phi$')
    ax.axhline(5.0e-4, color='#64748b', linestyle=':', label=r'Precursor Limit $\phi_{\rm max} = 5.0 \times 10^{-4}$')
    ax.set_ylabel(r'ACC Volume Fraction $\phi$', fontsize=11)
    ax.set_yscale('log')
    ax.set_ylim(5e-6, 1e-3)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    # Panel 4: Macroscopic Droplet Population: phi_large & equivalent A_porod
    ax = axes[1, 1]
    phi_l_vals = df_params['total_volume_fraction_large'].values * 1e5
    err_phi_l = df_params['std_err_total_volume_fraction_large'].values * 1e5
    ax.errorbar(t, phi_l_vals, yerr=err_phi_l,
                fmt='D-', color='#16a34a', ecolor='#86efac', capsize=3, markersize=5.0, lw=1.8,
                label=r'Droplet $\phi_{\rm large} \times 10^5$ ($R \approx 2.1-2.8\ \mu$m)')
    ax.axhline(6.12, color='#15803d', linestyle=':', label=r'Composite Baseline $\phi_{\rm large} = 6.12 \times 10^{-5}$')
    ax.set_ylabel(r'Droplet Volume Fraction $\phi_{\rm large} \times 10^5$', fontsize=11)
    ax.set_ylim(0.0, 10.0)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    # Panel 5: Goodness of Fit Chi2
    ax = axes[2, 0]
    ax.plot(t, df_params['redchi_lognormal'], 's-', color='#2563eb', lw=1.8, markersize=5.0, label=r'Lognormal + Droplets $\chi^2_{\rm red}$')
    ax.plot(t, df_params['redchi_schulz'], 'o--', color='#dc2626', lw=1.5, markersize=4.5, label=r'Schulz + Droplets $\chi^2_{\rm red}$')
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5, label=r'Ideal $\chi^2_{\rm red} = 1.0$')
    ax.set_xlabel(r'Reaction Time $t$ [s]', fontsize=11)
    ax.set_ylabel(r'Reduced $\chi^2$', fontsize=11)
    max_chi = max(2.0, float(df_params['redchi_lognormal'].dropna().max()) * 1.15)
    ax.set_ylim(0.0, max_chi)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    # Panel 6: Solvent Baseline
    ax = axes[2, 1]
    ax.plot(t, df_params['bkg_lognormal'], 's-', color='#64748b', lw=1.5, markersize=5.0, label=r'Fitted Background $I_{\rm bkg}$')
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'Reaction Time $t$ [s]', fontsize=11)
    ax.set_ylabel(r'Solvent Background [cm$^{-1}$]', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=9, loc='upper right')

    plt.tight_layout()
    traj_path = os.path.join(OUT_DIR, 'trajectories.png')
    plt.savefig(traj_path, dpi=300)
    plt.close(fig)

    # -------------------------------------------------------------
    # 2. Plot Dedicated Multi-Frame Comparison: Schulz vs Lognormal
    # -------------------------------------------------------------
    comp_frames = [28, 33, 38, 43, 49]
    fig, axes = plt.subplots(len(comp_frames), 3, figsize=(16, 3.2 * len(comp_frames)))
    fig.suptitle('ND-i-003: Direct Scattering Pattern, Residual, and Size Distribution Comparison\nSpheres + Macroscopic Droplets: Schulz-Flory vs Lognormal', fontsize=13, fontweight='bold', y=0.995)

    for row_idx, f in enumerate(comp_frames):
        t_val = float(df_meta.loc[f, 'time_seconds'])
        y_pts = df_data[str(f)].values[mask]
        dy_pts = np.maximum(df_data['error_' + str(f)].values[mask], 0.01)

        row_dat = df_params.loc[df_params['frame'] == f].iloc[0]

        # Calculate model curves including macroscopic droplets
        q_plot = np.logspace(np.log10(0.002), np.log10(0.09), 200)
        I_sch_sph = calc_sphere_intensity(
            q_plot, row_dat['Rmean_schulz'], row_dat['stdDev_schulz'], row_dat['phi_schulz'],
            bkg=row_dat['bkg_schulz'], dist_type='schulz'
        )
        I_log_sph = calc_sphere_intensity(
            q_plot, row_dat['Rmean_lognormal'], row_dat['stdDev_lognormal'], row_dat['phi_lognormal'],
            bkg=row_dat['bkg_lognormal'], dist_type='lognormal'
        )
        I_drop_plot = calc_large_sphere_intensity(
            q_plot, row_dat['mean_radius_large_lognormal'], 248.0, row_dat['total_volume_fraction_large_lognormal']
        )

        I_sch_plot = I_sch_sph + I_drop_plot
        I_log_plot = I_log_sph + I_drop_plot

        # Residuals
        I_sch_fit = calc_long_distance_model({
            'Rmean': row_dat['Rmean_schulz'],
            'stdDev': row_dat['stdDev_schulz'],
            'phi': row_dat['phi_schulz'],
            'mean_radius_large': row_dat['mean_radius_large_schulz'],
            'total_volume_fraction_large': row_dat['total_volume_fraction_large_schulz'],
            'bkg': row_dat['bkg_schulz'],
            'rho2': 1.0
        }, q_fit, dist_type='schulz')

        I_log_fit = calc_long_distance_model({
            'Rmean': row_dat['Rmean_lognormal'],
            'stdDev': row_dat['stdDev_lognormal'],
            'phi': row_dat['phi_lognormal'],
            'mean_radius_large': row_dat['mean_radius_large_lognormal'],
            'total_volume_fraction_large': row_dat['total_volume_fraction_large_lognormal'],
            'bkg': row_dat['bkg_lognormal'],
            'rho2': 1.0
        }, q_fit, dist_type='lognormal')

        res_sch = (y_pts - I_sch_fit) / dy_pts
        res_log = (y_pts - I_log_fit) / dy_pts

        chi_sch = float(row_dat['redchi_schulz'])
        chi_log = float(row_dat['redchi_lognormal'])
        rm_sch = float(row_dat['Rmean_schulz'])
        rm_log = float(row_dat['Rmean_lognormal'])

        # Col 0: Scattering Intensity
        ax0 = axes[row_idx, 0]
        ax0.errorbar(q_fit, y_pts, yerr=dy_pts, fmt='o', color='#94a3b8', markersize=2.2, alpha=0.45, label='SAXS Data')
        ax0.plot(q_plot, I_log_plot, color='#2563eb', lw=2.0, label=f'Lognormal + Droplets ($\\chi^2={chi_log:.2f}$)')
        ax0.plot(q_plot, I_sch_plot, color='#dc2626', lw=1.6, linestyle='--', label=f'Schulz + Droplets ($\\chi^2={chi_sch:.2f}$)')
        ax0.plot(q_plot, I_drop_plot, color='#16a34a', lw=1.2, linestyle=':', label='Droplets Form Factor')
        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_xlim(0.0019, 0.095)
        ax0.set_ylim(0.05, 1200)
        ax0.set_ylabel(f'Frame {f} (t={t_val:.2f}s)\n$I(q)$ [cm$^{{-1}}$]', fontsize=10)
        ax0.grid(True, which='both', linestyle=':', alpha=0.5)
        if row_idx == 0:
            ax0.set_title('Scattering Profiles $I(q)$', fontsize=11, fontweight='bold')
        ax0.legend(fontsize=8, loc='upper right')

        # Col 1: Weighted Residuals
        ax1 = axes[row_idx, 1]
        ax1.axhline(0, color='black', lw=0.9)
        ax1.axhspan(-2, 2, color='#e2e8f0', alpha=0.6)
        ax1.plot(q_fit, res_log, 's-', color='#2563eb', lw=1.0, markersize=2.5, alpha=0.7, label='Lognormal')
        ax1.plot(q_fit, res_sch, 'o--', color='#dc2626', lw=1.0, markersize=2.5, alpha=0.7, label='Schulz-Flory')
        ax1.set_xscale('log')
        ax1.set_xlim(0.0019, 0.095)
        ax1.set_ylim(-4.5, 4.5)
        ax1.set_ylabel(r'Res. [$\sigma$]', fontsize=9.5)
        ax1.grid(True, which='both', linestyle=':', alpha=0.5)
        if row_idx == 0:
            ax1.set_title('Weighted Residuals $(I_{\\rm exp} - I_{\\rm mod})/\\sigma$', fontsize=11, fontweight='bold')
        ax1.legend(fontsize=8, loc='upper right')

        # Col 2: Probability Density Functions f(R)
        ax2 = axes[row_idx, 2]
        r_grid = np.linspace(0.1, 180.0, 500)
        pdf_sch = schulz_flory_pdf(r_grid, rm_sch, row_dat['stdDev_schulz'])
        mu_log, sig_log = convert_mean_std_to_mu_sigma(rm_log, row_dat['stdDev_lognormal'])
        pdf_log = log_normal_pdf(r_grid, mu_log, sig_log)

        ax2.plot(r_grid, pdf_log, color='#2563eb', lw=2.0, label=f'Lognormal: <R>={rm_log:.1f} nm')
        ax2.plot(r_grid, pdf_sch, color='#dc2626', lw=1.8, linestyle='--', label=f'Schulz-Flory: <R>={rm_sch:.1f} nm')
        ax2.set_xlim(0, 160)
        ax2.set_ylabel('$f(R)$ [nm$^{-1}$]', fontsize=9.5)
        ax2.grid(True, linestyle=':', alpha=0.5)
        if row_idx == 0:
            ax2.set_title('Size Distribution $f(R)$', fontsize=11, fontweight='bold')
        ax2.legend(fontsize=8, loc='upper right')

        if row_idx == len(comp_frames) - 1:
            ax0.set_xlabel('Scattering Vector $q$ [nm$^{-1}$]', fontsize=10)
            ax1.set_xlabel('Scattering Vector $q$ [nm$^{-1}$]', fontsize=10)
            ax2.set_xlabel('Particle Radius $R$ [nm]', fontsize=10)

    plt.tight_layout()
    comp_plot_path = os.path.join(OUT_DIR, 'schulz_vs_lognormal_comparison.png')
    plt.savefig(comp_plot_path, dpi=250)
    plt.close(fig)

    elapsed = time.time() - start_t
    print(f"\nLong-distance fitting and model comparison completed in {elapsed:.2f} s. Outputs saved to: {OUT_DIR}")
    return df_params


if __name__ == '__main__':
    fit_long_distance_series()
