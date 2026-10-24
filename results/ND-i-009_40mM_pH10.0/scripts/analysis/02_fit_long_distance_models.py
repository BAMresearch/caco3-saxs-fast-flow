#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 02: Forward Model Fitting for Long-Distance Dataset (ND-i-009).
---------------------------------------------------------------------
Physical Formulation & Investigation:
1. Frame 0 Physics:
   - Baseline pre-reaction state with NO large spheres (has_large_spheres=False, phi_large=0, R_large=0).
   - xi -> infinity and constant.
2. Large Sphere Polydispersity:
   - Polydispersity sigma/R is fixed to a constant physical ratio (sigma/R = 0.065).
   - Prevents unphysical frame-to-frame fluctuations in droplet dispersion.
   - Eliminates parameter cross-talk between the droplet Porod decay and the mass fractal slope q^-D.
3. Fractal Dimension D Continuity:
   - Regularized to prevent unphysical frame-to-frame spikes (|dD| < 0.015).
4. Physical Real Time:
   - t = 0.045 + 0.015 * frame [s].

Outputs saved to:
  output/long_distance/option_c_forward/
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
from lmfit import minimize, Parameters

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import calc_model_intensity, decompose_model_components, RHO2_CALIBRATED
from models.fitting_engine import check_frame_signal_quality

DATA_PATH = os.path.join(PROJECT_DIR, 'data', 'dataNDi009_fitting_results_long_distance.csv')
SHORT_CSV = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'viscoelastic', 'parameters.csv')
BASELINE_JSON = os.path.join(PROJECT_DIR, 'output', 'merged', 'baseline_parameters.json')

CONST_POLY_RATIO = 0.065  # Fixed physical relative polydispersity sigma_large / R_large


def load_dataset_and_transferred_params():
    df_data = pd.read_csv(DATA_PATH)
    q_data = df_data['q'].values
    frame_cols = [c for c in df_data.columns if c not in ['Unnamed: 0', 'q'] and not c.startswith('error_')]
    frame_indices = sorted([int(c) for c in frame_cols])

    transferred_Rmean = 8.51
    transferred_stdDev = 3.34
    if os.path.exists(SHORT_CSV):
        df_short = pd.read_csv(SHORT_CSV)
        transferred_Rmean = float(df_short['Rmean'].mean())
        transferred_stdDev = float(df_short['stdDev'].mean())
    elif os.path.exists(BASELINE_JSON):
        with open(BASELINE_JSON, 'r') as f:
            b_json = json.load(f)
        transferred_Rmean = float(b_json.get('Rmean', 8.51))
        transferred_stdDev = float(b_json.get('stdDev', 3.34))

    return df_data, q_data, frame_indices, transferred_Rmean, transferred_stdDev


def create_forward_params(Rmean, stdDev, seed=None):
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)
    p.add('Rmean', value=Rmean, vary=False)
    p.add('stdDev', value=stdDev, vary=False)
    p.add('xi', value=1e7, vary=False)  # Infinite and constant!
    p.add('bkg', value=0.0, vary=False)

    phi_val = 0.001 if seed is None else seed['phi'].value
    d_val = 1.45 if seed is None else seed['D'].value
    r_val = 2800.0 if seed is None else seed['mean_radius_large'].value
    phil_val = 4e-5 if seed is None else seed['total_volume_fraction_large'].value

    p.add('phi', value=phi_val, min=1e-5, max=0.05, vary=True)
    p.add('D', value=d_val, min=1.2, max=2.3, vary=True)
    p.add('mean_radius_large', value=r_val, min=2000.0, max=3800.0, vary=True)
    # Lock polydispersity to constant physical ratio
    p.add('poly_ratio', value=CONST_POLY_RATIO, vary=False)
    p.add('std_radius_large', expr=f'mean_radius_large * {CONST_POLY_RATIO}')
    p.add('total_volume_fraction_large', value=phil_val, min=1e-8, max=0.005, vary=True)
    return p


def residuals_forward_gradual(params, q, y, dy, prev_p=None, f=1):
    calc = calc_model_intensity(params, q, include_large_spheres=True, include_oz=False, integration_steps=90)
    res = (y - calc) / dy

    if prev_p is not None:
        w_R = 1.0 if f <= 30 else 0.5
        dR = (params['mean_radius_large'].value - prev_p['mean_radius_large'].value) / 35.0
        # Strict continuity on D: eliminates unphysical spikes
        dD = (params['D'].value - prev_p['D'].value) / 0.012
        dphi = (params['phi'].value - prev_p['phi'].value) / max(1e-4, prev_p['phi'].value * 0.15)
        dphil = (params['total_volume_fraction_large'].value - prev_p['total_volume_fraction_large'].value) / max(1e-6, prev_p['total_volume_fraction_large'].value * 0.25)

        res = np.append(res, [w_R * dR, 2.2 * dD, 0.4 * dphi, 0.3 * dphil])
    return res


def plot_summary_trajectories(df_params, save_path):
    time_sec = df_params['time_seconds'].values

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("ND-i-009 Long Distance (31.0 m): Option C Forward Kinetics (Constant Polydispersity)", fontsize=13, fontweight='bold')

    df_active = df_params[df_params['frame'] >= 1].copy()
    t_active = df_active['time_seconds'].values

    # Panel A: Large Sphere Kinetics
    ax_a = axes[0, 0]
    ax_a.plot(t_active, df_active['mean_radius_large'], 'o-', color='#1f77b4', markersize=2.5, linewidth=1.8, label=r'Droplet Radius $\langle R_{\mathrm{large}} \rangle$')
    ax_a.fill_between(t_active,
                      df_active['mean_radius_large'] - df_active['std_err_mean_radius_large'],
                      df_active['mean_radius_large'] + df_active['std_err_mean_radius_large'],
                      color='#1f77b4', alpha=0.35, label=r'Fit Uncertainty $\pm 1\sigma$')
    ax_a.fill_between(t_active,
                      df_active['mean_radius_large'] - df_active['std_radius_large'],
                      df_active['mean_radius_large'] + df_active['std_radius_large'],
                      color='#1f77b4', alpha=0.12, label=r'Polydispersity ($\sigma/R = 0.065$)')
    ax_a.set_ylabel(r'Radius $R_{\mathrm{large}}$ [nm]', fontsize=11)
    ax_a.set_xlabel('Reaction Time [s]', fontsize=11)
    ax_a.set_title(r'(A) Large Sphere Population Kinetics $\pm 1\sigma$', fontsize=11, fontweight='bold')
    ax_a.grid(True, linestyle='--', alpha=0.5)
    ax_a.legend(loc='best', fontsize=9.0)

    # Panel B: Volume Fractions
    ax_b = axes[0, 1]
    ax_b.plot(time_sec, df_params['phi'] * 1e3, 's-', color='#2ca02c', markersize=2.5, linewidth=1.8, label=r'Mass Fractal $\phi_{\mathrm{agg}} \times 10^3$')
    ax_b.fill_between(time_sec,
                      (df_params['phi'] - df_params['std_err_phi']) * 1e3,
                      (df_params['phi'] + df_params['std_err_phi']) * 1e3,
                      color='#2ca02c', alpha=0.22, label=r'$\phi_{\mathrm{agg}} \pm 1\sigma$ uncertainty')
    ax_b.plot(time_sec, df_params['total_volume_fraction_large'] * 1e3, '^-', color='#d62728', markersize=2.5, linewidth=1.8, label=r'Large Spheres $\phi_{\mathrm{large}} \times 10^3$')
    ax_b.fill_between(time_sec,
                      (df_params['total_volume_fraction_large'] - df_params['std_err_total_volume_fraction_large']) * 1e3,
                      (df_params['total_volume_fraction_large'] + df_params['std_err_total_volume_fraction_large']) * 1e3,
                      color='#d62728', alpha=0.22, label=r'$\phi_{\mathrm{large}} \pm 1\sigma$ uncertainty')
    ax_b.set_ylabel(r'Volume Fraction $\phi \times 10^3$', fontsize=11)
    ax_b.set_xlabel('Reaction Time [s]', fontsize=11)
    ax_b.set_title(r'(B) Phase Volume Partitioning $\pm 1\sigma$', fontsize=11, fontweight='bold')
    ax_b.grid(True, linestyle='--', alpha=0.5)
    ax_b.legend(loc='best', fontsize=9.0)

    # Panel C: Mesoscale Structure (D and xi)
    ax_c = axes[1, 0]
    ax_c.plot(time_sec, df_params['D'], 'd-', color='#9467bd', markersize=2.5, linewidth=1.8, label=r'Fractal Dimension $D_f$')
    ax_c.fill_between(time_sec,
                      df_params['D'] - df_params['std_err_D'],
                      df_params['D'] + df_params['std_err_D'],
                      color='#9467bd', alpha=0.25, label=r'$D_f \pm 1\sigma$ uncertainty')
    ax_c.set_ylabel(r'Fractal Dimension $D_f$', fontsize=11, color='#9467bd')
    ax_c.tick_params(axis='y', labelcolor='#9467bd')
    ax_c.set_ylim(1.2, 2.3)
    ax_c.set_xlabel('Reaction Time [s]', fontsize=11)
    ax_c.set_title(r'(C) Mesoscale Network Structure ($D_f \pm 1\sigma, \xi$)', fontsize=11, fontweight='bold')
    ax_c.grid(True, linestyle='--', alpha=0.5)
    ax_c.legend(loc='lower right', fontsize=9.0)

    ax_c2 = ax_c.twinx()
    ax_c2.axhline(1e4, color='#ff7f0e', linestyle='--', linewidth=1.8, label=r'Infinite Cutoff $\xi \to \infty$')
    ax_c2.set_ylabel(r'Correlation Cutoff $\xi$ [nm]', fontsize=11, color='#ff7f0e')
    ax_c2.tick_params(axis='y', labelcolor='#ff7f0e')
    ax_c2.set_yscale('log')
    ax_c2.set_ylim(1e2, 2e4)

    # Panel D: Fit Quality
    ax_d = axes[1, 1]
    ax_d.plot(time_sec, df_params['redchi'], 'k.-', markersize=3, linewidth=1.5, label=r'$\chi^2_{\mathrm{red}}$')
    ax_d.axhline(1.0, color='gray', linestyle=':', label=r'Ideal Fit $\chi^2_{\mathrm{red}}=1.0$')
    ax_d.set_ylabel(r'Reduced $\chi^2_{\mathrm{red}}$', fontsize=11)
    ax_d.set_xlabel('Reaction Time [s]', fontsize=11)
    ax_d.set_title('(D) Reduced Chi-Square Trajectory', fontsize=11, fontweight='bold')
    ax_d.grid(True, linestyle='--', alpha=0.5)
    ax_d.legend(loc='best', fontsize=9.0)

    plt.tight_layout()
    plt.savefig(save_path, dpi=160)
    plt.close(fig)
    print(f"  Saved trajectory plot: {save_path}")


def run_forward_fitting():
    df_data, q_data, frame_indices, Rmean, stdDev = load_dataset_and_transferred_params()

    out_dir = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'option_c_forward')
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("  LONG-DISTANCE FORWARD FITTING: CONSTANT POLYDISPERSITY (sigma/R = 0.065)")
    print(f"  Frames: {len(frame_indices)} | Rmean={Rmean:.2f} nm, stdDev={stdDev:.2f} nm (fixed)")
    print("=" * 80)

    start_t = time.time()
    sel = (q_data >= 0.002) & (q_data <= 0.22)
    q_fit = q_data[sel]

    # --- Frame 0: Baseline state with NO large spheres ---
    y0 = df_data['0'].values[sel]
    dy0 = np.maximum(df_data['error_0'].values[sel], 1e-4)

    p0 = Parameters()
    p0.add('rho2', value=RHO2_CALIBRATED, vary=False)
    p0.add('phi', value=0.005, min=1e-5, max=0.05, vary=True)
    p0.add('Rmean', value=Rmean, vary=False)
    p0.add('stdDev', value=stdDev, vary=False)
    p0.add('D', value=1.20, min=1.1, max=2.3, vary=True)
    p0.add('xi', value=1e7, vary=False)
    p0.add('bkg', value=0.0, vary=False)
    p0.add('mean_radius_large', value=0.0, vary=False)
    p0.add('std_radius_large', value=0.0, vary=False)
    p0.add('total_volume_fraction_large', value=0.0, vary=False)

    out0 = minimize(lambda p: (y0 - calc_model_intensity(p, q_fit, include_large_spheres=False, include_oz=False)) / dy0, p0, method='least_squares')
    fit_curve_0 = calc_model_intensity(out0.params, q_data, include_large_spheres=False, include_oz=False)

    param_records = []
    fitted_curves = {'q': q_data, '0': fit_curve_0}

    row_0 = {
        'frame': 0,
        'time_seconds': 0.045,
        'phi': float(out0.params['phi'].value),
        'std_err_phi': float(out0.params['phi'].stderr) if out0.params['phi'].stderr else 0.0,
        'D': float(out0.params['D'].value),
        'std_err_D': float(out0.params['D'].stderr) if out0.params['D'].stderr else 0.0,
        'xi': 1e7,
        'std_err_xi': 0.0,
        'has_large_spheres': False,
        'mean_radius_large': 0.0,
        'std_err_mean_radius_large': 0.0,
        'std_radius_large': 0.0,
        'std_err_std_radius_large': 0.0,
        'total_volume_fraction_large': 0.0,
        'std_err_total_volume_fraction_large': 0.0,
        'redchi': float(out0.redchi),
        'aic': float(out0.aic),
        'snr': 0.0,
        'is_valid_signal': False
    }
    param_records.append(row_0)

    # --- Frames > 0: Forward propagation with constant polydispersity ---
    active_frames = [f for f in frame_indices if f > 0]
    p_init = create_forward_params(Rmean, stdDev)
    current_seed = p_init.copy()
    prev_fitted = None

    for idx, f in enumerate(active_frames):
        y_frame = df_data[str(f)].values[sel]
        dy_frame = np.maximum(df_data['error_' + str(f)].values[sel], 1e-4)
        is_valid, snr, _ = check_frame_signal_quality(q_fit, y_frame, dy_frame)

        p_fit_obj = create_forward_params(Rmean, stdDev, seed=current_seed)

        try:
            out = minimize(
                residuals_forward_gradual,
                p_fit_obj,
                args=(q_fit, y_frame, dy_frame, prev_fitted, f),
                method='least_squares',
                max_nfev=250
            )
            current_seed = out.params.copy()
            prev_fitted = out.params.copy()
            redchi = float(out.redchi)
            aic_val = float(out.aic)
        except Exception as e:
            print(f"[Option C Forward] Frame {f} warning: {e}")
            redchi = 99.0
            aic_val = 9999.0

        p = current_seed
        fit_curve = calc_model_intensity(p, q_data, include_large_spheres=True, include_oz=False)
        fitted_curves[str(f)] = fit_curve

        t_sec = 0.045 + 0.015 * f
        row_res = {
            'frame': f,
            'time_seconds': t_sec,
            'phi': float(p['phi'].value),
            'std_err_phi': float(p['phi'].stderr) if p['phi'].stderr else 0.0,
            'D': float(p['D'].value),
            'std_err_D': float(p['D'].stderr) if p['D'].stderr else 0.0,
            'xi': 1e7,
            'std_err_xi': 0.0,
            'has_large_spheres': True,
            'mean_radius_large': float(p['mean_radius_large'].value),
            'std_err_mean_radius_large': float(p['mean_radius_large'].stderr) if p['mean_radius_large'].stderr else 0.0,
            'std_radius_large': float(p['std_radius_large'].value),
            'std_err_std_radius_large': float(p['std_radius_large'].stderr) if p['std_radius_large'].stderr else 0.0,
            'total_volume_fraction_large': float(p['total_volume_fraction_large'].value),
            'std_err_total_volume_fraction_large': float(p['total_volume_fraction_large'].stderr) if p['total_volume_fraction_large'].stderr else 0.0,
            'redchi': redchi,
            'aic': aic_val,
            'snr': float(snr),
            'is_valid_signal': bool(is_valid)
        }
        param_records.append(row_res)

        if (idx + 1) % 50 == 0 or idx == len(active_frames) - 1:
            elapsed = time.time() - start_t
            print(f"  [Option C Forward] Frame {f:03d}/{active_frames[-1]} (t={t_sec:.3f}s) | chi2={redchi:.3f}, D={p['D'].value:.3f}, R_large={p['mean_radius_large'].value:.0f} nm, phi_large={p['total_volume_fraction_large'].value:.2e} [{elapsed:.1f}s]")

    elapsed_total = time.time() - start_t
    print(f"\n[Option C Forward] Completed forward fitting in {elapsed_total:.1f} s.")

    df_params = pd.DataFrame(param_records)
    df_curves = pd.DataFrame(fitted_curves)

    params_csv = os.path.join(out_dir, 'parameters.csv')
    curves_csv = os.path.join(out_dir, 'fitted_curves.csv')
    df_params.to_csv(params_csv, index=False)
    df_curves.to_csv(curves_csv, index=False)
    print(f"  Saved parameters to: {params_csv}")
    print(f"  Saved fitted curves to: {curves_csv}")

    # Generate summary trajectory plot
    traj_path = os.path.join(out_dir, 'trajectories.png')
    plot_summary_trajectories(df_params, traj_path)

    return df_params, df_curves


if __name__ == '__main__':
    run_forward_fitting()

