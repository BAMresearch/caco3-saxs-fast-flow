#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Fit Composite Baseline Scattering Curve for ND-i-001.
----------------------------------------------------------------
Establishes the global structural baseline across q in [0.002, 2.5] nm^-1:
- Primary building blocks: Rmean, stdDev, phi, D, xi
- Large droplet population: mean_radius_large, std_radius_large, total_volume_fraction_large
- Viscoelastic network fluctuations (Fisher-Burford / Tanaka): A_fluc, xi_fluc, mu_fluc = 4.0
- Incoherent background: bkg

Outputs:
  output/composite/baseline_parameters.json
  output/composite/baseline_parameters.csv
  output/composite/composite_baseline_fit.png
  output/composite/composite_baseline_fit.pdf
"""

import os
import sys
import json
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
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_fb_intensity,
    RHO2_CALIBRATED
)
from models.fitting_engine import residuals_func


def fit_composite_baseline():
    print("=" * 80)
    print("  COMPOSITE BASELINE SCATTERING CURVE FITTING (ND-i-001, 40 mM CaCO3, pH 11.3, Pure)")
    print("=" * 80)

    data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi001_fitting_results_composite.csv')
    out_dir = os.path.join(PROJECT_DIR, 'output', 'composite')
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(data_path)
    mask = (df['q'] >= 0.002) & (df['q'] <= 2.5) & (df['I'] > 0)
    q_data = df.loc[mask, 'q'].values
    y_data = df.loc[mask, 'I'].values
    dy_data = np.maximum(df.loc[mask, 'dI'].values, 1e-4)

    print(f"Loaded {len(q_data)} points from composite dataset across q: [{q_data.min():.5f}, {q_data.max():.5f}] nm^-1")

    # Load short-distance kinetics endpoint (Frame 49) to ground high-q parameters
    sd_params_path = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    df_sd = pd.read_csv(sd_params_path)
    row49 = df_sd.iloc[-1]
    print(f"Grounded to short-distance kinetics frame {int(row49['frame_index'])} ({row49['filename']}) at t = {row49['time_seconds']:.2f} s")

    # Initial parameter object grounded to short-distance kinetics
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)
    # Ground high-q building block and fluctuation components to Frame 49
    p.add('Rmean', value=float(row49['Rmean']), vary=False)
    p.add('stdDev', value=float(row49['stdDev']), vary=False)
    p.add('D', value=float(row49['D']), vary=False)
    p.add('A_fluc', value=float(row49['A_fluc']), vary=False)
    p.add('xi_fluc', value=float(row49['xi_fluc']), vary=False)
    p.add('mu_fluc', value=4.0, vary=False)
    p.add('bkg', value=float(row49['bkg']), vary=False)

    # Fit low-q and multiscale network parameters on composite curve
    p.add('phi', value=float(row49['phi']), min=1e-5, max=0.01, vary=True)
    p.add('xi', value=50000.0, min=100.0, max=100000.0, vary=True)
    p.add('mean_radius_large', value=3195.0, min=500.0, max=10000.0, vary=True)
    p.add('std_radius_large', value=352.0, min=10.0, max=1500.0, vary=True)
    p.add('total_volume_fraction_large', value=2.6e-5, min=1e-7, max=0.001, vary=True)

    print("\n--- Running Global Refinement for Low-q Droplet and Network Parameters ---")
    out_cov = minimize(residuals_func, p, args=(q_data, y_data, dy_data, True, False, 1.0, 100), method='least_squares', calc_covar=True)
    final_params = out_cov.params
    red_chi2 = out_cov.redchi

    print("\n" + "=" * 80)
    print(f"  BASELINE FIT CONVERGED: Reduced Chi-Square = {red_chi2:.4f}")
    print("=" * 80)

    param_dict = {}
    for p_name in final_params:
        p_val = float(final_params[p_name].value)
        p_err = float(final_params[p_name].stderr) if final_params[p_name].stderr is not None else 0.0
        param_dict[p_name] = p_val
        param_dict[f"std_err_{p_name}"] = p_err

    # Transfer exact standard errors for short-distance grounded parameters
    for k in ['Rmean', 'stdDev', 'D', 'A_fluc', 'xi_fluc', 'bkg']:
        err_col = f"std_err_{k}"
        if err_col in row49:
            param_dict[err_col] = float(row49[err_col])

    for p_name in final_params:
        print(f"  {p_name:<28}: {param_dict[p_name]:12.6f} ± {param_dict[f'std_err_{p_name}']:10.6f}")

    param_dict['redchi'] = float(red_chi2)
    param_dict['aic'] = float(out_cov.aic)
    param_dict['bic'] = float(out_cov.bic)

    json_path = os.path.join(out_dir, 'baseline_parameters.json')
    with open(json_path, 'w') as f:
        json.dump(param_dict, f, indent=4)
    print(f"\nSaved baseline parameters JSON: {json_path}")

    csv_path = os.path.join(out_dir, 'baseline_parameters.csv')
    pd.DataFrame([param_dict]).to_csv(csv_path, index=False)
    print(f"Saved baseline parameters CSV: {csv_path}")

    # --- Plot Baseline Fit & Component Decomposition ---
    print("\nGenerating composite baseline fit figure...")
    q_dense = np.logspace(np.log10(q_data.min()), np.log10(q_data.max()), 800)

    I_prim = calc_primary_fractal_intensity(
        q_dense, final_params['Rmean'].value, final_params['stdDev'].value,
        final_params['phi'].value, final_params['D'].value, final_params['xi'].value,
        rho2=RHO2_CALIBRATED, integration_steps=200
    )
    I_large = calc_large_sphere_intensity(
        q_dense, final_params['mean_radius_large'].value, final_params['std_radius_large'].value,
        final_params['total_volume_fraction_large'].value, rho2=RHO2_CALIBRATED, integration_steps=200
    )
    I_fluc = calc_fb_intensity(
        q_dense, final_params['A_fluc'].value, final_params['xi_fluc'].value,
        mu_fluc=final_params.get('mu_fluc', 4.0).value if 'mu_fluc' in final_params else 4.0
    )
    I_bkg = np.full_like(q_dense, final_params['bkg'].value)
    I_tot = I_prim + I_large + I_fluc + I_bkg

    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.12)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    # Top Panel: Scattering Curve and Model Components
    ax_fit.loglog(q_data, y_data, 'o', color='#334155', markersize=2.5, alpha=0.45, label='Experimental Composite Data')
    ax_fit.loglog(q_dense, I_tot, color='#dc2626', linestyle='-', linewidth=2.2, label=f'Full Multiscale Fit ($\\chi^2_{{\\mathrm{{red}}}} = {red_chi2:.2f}$)')
    ax_fit.loglog(q_dense, I_prim, color='#059669', linestyle='--', linewidth=1.8,
                  label=f'Primary Fractal Network ($\\langle R \\rangle = {final_params["Rmean"].value:.2f} \\pm {final_params["stdDev"].value:.2f}$ nm, $D_f = {final_params["D"].value:.2f}$)')
    ax_fit.loglog(q_dense, I_large, color='#d97706', linestyle='-.', linewidth=1.8,
                  label=f'Secondary Droplets ($\\langle R \\rangle = {final_params["mean_radius_large"].value:.0f}$ nm, $\\phi = {final_params["total_volume_fraction_large"].value:.2e}$)')
    ax_fit.loglog(q_dense, I_fluc, color='#9333ea', linestyle=':', linewidth=1.8,
                  label=f'Tanaka Viscoelastic Fluctuations ($\\xi_{{\\mathrm{{fluc}}}} = {final_params["xi_fluc"].value:.2f}$ nm, $A_{{\\mathrm{{fluc}}}} = {final_params["A_fluc"].value:.3f}$ cm$^{{-1}}$)')
    ax_fit.axhline(final_params['bkg'].value, color='#64748b', linestyle=':', linewidth=1.2,
                   label=f'Incoherent Background ({final_params["bkg"].value:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-001: Full-Range Multiscale Composite Baseline Fit (40 mM CaCO$_3$, pH 11.3, Pure)\n'
                     r'$q \in [0.002, 2.5]\ \mathrm{nm}^{-1}$ across 4 length scale decades', fontsize=12, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=9.0, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[y_data > 0].min() * 0.5), top=y_data.max() * 3.0)

    # Bottom Panel: Normalized Residuals
    y_calc_data = calc_model_intensity(final_params, q_data, include_large_spheres=True, integration_steps=120)
    residuals = (y_data - y_calc_data) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color='#d90429', markersize=2.2, alpha=0.6)
    ax_res.set_xscale('log')
    ax_res.set_xlabel('Scattering Vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Normalized Res. [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plot_path_png = os.path.join(out_dir, 'composite_baseline_fit.png')
    plot_path_pdf = os.path.join(out_dir, 'composite_baseline_fit.pdf')
    plt.savefig(plot_path_png, dpi=300)
    plt.savefig(plot_path_pdf)
    plt.close()
    print(f"Saved baseline fit figures: {plot_path_png} and {plot_path_pdf}")
    print("=" * 80)


if __name__ == '__main__':
    fit_composite_baseline()
