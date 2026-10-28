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

    # Initial parameter object grounded to short-distance kinetics
    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)
    p.add('phi', value=0.00024, min=1e-7, max=0.1, vary=True)
    p.add('Rmean', value=6.24, min=2.0, max=15.0, vary=True)
    p.add('stdDev', value=2.47, min=0.2, max=10.0, vary=True)
    p.add('D', value=1.86, min=1.1, max=2.999, vary=True)
    p.add('xi', value=735.0, min=50.0, max=50000.0, vary=False)
    p.add('bkg', value=0.00087, min=1e-6, max=1.0, vary=False)
    p.add('mean_radius_large', value=3042.0, min=500.0, max=10000.0, vary=False)
    p.add('std_radius_large', value=278.0, min=10.0, max=1500.0, vary=False)
    p.add('total_volume_fraction_large', value=7.7e-6, min=1e-8, max=0.01, vary=False)
    p.add('A_fluc', value=0.0054, min=1e-5, max=5.0, vary=False)
    p.add('xi_fluc', value=3.05395, min=0.5, max=15.0, vary=False)
    p.add('mu_fluc', value=4.0, vary=False)

    # --- 4-Stage Hierarchical Fit ---
    print("\n--- Running Stage 1: Primary Building Block & Fractal Fit (q in [0.015, 0.4] nm^-1) ---")
    s1_mask = (q_data >= 0.015) & (q_data <= 0.4)
    out1 = minimize(residuals_func, p, args=(q_data[s1_mask], y_data[s1_mask], dy_data[s1_mask], False, False, 1.0, 80), method='least_squares')
    p = out1.params.copy()
    print(f"Stage 1 Complete: Rmean = {p['Rmean'].value:.3f} nm, D = {p['D'].value:.3f}, phi = {p['phi'].value:.6e}")

    print("\n--- Running Stage 2: Large Droplet Population Fit at Low-q (q < 0.01 nm^-1) ---")
    p['phi'].vary = False
    p['Rmean'].vary = False
    p['stdDev'].vary = False
    p['D'].vary = False
    p['xi'].vary = True
    p['mean_radius_large'].vary = True
    p['std_radius_large'].vary = True
    p['total_volume_fraction_large'].vary = True
    s2_mask = (q_data < 0.01)
    out2 = minimize(residuals_func, p, args=(q_data[s2_mask], y_data[s2_mask], dy_data[s2_mask], True, False, 1.0, 80), method='least_squares')
    p = out2.params.copy()
    print(f"Stage 2 Complete: R_large = {p['mean_radius_large'].value:.1f} nm, phi_large = {p['total_volume_fraction_large'].value:.6e}, xi = {p['xi'].value:.1f} nm")

    print("\n--- Running Stage 3: Viscoelastic Network Fluctuations Fit at High-q (q > 0.3 nm^-1) ---")
    p['xi'].vary = False
    p['mean_radius_large'].vary = False
    p['std_radius_large'].vary = False
    p['total_volume_fraction_large'].vary = False
    p['A_fluc'].vary = True
    p['xi_fluc'].vary = False  # Grounded to short-distance kinetics baseline xi_fluc = 3.05 nm
    p['bkg'].vary = True
    s3_mask = (q_data > 0.3)
    out3 = minimize(residuals_func, p, args=(q_data[s3_mask], y_data[s3_mask], dy_data[s3_mask], False, False, 1.0, 80), method='least_squares')
    p = out3.params.copy()
    print(f"Stage 3 Complete: A_fluc = {p['A_fluc'].value:.4f} cm^-1, xi_fluc = {p['xi_fluc'].value:.3f} nm, bkg = {p['bkg'].value:.6e}")

    print("\n--- Running Stage 4: Full Multi-Component Global Refinement Across All Length Scales ---")
    for k in p:
        if k not in ['rho2', 'mu_fluc', 'xi_fluc']:
            p[k].vary = True
    p['xi_fluc'].vary = False  # Maintain grounded correlation length (3.05 nm)
    p['Rmean'].min = 4.0
    p['Rmean'].max = 8.5
    p['stdDev'].min = 1.5
    p['stdDev'].max = 4.0
    p['bkg'].min = 0.0005
    p['bkg'].max = 0.0015
    out_final = minimize(residuals_func, p, args=(q_data, y_data, dy_data, True, False, 1.0, 100), method='least_squares', max_nfev=600)
    
    # Compute exact covariance matrix and standard errors around converged solution
    p_cov = Parameters()
    p_cov.add('rho2', value=RHO2_CALIBRATED, vary=False)
    for k in ['phi', 'D', 'xi', 'mean_radius_large', 'std_radius_large', 'total_volume_fraction_large', 'A_fluc']:
        p_cov.add(k, value=out_final.params[k].value, min=1e-8, max=1e5, vary=True)
    p_cov.add('Rmean', value=out_final.params['Rmean'].value, min=4.0, max=8.5, vary=True)
    p_cov.add('stdDev', value=out_final.params['stdDev'].value, min=1.5, max=4.0, vary=True)
    p_cov.add('bkg', value=out_final.params['bkg'].value, min=0.0005, max=0.0015, vary=True)
    p_cov.add('xi_fluc', value=out_final.params['xi_fluc'].value, vary=False)
    p_cov.add('mu_fluc', value=4.0, vary=False)
    out_cov = minimize(residuals_func, p_cov, args=(q_data, y_data, dy_data, True, False, 1.0, 100), method='least_squares', calc_covar=True)
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
        print(f"  {p_name:<28}: {p_val:12.6f} ± {p_err:10.6f}")

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
