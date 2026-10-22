#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Fit Composite Baseline Scattering Curve for ND-i-006.
----------------------------------------------------------------
Establishes the global structural baseline across q in [0.002, 2.5] nm^-1
for 40.0 mM CaCO3 (pH 9.75):
- Primary Building Blocks in Mass Fractal Network: R_prim, std_prim, phi_net, Dm, xi
- Secondary Large Droplets / Large Spheres: mean_radius_large, std_radius_large, total_volume_fraction_large
- Macroscopic Porod Surface Tail: A_porod / q^4
- Incoherent background: bkg
- NO ACC PARTICLES (experimentally absent in this condition)

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

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_primary_fractal_intensity,
    calc_large_sphere_intensity,
    calc_porod_intensity,
    RHO2_CALIBRATED
)


def fit_composite_baseline():
    print("=" * 80)
    print("  COMPOSITE BASELINE SCATTERING CURVE FITTING (ND-i-006, 40.0 mM CaCO3, pH 9.75)")
    print("  Model: Spinodal Mass Fractal Network + Large Spheres / Droplets (No ACC Particles)")
    print("=" * 80)

    data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi006_fitting_results_composite.csv')
    out_dir = os.path.join(PROJECT_DIR, 'output', 'composite')
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(data_path)
    mask = (df['q'] >= 0.002) & (df['q'] <= 2.5) & (df['I'] > 0)
    q_data = df.loc[mask, 'q'].values
    y_data = df.loc[mask, 'I'].values
    dy_data = np.maximum(df.loc[mask, 'dI'].values, 1e-4)

    print(f"Loaded {len(q_data)} points from composite dataset across q: [{q_data.min():.5f}, {q_data.max():.5f}] nm^-1")

    dy_fit = np.maximum(dy_data, 0.04 * y_data)

    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)

    # 1. Mesoscale Spinodal Fractal Network & Building Blocks
    p.add('R_prim', value=3.0, min=1.5, max=15.0, vary=True)
    p.add('std_prim', value=2.8, min=0.5, max=10.0, vary=True)
    p.add('phi_net', value=2.5e-5, min=1e-7, max=0.01, vary=True)
    p.add('Dm', value=2.00, min=1.2, max=2.90, vary=True)
    p.add('xi', value=800.0, min=50.0, max=30000.0, vary=True)

    # 2. Secondary Large Spheres (Macroscopic droplets / large spheres)
    p.add('mean_radius_large', value=3100.0, min=1000.0, max=8000.0, vary=True)
    p.add('std_radius_large', value=270.0, min=20.0, max=1500.0, vary=True)
    p.add('total_volume_fraction_large', value=1.1e-4, min=1e-8, max=0.005, vary=True)

    # 3. Macroscopic Porod Surface Tail & Incoherent Background
    p.add('A_porod', value=1.0e-6, min=0.0, max=0.01, vary=True)
    p.add('bkg', value=2.3e-4, min=1e-6, max=0.01, vary=True)

    def model_eval(pars, q_eval):
        I_net = calc_primary_fractal_intensity(
            q_eval, pars['R_prim'].value, pars['std_prim'].value,
            pars['phi_net'].value, pars['Dm'].value, pars['xi'].value,
            rho2=RHO2_CALIBRATED, integration_steps=80
        )
        I_large = calc_large_sphere_intensity(
            q_eval, pars['mean_radius_large'].value, pars['std_radius_large'].value,
            pars['total_volume_fraction_large'].value,
            rho2=RHO2_CALIBRATED, integration_steps=80
        )
        I_por = calc_porod_intensity(q_eval, pars['A_porod'].value, 4.0)
        return I_net + I_large + I_por + pars['bkg'].value

    print("\n--- Running Stage 1: Robust Log-Space Global Pre-alignment ---")
    def res_log(pars):
        c = model_eval(pars, q_data)
        return np.log10(np.maximum(y_data, 1e-6)) - np.log10(np.maximum(c, 1e-10))

    out1 = minimize(res_log, p, method='least_squares', max_nfev=150)
    p = out1.params.copy()
    print(f"Stage 1 Complete: R_prim = {p['R_prim'].value:.2f} nm, Dm = {p['Dm'].value:.3f}, xi = {p['xi'].value:.1f} nm, R_large = {p['mean_radius_large'].value:.1f} nm")

    print("\n--- Running Stage 2: Weighted Physical Least-Squares Refinement ---")
    def res_phys(pars):
        c = model_eval(pars, q_data)
        return (y_data - c) / dy_fit

    out2 = minimize(res_phys, p, method='least_squares', max_nfev=250)
    p = out2.params.copy()
    print(f"Stage 2 Complete: R_prim = {p['R_prim'].value:.2f} nm, Dm = {p['Dm'].value:.3f}, xi = {p['xi'].value:.1f} nm, R_large = {p['mean_radius_large'].value:.1f} nm")

    print("\n--- Running Stage 3: Exact Covariance Matrix & Standard Errors ---")
    p_cov = p.copy()
    for k in ['R_prim', 'std_prim', 'phi_net', 'Dm', 'xi', 'mean_radius_large', 'std_radius_large', 'total_volume_fraction_large', 'A_porod', 'bkg']:
        p_cov[k].vary = True

    try:
        out_cov = minimize(res_phys, p_cov, method='leastsq', calc_covar=True)
        final_params = out_cov.params
    except Exception as e:
        print(f"Covariance warning: {e}, falling back to least_squares parameters")
        final_params = out2.params

    c_final = model_eval(final_params, q_data)
    dof = max(1, len(q_data) - 9)
    red_chi2_phys = np.sum(((y_data - c_final) / dy_fit) ** 2) / dof
    red_chi2_raw = np.sum(((y_data - c_final) / dy_data) ** 2) / dof

    print("\n" + "=" * 80)
    print(f"  BASELINE FIT CONVERGED: Reduced Chi-Square (phys) = {red_chi2_phys:.4f}")
    print("=" * 80)

    param_dict = {}
    for p_name in final_params:
        p_val = float(final_params[p_name].value)
        p_err = float(final_params[p_name].stderr) if final_params[p_name].stderr is not None else 0.0
        param_dict[p_name] = p_val
        param_dict[f"std_err_{p_name}"] = p_err
        print(f"  {p_name:<28}: {p_val:12.6f} ± {p_err:10.6f}")

    # Aliases
    param_dict['Rmean'] = param_dict['R_prim']
    param_dict['stdDev'] = param_dict['std_prim']
    param_dict['D'] = param_dict['Dm']
    param_dict['phi'] = param_dict['phi_net']
    param_dict['redchi'] = float(red_chi2_phys)
    param_dict['redchi_raw'] = float(red_chi2_raw)

    json_path = os.path.join(out_dir, 'baseline_parameters.json')
    with open(json_path, 'w') as f:
        json.dump(param_dict, f, indent=4)
    print(f"\nSaved baseline parameters JSON: {json_path}")

    csv_path = os.path.join(out_dir, 'baseline_parameters.csv')
    pd.DataFrame([param_dict]).to_csv(csv_path, index=False)
    print(f"Saved baseline parameters CSV: {csv_path}")

    # Plot
    print("\nGenerating composite baseline fit figure...")
    q_dense = np.logspace(np.log10(q_data.min()), np.log10(q_data.max()), 600)

    I_net = calc_primary_fractal_intensity(
        q_dense, final_params['R_prim'].value, final_params['std_prim'].value,
        final_params['phi_net'].value, final_params['Dm'].value, final_params['xi'].value,
        rho2=RHO2_CALIBRATED, integration_steps=80
    )
    I_large = calc_large_sphere_intensity(
        q_dense, final_params['mean_radius_large'].value, final_params['std_radius_large'].value,
        final_params['total_volume_fraction_large'].value, rho2=RHO2_CALIBRATED, integration_steps=80
    )
    I_por = calc_porod_intensity(q_dense, final_params['A_porod'].value, 4.0)
    I_bkg = np.full_like(q_dense, final_params['bkg'].value)
    I_tot = I_net + I_large + I_por + I_bkg

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color='#334155', markersize=2.5, alpha=0.40, label='Experimental Composite Data')
    ax_fit.loglog(q_dense, I_tot, 'r-', linewidth=2.2, label=f'Full Unified Fit ($\\chi^2_{{\\mathrm{{red}}}} = {red_chi2_phys:.2f}$)')
    ax_fit.loglog(q_dense, I_net, 'g--', linewidth=1.6,
                  label=f'Spinodal Fractal Network ($R_{{\\rm prim}} = {final_params["R_prim"].value:.1f}$ nm, $D_m = {final_params["Dm"].value:.2f}$, $\\xi = {final_params["xi"].value:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_large, color='#0284c7', linestyle='-.', linewidth=1.6,
                  label=f'Large Spheres / Droplets ($\\langle R \\rangle = {final_params["mean_radius_large"].value:.0f}$ nm, $\\phi = {final_params["total_volume_fraction_large"].value:.2e}$)')
    ax_fit.loglog(q_dense, I_por, color='#d97706', linestyle=':', linewidth=1.5,
                  label=f'Porod Surface Tail ($A = {final_params["A_porod"].value*1e6:.2f}\\times 10^{{-6}}$)')
    ax_fit.axhline(final_params['bkg'].value, color='#8d99ae', linestyle=':', linewidth=1.2,
                   label=f'Incoherent Background ({final_params["bkg"].value:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-006: Multiscale Composite Baseline Fit (40.0 mM CaCO$_3$, pH 9.75)\n'
                     r'Spinodal Mass Fractal Network + Large Droplets / Spheres ($q \in [0.002, 2.5]\ \mathrm{nm}^{-1}$)', fontsize=11, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.set_xlim(0.0018, 2.65)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[y_data > 0].min() * 0.5), top=y_data.max() * 3.0)

    residuals = (y_data - c_final) / dy_fit
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color='#dc2626', markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.65)
    ax_res.set_xlabel('Scattering Vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Res. [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'composite_baseline_fit.png')
    pdf_path = os.path.join(out_dir, 'composite_baseline_fit.pdf')
    plt.savefig(plot_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved baseline fit figure: {plot_path} & {pdf_path}")
    print("=" * 80)


if __name__ == '__main__':
    fit_composite_baseline()
