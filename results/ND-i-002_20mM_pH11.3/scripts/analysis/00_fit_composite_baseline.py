#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Fit Composite Baseline Scattering Curve for ND-i-002.
----------------------------------------------------------------
Establishes the global structural baseline across q in [0.002, 2.54] nm^-1:
- Primary building blocks in Spinodal Fractal Network: Rmean, stdDev, phi, D, xi
- Macroscopic Porod surface tail: A_porod / q^4
- Sub-nanometer clusters / high-q excess: R_pnc, std_pnc, phi_pnc
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
    calc_primary_fractal_intensity,
    calc_pnc_intensity,
    RHO2_CALIBRATED
)


def fit_composite_baseline():
    print("=" * 80)
    print("  COMPOSITE BASELINE SCATTERING CURVE FITTING (ND-i-002, 20 mM CaCO3, pH 11.3)")
    print("=" * 80)

    data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi002_fitting_results_composite.csv')
    out_dir = os.path.join(PROJECT_DIR, 'output', 'composite')
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(data_path)
    mask = (df['q'] >= 0.002) & (df['q'] <= 2.54) & (df['I'] > -0.01)
    q_data = df.loc[mask, 'q'].values
    y_data = df.loc[mask, 'I'].values
    dy_data = np.maximum(df.loc[mask, 'dI'].values, 1e-4)

    print(f"Loaded {len(q_data)} points from composite dataset across q: [{q_data.min():.5f}, {q_data.max():.5f}] nm^-1")

    # 3% systematic uncertainty floor across camera configurations
    dy_fit = np.maximum(dy_data, 0.03 * np.maximum(y_data, 0))

    p = Parameters()
    p.add('rho2', value=RHO2_CALIBRATED, vary=False)

    # 1. Primary building block & spinodal network
    p.add('phi', value=4.5e-5, min=1e-6, max=0.01, vary=True)
    p.add('Rmean', value=11.2, min=5.0, max=25.0, vary=True)
    p.add('stdDev', value=4.0, min=0.5, max=10.0, vary=True)
    p.add('D', value=1.88, min=1.4, max=2.5, vary=True)
    p.add('xi', value=300.0, min=50.0, max=5000.0, vary=True)

    # 2. Macroscopic Porod surface tail
    p.add('A_porod', value=5.5e-8, min=1e-11, max=0.01, vary=True)

    # 3. Sub-nanometer clusters / high-q excess
    p.add('R_pnc', value=1.5, min=0.2, max=5.0, vary=True)
    p.add('std_pnc', value=1.0, min=0.1, max=3.0, vary=True)
    p.add('phi_pnc', value=1.4e-5, min=1e-7, max=0.001, vary=True)

    # 4. Incoherent flat background
    p.add('bkg', value=1.3e-4, min=1e-6, max=0.01, vary=True)

    def model_eval(pars, q_eval):
        I_net = calc_primary_fractal_intensity(
            q_eval, pars['Rmean'].value, pars['stdDev'].value,
            pars['phi'].value, pars['D'].value, pars['xi'].value,
            rho2=RHO2_CALIBRATED, integration_steps=80
        )
        I_por = pars['A_porod'].value / (q_eval ** 4)
        I_pnc = calc_pnc_intensity(
            q_eval, pars['R_pnc'].value, pars['std_pnc'].value,
            pars['phi_pnc'].value, rho2=RHO2_CALIBRATED
        )
        return I_net + I_por + I_pnc + pars['bkg'].value

    def res(pars):
        c = model_eval(pars, q_data)
        return (y_data - c) / dy_fit

    print("\n--- Running Multi-Component Least-Squares Fit ---")
    out = minimize(res, p, method='least_squares', max_nfev=300)
    final_params = out.params

    c_final = model_eval(final_params, q_data)
    dof = max(1, len(q_data) - 9)
    red_chi2 = np.sum(((y_data - c_final) / dy_fit) ** 2) / dof
    red_chi2_raw = np.sum(((y_data - c_final) / dy_data) ** 2) / dof

    print("\n" + "=" * 80)
    print(f"  BASELINE FIT CONVERGED: Reduced Chi-Square = {red_chi2:.4f} (raw: {red_chi2_raw:.4f})")
    print("=" * 80)

    param_dict = {}
    for p_name in final_params:
        p_val = float(final_params[p_name].value)
        p_err = float(final_params[p_name].stderr) if final_params[p_name].stderr is not None else 0.0
        param_dict[p_name] = p_val
        param_dict[f"std_err_{p_name}"] = p_err
        print(f"  {p_name:<28}: {p_val:12.6f} ± {p_err:10.6f}")

    param_dict['redchi'] = float(red_chi2)
    param_dict['redchi_raw'] = float(red_chi2_raw)

    json_path = os.path.join(out_dir, 'baseline_parameters.json')
    with open(json_path, 'w') as f:
        json.dump(param_dict, f, indent=4)
    print(f"\nSaved baseline parameters JSON: {json_path}")

    csv_path = os.path.join(out_dir, 'baseline_parameters.csv')
    pd.DataFrame([param_dict]).to_csv(csv_path, index=False)
    print(f"Saved baseline parameters CSV: {csv_path}")

    # --- Plot Publication Figure ---
    print("\nGenerating composite baseline fit figure...")
    q_dense = np.geomspace(q_data.min(), q_data.max(), 800)

    I_net_dense = calc_primary_fractal_intensity(
        q_dense, final_params['Rmean'].value, final_params['stdDev'].value,
        final_params['phi'].value, final_params['D'].value, final_params['xi'].value,
        rho2=RHO2_CALIBRATED, integration_steps=120
    )
    I_por_dense = final_params['A_porod'].value / (q_dense ** 4)
    I_pnc_dense = calc_pnc_intensity(
        q_dense, final_params['R_pnc'].value, final_params['std_pnc'].value,
        final_params['phi_pnc'].value, rho2=RHO2_CALIBRATED
    )
    I_bkg_dense = np.full_like(q_dense, final_params['bkg'].value)
    I_tot_dense = I_net_dense + I_por_dense + I_pnc_dense + I_bkg_dense

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    # Top Panel
    mask_pos = y_data > 0
    ax_fit.loglog(q_data[mask_pos], y_data[mask_pos], 'o', color='#334155', markersize=2.5, alpha=0.40,
                  label='Experimental Scattering Data')
    ax_fit.loglog(q_dense, I_tot_dense, color='#dc2626', linestyle='-', linewidth=2.2,
                  label=f'Total Unified Fit ($\\chi^2_{{\\rm red}} = {red_chi2:.2f}$)')
    ax_fit.loglog(q_dense, I_net_dense, color='#059669', linestyle='--', linewidth=1.8,
                  label=f'Spinodal Fractal Network ($R = {final_params["Rmean"].value:.1f}$ nm, $D_m = {final_params["D"].value:.2f}$, $\\xi = {final_params["xi"].value:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_por_dense, color='#d97706', linestyle='-.', linewidth=1.8,
                  label=f'Macroscopic Droplets / Porod Tail ($A = {final_params["A_porod"].value*1e8:.2f}\\times 10^{{-8}}$ cm$^{{-1}}$nm$^{{-4}}$)')
    ax_fit.loglog(q_dense, I_pnc_dense, color='#9333ea', linestyle=':', linewidth=1.8,
                  label=f'Tanaka PNC Mesh ($R_{{\\rm pnc}} = {final_params["R_pnc"].value:.1f}$ nm, $\\phi = {final_params["phi_pnc"].value:.2e}$)')
    ax_fit.axhline(final_params['bkg'].value, color='#64748b', linestyle=':', linewidth=1.2,
                   label=f'Solvent Background ({final_params["bkg"].value:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-002: Multiscale Composite Baseline Fit (20.0 mM CaCO$_3$, pH 11.3, Spinodal Gelation)\n'
                     r'Spinodal Fractal Network + Porod Tail + Tanaka PNC Mesh ($q \in [0.002, 2.54]\ \mathrm{nm}^{-1}$)',
                     fontsize=11, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.set_xlim(0.0018, 2.70)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[mask_pos].min() * 0.5), top=y_data.max() * 3.0)

    # Bottom Panel: Residuals
    residuals = (y_data - c_final) / dy_fit
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color='#475569', markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.70)
    ax_res.set_xlabel('Scattering Vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Res. [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'composite_baseline_fit.png')
    plt.savefig(plot_path, dpi=300)
    pdf_path = os.path.join(out_dir, 'composite_baseline_fit.pdf')
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved baseline fit figure: {plot_path} & {pdf_path}")
    print("=" * 80)


if __name__ == '__main__':
    fit_composite_baseline()
