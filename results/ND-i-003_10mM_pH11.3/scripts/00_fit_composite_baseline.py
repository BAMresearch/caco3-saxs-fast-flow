#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 00: Fit Composite Baseline Scattering Curve for ND-i-003 (10 mM CaCO3, pH 11.3).
---------------------------------------------------------------------------------------
Fits the full-range stitched composite curve (q in [0.002, 2.54] nm^-1):
Model: Discrete Polydisperse ACC Nanospheres (lognormal distribution) + Macroscopic Porod Interface + Flat Solvent Background.

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
from lmfit import minimize, Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_DIR)

from models.scattering_models import (
    calc_composite_model,
    calc_sphere_intensity,
    calc_large_sphere_intensity,
    RHO2_CALIBRATED
)
from models.fitting_engine import create_composite_parameters


def fit_composite_baseline():
    print("=" * 80)
    print("  COMPOSITE BASELINE SCATTERING CURVE FIT (ND-i-003, 10 mM CaCO3, pH 11.3)")
    print("  Model: Bimodal ACC Nanospheres + Macroscopic Droplets + Flat Background")
    print("=" * 80)

    out_dir = os.path.join(PROJECT_DIR, 'output', 'composite')
    os.makedirs(out_dir, exist_ok=True)

    data_path = os.path.join(PROJECT_DIR, 'data', 'dataNDi003_fitting_results_composite.csv')
    df = pd.read_csv(data_path)
    q = df['q'].values
    y = df['0'].values
    dy = np.maximum(df['error_0'].values, 1e-4)

    # Fit across full experimental range
    mask = (q >= 0.002) & (q <= 2.54) & (y > -0.01)
    q_fit = q[mask]
    y_fit = y[mask]
    dy_fit = np.maximum(dy[mask], 0.03 * np.maximum(y_fit, 0))

    print(f"Loaded composite dataset ({len(df)} points). Fitting across {len(q_fit)} points (q: [{q_fit.min():.5f}, {q_fit.max():.5f}] nm^-1)")

    params = create_composite_parameters()

    def residual(p, q_pts, y_pts, dy_pts):
        ymod = calc_composite_model(p, q_pts)
        return (y_pts - ymod) / dy_pts

    print("\nRunning least_squares optimization...")
    out = minimize(residual, params, args=(q_fit, y_fit, dy_fit), method='least_squares', max_nfev=600)

    # Estimate covariance matrix and standard errors
    p_cov = out.params.copy()
    for k in p_cov:
        if k not in ['rho2']:
            p_cov[k].vary = True
    try:
        out_cov = minimize(residual, p_cov, args=(q_fit, y_fit, dy_fit), method='leastsq', calc_covar=True)
        final_params = out_cov.params
    except Exception as e:
        print(f"Covariance warning: {e}, falling back to least_squares parameters.")
        final_params = out.params

    c_final = calc_composite_model(final_params, q_fit)
    dof = max(1, len(q_fit) - 10)
    red_chi2 = np.sum(((y_fit - c_final) / dy_fit) ** 2) / dof
    red_chi2_raw = np.sum(((y_fit - c_final) / dy[mask]) ** 2) / dof

    print("\nFit Converged!")
    for k in final_params:
        err = final_params[k].stderr if final_params[k].stderr is not None else 0.0
        print(f"  {k:28s} = {final_params[k].value:10.4e} +/- {err:10.4e}")
    print(f"  Reduced chi2 = {red_chi2:.3f} (raw: {red_chi2_raw:.3f})")

    base_dict = {
        'Rmean': float(final_params['Rmean'].value),
        'std_err_Rmean': float(final_params['Rmean'].stderr or 0.0),
        'stdDev': float(final_params['stdDev'].value),
        'std_err_stdDev': float(final_params['stdDev'].stderr or 0.0),
        'phi': float(final_params['phi'].value),
        'std_err_phi': float(final_params['phi'].stderr or 0.0),
        'R2': float(final_params['R2'].value),
        'std_err_R2': float(final_params['R2'].stderr or 0.0),
        'stdDev2': float(final_params['stdDev2'].value),
        'std_err_stdDev2': float(final_params['stdDev2'].stderr or 0.0),
        'phi2': float(final_params['phi2'].value),
        'std_err_phi2': float(final_params['phi2'].stderr or 0.0),
        'mean_radius_large': float(final_params['mean_radius_large'].value),
        'std_err_mean_radius_large': float(final_params['mean_radius_large'].stderr or 0.0),
        'std_radius_large': float(final_params['std_radius_large'].value),
        'std_err_std_radius_large': float(final_params['std_radius_large'].stderr or 0.0),
        'total_volume_fraction_large': float(final_params['total_volume_fraction_large'].value),
        'std_err_total_volume_fraction_large': float(final_params['total_volume_fraction_large'].stderr or 0.0),
        'A_porod': 0.0,
        'std_err_A_porod': 0.0,
        'bkg': float(final_params['bkg'].value),
        'std_err_bkg': float(final_params['bkg'].stderr or 0.0),
        'chisqr': float(out.chisqr),
        'redchi': float(red_chi2),
        'redchi_raw': float(red_chi2_raw)
    }

    with open(os.path.join(out_dir, 'baseline_parameters.json'), 'w') as f:
        json.dump(base_dict, f, indent=4)

    pd.DataFrame([base_dict]).to_csv(os.path.join(out_dir, 'baseline_parameters.csv'), index=False)

    # Plot publication figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
    fig.suptitle('ND-i-003: Full-Range Composite Baseline Fit (10 mM $\\mathrm{CaCO}_3$, pH 11.3)\n'
                 'Bimodal ACC Nanospheres + Macroscopic Droplets', fontsize=12, fontweight='bold')

    q_dense = np.geomspace(q.min(), q.max(), 1000)
    y_dense = calc_composite_model(final_params, q_dense)
    y_acc1_dense = calc_sphere_intensity(q_dense, base_dict['Rmean'], base_dict['stdDev'], base_dict['phi'], bkg=0.0, rho2=1.0, dist_type='lognormal')
    y_acc2_dense = calc_sphere_intensity(q_dense, base_dict['R2'], base_dict['stdDev2'], base_dict['phi2'], bkg=0.0, rho2=1.0, dist_type='lognormal')
    y_drop_dense = calc_large_sphere_intensity(q_dense, base_dict['mean_radius_large'], base_dict['std_radius_large'], base_dict['total_volume_fraction_large'], rho2=1.0)

    mask_pos = y > 0
    ax1.errorbar(q[mask_pos], y[mask_pos], yerr=dy[mask_pos], fmt='o', markersize=2.5,
                 color='#334155', alpha=0.45, label='Experimental Scattering Data')
    ax1.plot(q_dense, y_dense, color='#dc2626', linewidth=2.2,
             label=f'Total Unified Fit ($\\chi^2_{{\\rm red}} = {base_dict["redchi"]:.2f}$)')
    ax1.plot(q_dense, y_acc1_dense, color='#0284c7', linestyle='-.', linewidth=1.8,
             label=f'Discrete ACC Spheres (Mode 1) ($\\langle R_1 \\rangle = {base_dict["Rmean"]:.1f}$ nm, $\\sigma_1 = {base_dict["stdDev"]:.1f}$ nm, $\\phi_1 = {base_dict["phi"]:.2e}$)')
    ax1.plot(q_dense, y_acc2_dense, color='#0ea5e9', linestyle=':', linewidth=1.6,
             label=f'Discrete ACC Spheres (Mode 2) ($\\langle R_2 \\rangle = {base_dict["R2"]:.1f}$ nm, $\\sigma_2 = {base_dict["stdDev2"]:.1f}$ nm, $\\phi_2 = {base_dict["phi2"]:.2e}$)')
    ax1.plot(q_dense, y_drop_dense, color='#d97706', linestyle='-.', linewidth=1.8,
             label=f'Macroscopic Droplets / Porod Tail ($\\langle R_{{\\rm large}} \\rangle = {base_dict["mean_radius_large"]:.0f}$ nm, $\\phi_{{\\rm large}} = {base_dict["total_volume_fraction_large"]:.2e}$)')
    ax1.axhline(base_dict['bkg'], color='#64748b', linestyle=':', linewidth=1.2,
                label=f'Solvent Background ({base_dict["bkg"]:.2e} cm$^{{-1}}$)')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(q.min() * 0.9, q.max() * 1.1)
    ax1.set_ylim(1e-4, 5e3)
    ax1.set_ylabel(r'Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax1.legend(loc='upper right', fontsize=9.0, framealpha=0.92)

    # Residuals
    residuals = (y_fit - c_final) / dy_fit
    ax2.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax2.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax2.plot(q_fit, residuals, 'o', markersize=2, color='#475569', alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_xlim(q.min() * 0.9, q.max() * 1.1)
    ax2.set_ylim(-5, 5)
    ax2.set_xlabel(r'Scattering Vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax2.set_ylabel(r'Res. [$\sigma$]', fontsize=10, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'composite_baseline_fit.png')
    pdf_path = os.path.join(out_dir, 'composite_baseline_fit.pdf')
    plt.savefig(plot_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close(fig)

    print(f"Saved baseline outputs to: {out_dir}")
    return base_dict


if __name__ == '__main__':
    fit_composite_baseline()
