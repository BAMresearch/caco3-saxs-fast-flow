#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Standardized Composite Baseline Plot Generator
----------------------------------------------
Generates publication-quality baseline fit figures and component decompositions
for all 6 experimental reaction series using the standardized project palette.

Fixed, unique component colors:
- Data: #334155 (slate dark charcoal markers)
- Total unified fit: #dc2626 (crimson red solid line)
- Spinodal fractal network: #059669 (emerald green dashed line)
- Discrete ACC nanospheres: #0284c7 (ocean blue dash-dot line)
- Secondary ACC mode (10 mM): #0ea5e9 (sky blue dotted line)
- Giant droplets / Porod tail: #d97706 (warm amber dash-dot line)
- Viscoelastic mesh / PNC fluctuations: #9333ea (violet purple dotted line)
- Incoherent background: #64748b (slate grey dotted line)
- Residuals: #475569 (slate grey markers)
"""

import os
import sys
import json
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, RESULTS_DIR)
from saxs_colors import SAXS_COLORS, SAXS_STYLES


def load_scattering_models(sdir, module_alias):
    mod_path = os.path.join(sdir, 'models', 'scattering_models.py')
    spec = importlib.util.spec_from_file_location(module_alias, mod_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_alias] = mod
    spec.loader.exec_module(mod)
    return mod


def load_data(sdir, csv_name):
    path = os.path.join(sdir, 'data', csv_name)
    df = pd.read_csv(path)
    q = df['q'].values
    if '0' in df.columns:
        y = df['0'].values
        dy = np.maximum(df['error_0'].values, 1e-4)
    else:
        y = df['I'].values
        dy = np.maximum(df['dI'].values, 1e-4)
    return q, y, dy


def plot_nd_i_003():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-003_10mM_pH11.3')
    m03 = load_scattering_models(sdir, 'm03_models')
    from lmfit import Parameters

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi003_fitting_results_composite.csv')
    mask = (q >= 0.002) & (q <= 2.54) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    q_dense = np.geomspace(q_data.min(), q_data.max(), 1000)

    fit_params = Parameters()
    for k in ['Rmean', 'stdDev', 'phi', 'R2', 'stdDev2', 'phi2',
              'mean_radius_large', 'std_radius_large', 'total_volume_fraction_large', 'bkg']:
        fit_params.add(k, value=p[k])
    fit_params.add('rho2', value=1.0)

    y_dense = m03.calc_composite_model(fit_params, q_dense)
    y_acc1_dense = m03.calc_sphere_intensity(q_dense, p['Rmean'], p['stdDev'], p['phi'], bkg=0.0, rho2=1.0, dist_type='lognormal')
    y_acc2_dense = m03.calc_sphere_intensity(q_dense, p['R2'], p['stdDev2'], p['phi2'], bkg=0.0, rho2=1.0, dist_type='lognormal')
    y_drop_dense = m03.calc_large_sphere_intensity(q_dense, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'], rho2=1.0)
    c_final = m03.calc_composite_model(fit_params, q_data)

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    # Top panel
    ax_fit.errorbar(q_data, y_data, yerr=dy_data, fmt='o', markersize=2.5,
                    color=SAXS_COLORS['data'], alpha=0.45, label='Composite scattering data (ND-i-003, 10 mM)')
    ax_fit.plot(q_dense, y_dense, color=SAXS_COLORS['total_fit'], linewidth=2.2,
                label=f'Total unified fit ($\\chi^2_{{\\rm red}} = {p["redchi"]:.2f}$)')
    ax_fit.plot(q_dense, y_acc1_dense, color=SAXS_COLORS['acc_spheres'], linestyle='-.', linewidth=1.8,
                label=f'Discrete ACC spheres (mode 1) ($\\langle R_1 \\rangle = {p["Rmean"]:.1f}$ nm, $\\sigma_1 = {p["stdDev"]:.1f}$ nm)')
    ax_fit.plot(q_dense, y_acc2_dense, color=SAXS_COLORS['acc_spheres_mode2'], linestyle=':', linewidth=1.6,
                label=f'Discrete ACC spheres (mode 2) ($\\langle R_2 \\rangle = {p["R2"]:.1f}$ nm, $\\sigma_2 = {p["stdDev2"]:.1f}$ nm)')
    ax_fit.plot(q_dense, y_drop_dense, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.8,
                label=f'Macroscopic droplets / Porod tail ($\\langle R_{{\\rm lg}} \\rangle = {p["mean_radius_large"]:.0f}$ nm, $\\phi_{{\\rm lg}} = {p["total_volume_fraction_large"]:.2e}$)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_xscale('log')
    ax_fit.set_yscale('log')
    ax_fit.set_xlim(0.0018, 2.70)
    ax_fit.set_ylim(1e-4, 5e3)
    ax_fit.set_ylabel(r'Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-003: Multi-scale composite baseline fit (10.0 mM CaCO$_3$, pH 11.3, sub-spinodal)\n'
                     r'Bimodal ACC nanospheres + dense liquid droplets ($q \in [0.002, 2.54]\ \mathrm{nm}^{-1}$)',
                     fontsize=11, fontweight='bold')
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)

    # Residuals
    residuals = (y_data - c_final) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', markersize=2.2, color=SAXS_COLORS['residuals'], alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.70)
    ax_res.set_ylim(-5, 5)
    ax_res.set_xlabel(r'Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel(r'Residuals [$\sigma$]', fontsize=10, fontweight='bold')
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    merged_dir = os.path.join(sdir, 'output', 'merged')
    if os.path.exists(merged_dir):
        plt.savefig(os.path.join(merged_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(merged_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.close()
    print("ND-i-003 standardized baseline figure generated.")


def plot_nd_i_004():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-004_15mM_pH11.3')
    m04 = load_scattering_models(sdir, 'm04_models')

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi004_fitting_results_composite.csv')
    mask = (q >= 0.002) & (q <= 2.50) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    q_dense = np.geomspace(q_data.min(), q_data.max(), 600)

    I_net = m04.calc_primary_fractal_intensity(
        q_dense, p['R_prim'], p['std_prim'], p['phi_net'], p['Dm'], p['xi'],
        rho2=m04.RHO2_CALIBRATED, integration_steps=120
    )
    I_acc = m04.calc_large_sphere_intensity(
        q_dense, p['R_acc'], p['std_acc'], p['phi_acc'],
        rho2=m04.RHO2_CALIBRATED, integration_steps=120
    )
    I_por = m04.calc_porod_intensity(q_dense, p['A_porod'], 4.0)
    I_bkg = np.full_like(q_dense, p['bkg'])
    I_tot = I_net + I_acc + I_por + I_bkg

    I_net_d = m04.calc_primary_fractal_intensity(q_data, p['R_prim'], p['std_prim'], p['phi_net'], p['Dm'], p['xi'], rho2=m04.RHO2_CALIBRATED)
    I_acc_d = m04.calc_large_sphere_intensity(q_data, p['R_acc'], p['std_acc'], p['phi_acc'], rho2=m04.RHO2_CALIBRATED)
    I_por_d = m04.calc_porod_intensity(q_data, p['A_porod'], 4.0)
    c_final = I_net_d + I_acc_d + I_por_d + p['bkg']

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color=SAXS_COLORS['data'], markersize=2.5, alpha=0.40, label='Composite scattering data (ND-i-004, 15 mM)')
    ax_fit.loglog(q_dense, I_tot, color=SAXS_COLORS['total_fit'], linestyle='-', linewidth=2.2, label=f'Total unified fit ($\\chi^2_{{\\mathrm{{red}}}} = {p["redchi"]:.2f}$)')
    ax_fit.loglog(q_dense, I_net, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.8,
                  label=f'Spinodal fractal network ($R_{{\\rm prim}} = {p["R_prim"]:.1f}$ nm, $D_m = {p["Dm"]:.2f}$, $\\xi = {p["xi"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_acc, color=SAXS_COLORS['acc_spheres'], linestyle='-.', linewidth=1.8,
                  label=f'Discrete ACC spheres ($\\langle R \\rangle = {p["R_acc"]:.1f} \\pm {p["std_acc"]:.1f}$ nm, $\\phi = {p["phi_acc"]:.2e}$)')
    ax_fit.loglog(q_dense, I_por, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.6,
                  label=f'Macroscopic droplets / Porod tail ($A = {p["A_porod"]*1e6:.2f}\\times 10^{{-6}}$)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-004: Multi-scale composite baseline fit (15.0 mM CaCO$_3$, pH 11.3, crossover boundary)\n'
                     r'Spinodal mass fractal network + ACC nanosphere form factor ($q \in [0.002, 2.5]\ \mathrm{nm}^{-1}$)', fontsize=11, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.set_xlim(0.0018, 2.65)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[y_data > 0].min() * 0.5), top=y_data.max() * 3.0)

    residuals = (y_data - c_final) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color=SAXS_COLORS['residuals'], markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.65)
    ax_res.set_xlabel('Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Residuals [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.close()
    print("ND-i-004 standardized baseline figure generated.")


def plot_nd_i_002():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-002_20mM_pH11.3')
    m02 = load_scattering_models(sdir, 'm02_models')

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi002_fitting_results_composite.csv')
    mask = (q >= 0.002) & (q <= 2.54) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    q_dense = np.geomspace(q_data.min(), q_data.max(), 600)

    I_net_dense = m02.calc_primary_fractal_intensity(
        q_dense, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'],
        rho2=m02.RHO2_CALIBRATED, integration_steps=150
    )
    I_por_dense = m02.calc_porod_intensity(q_dense, p['A_porod'], 4.0)
    I_fluc_dense = m02.calc_fb_intensity(
        q_dense, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0)
    )
    I_bkg_dense = np.full_like(q_dense, p['bkg'])
    I_tot_dense = I_net_dense + I_por_dense + I_fluc_dense + I_bkg_dense

    I_net_d = m02.calc_primary_fractal_intensity(q_data, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'], rho2=m02.RHO2_CALIBRATED)
    I_por_d = m02.calc_porod_intensity(q_data, p['A_porod'], 4.0)
    I_fluc_d = m02.calc_fb_intensity(q_data, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0))
    c_final = I_net_d + I_por_d + I_fluc_d + p['bkg']

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color=SAXS_COLORS['data'], markersize=2.5, alpha=0.40,
                  label='Composite scattering data (ND-i-002, 20 mM)')
    ax_fit.loglog(q_dense, I_tot_dense, color=SAXS_COLORS['total_fit'], linestyle='-', linewidth=2.2,
                  label=f'Total unified fit ($\\chi^2_{{\\rm red}} = {p["redchi"]:.2f}$)')
    ax_fit.loglog(q_dense, I_net_dense, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.8,
                  label=f'Spinodal fractal network ($R = {p["Rmean"]:.1f}$ nm, $D_m = {p["D"]:.2f}$, $\\xi = {p["xi"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_por_dense, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.8,
                  label=f'Macroscopic droplets / Porod tail ($A = {p["A_porod"]*1e8:.2f}\\times 10^{{-8}}$ cm$^{{-1}}$nm$^{{-4}}$)')
    ax_fit.loglog(q_dense, I_fluc_dense, color=SAXS_COLORS['pnc_fluctuations'], linestyle=':', linewidth=1.8,
                  label=f'Tanaka viscoelastic fluctuations ($\\xi_{{\\rm fluc}} = {p["xi_fluc"]:.1f}$ nm, $A_{{\\rm fluc}} = {p["A_fluc"]:.3f}$ cm$^{{-1}}$)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-002: Multi-scale composite baseline fit (20.0 mM CaCO$_3$, pH 11.3, spinodal gelation)\n'
                     r'Spinodal fractal network + Porod tail + sub-nm clusters ($q \in [0.002, 2.54]\ \mathrm{nm}^{-1}$)',
                     fontsize=11, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.set_xlim(0.0018, 2.70)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[y_data > 0].min() * 0.5), top=y_data.max() * 3.0)

    residuals = (y_data - c_final) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color=SAXS_COLORS['residuals'], markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.70)
    ax_res.set_xlabel('Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Residuals [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.close()
    print("ND-i-002 standardized baseline figure generated.")


def plot_nd_i_009():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-009_40mM_pH10.0')
    m09 = load_scattering_models(sdir, 'm09_models')
    from lmfit import Parameters

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi009_fitting_results_merged.csv')
    mask = (q >= 0.002) & (q <= 2.50) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    params = Parameters()
    for k, v in p.items():
        if not k.startswith('std_err_'):
            params.add(k, value=v)

    q_dense = np.geomspace(q_data.min(), q_data.max(), 600)
    I_tot, I_prim, I_large, I_oz, I_bkg = m09.decompose_model_components(params, q_dense)
    y_calc_data = m09.calc_model_intensity(params, q_data)

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color=SAXS_COLORS['data'], markersize=2.5, alpha=0.40, label='Composite scattering data (ND-i-009, 40 mM, pH 10.0)')
    ax_fit.loglog(q_dense, I_tot, color=SAXS_COLORS['total_fit'], linestyle='-', linewidth=2.2, label=f'Total unified fit ($\\chi^2_{{\\mathrm{{red}}}} = {p["redchi"]:.2f}$)')
    ax_fit.loglog(q_dense, I_prim, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.8,
                  label=f'Spinodal fractal network ($R = {p["Rmean"]:.2f}$ nm, $D_m = {p["D"]:.2f}$, $\\xi = {p["xi"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_large, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.8,
                  label=f'Macroscopic droplets / Porod tail ($\\langle R \\rangle = {p["mean_radius_large"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_oz, color=SAXS_COLORS['pnc_fluctuations'], linestyle=':', linewidth=1.8,
                  label=f'Tanaka PNC fluctuation mesh ($\\xi_{{\\rm fluc}} = {p["xi_OZ"]:.2f}$ nm)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-009: Multi-scale composite baseline fit (40.0 mM CaCO$_3$, pH 10.0, pure benchmark)\n'
                     r'Tanaka VPS network + droplet inversion + viscoelastic mesh ($q \in [0.002, 2.5]\ \mathrm{nm}^{-1}$)',
                     fontsize=11, fontweight='bold')
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.legend(loc='upper right', frameon=True, fontsize=8.8, framealpha=0.92)
    ax_fit.set_xlim(0.0018, 2.65)
    ax_fit.set_ylim(bottom=1e-1, top=1e7)

    residuals = (y_data - y_calc_data) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color=SAXS_COLORS['residuals'], markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.65)
    ax_res.set_xlabel('Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Residuals [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    merged_dir = os.path.join(sdir, 'output', 'merged')
    if os.path.exists(merged_dir):
        plt.savefig(os.path.join(merged_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(merged_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(merged_dir, 'merged_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("ND-i-009 standardized baseline figure generated.")


def plot_nd_i_001():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-001_40mM_pH11.3')
    m01 = load_scattering_models(sdir, 'm01_models')

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi001_fitting_results_composite.csv')
    mask = (q >= 0.002) & (q <= 2.50) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    q_dense = np.geomspace(q_data.min(), q_data.max(), 800)

    I_prim = m01.calc_primary_fractal_intensity(
        q_dense, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'],
        rho2=m01.RHO2_CALIBRATED, integration_steps=200
    )
    I_large = m01.calc_large_sphere_intensity(
        q_dense, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'],
        rho2=m01.RHO2_CALIBRATED, integration_steps=200
    )
    I_fluc = m01.calc_fb_intensity(
        q_dense, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0)
    )
    I_bkg = np.full_like(q_dense, p['bkg'])
    I_tot = I_prim + I_large + I_fluc + I_bkg

    I_prim_d = m01.calc_primary_fractal_intensity(q_data, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'], rho2=m01.RHO2_CALIBRATED)
    I_large_d = m01.calc_large_sphere_intensity(q_data, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'], rho2=m01.RHO2_CALIBRATED)
    I_fluc_d = m01.calc_fb_intensity(q_data, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0))
    c_final = I_prim_d + I_large_d + I_fluc_d + p['bkg']

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color=SAXS_COLORS['data'], markersize=2.5, alpha=0.45, label='Composite scattering data (ND-i-001, 40 mM, pH 11.3)')
    ax_fit.loglog(q_dense, I_tot, color=SAXS_COLORS['total_fit'], linestyle='-', linewidth=2.2, label=f'Total unified fit ($\\chi^2_{{\\mathrm{{red}}}} = {p["redchi"]:.2f}$)')
    ax_fit.loglog(q_dense, I_prim, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.8,
                  label=f'Spinodal fractal network ($\\langle R \\rangle = {p["Rmean"]:.2f}$ nm, $D_f = {p["D"]:.2f}$)')
    ax_fit.loglog(q_dense, I_large, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.8,
                  label=f'Macroscopic droplets / Porod tail ($\\langle R \\rangle = {p["mean_radius_large"]:.0f}$ nm, $\\phi = {p["total_volume_fraction_large"]:.2e}$)')
    ax_fit.loglog(q_dense, I_fluc, color=SAXS_COLORS['pnc_fluctuations'], linestyle=':', linewidth=1.8,
                  label=f'Tanaka viscoelastic fluctuations ($\\xi_{{\\rm fluc}} = {p["xi_fluc"]:.1f}$ nm, $A_{{\\rm fluc}} = {p["A_fluc"]:.3f}$ cm$^{{-1}}$)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-001: Full-range multi-scale composite baseline fit (40.0 mM CaCO$_3$, pH 11.3, pure instant quench)\n'
                     r'$q \in [0.002, 2.5]\ \mathrm{nm}^{-1}$ across 4 length scale decades', fontsize=11, fontweight='bold')
    ax_fit.legend(loc='upper right', fontsize=8.8, framealpha=0.92)
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.set_xlim(0.0018, 2.65)
    ax_fit.set_ylim(bottom=max(1e-4, y_data[y_data > 0].min() * 0.5), top=y_data.max() * 3.0)

    residuals = (y_data - c_final) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color=SAXS_COLORS['residuals'], markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.65)
    ax_res.set_xlabel('Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Residuals [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.close()
    print("ND-i-001 standardized baseline figure generated.")


def plot_nd_i_005():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-005_40mM_pH11.3_CoCl2')
    m05 = load_scattering_models(sdir, 'm05_models')

    json_path = os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')
    out_dir = os.path.join(sdir, 'output', 'composite')

    with open(json_path) as f:
        p = json.load(f)

    q, y, dy = load_data(sdir, 'dataNDi005_fitting_results_composite.csv')
    mask = (q >= 0.002) & (q <= 2.50) & (y > 0)
    q_data, y_data, dy_data = q[mask], y[mask], dy[mask]

    q_dense = np.geomspace(q_data.min(), q_data.max(), 800)

    I_prim = m05.calc_primary_fractal_intensity(
        q_dense, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'],
        rho2=m05.RHO2_CALIBRATED, integration_steps=200
    )
    I_large = m05.calc_large_sphere_intensity(
        q_dense, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'],
        rho2=m05.RHO2_CALIBRATED, integration_steps=200
    )
    I_fluc = m05.calc_fb_intensity(
        q_dense, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0)
    )
    I_bkg = np.full_like(q_dense, p['bkg'])
    I_tot = I_prim + I_large + I_fluc + I_bkg

    I_prim_d = m05.calc_primary_fractal_intensity(q_data, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'], rho2=m05.RHO2_CALIBRATED)
    I_large_d = m05.calc_large_sphere_intensity(q_data, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'], rho2=m05.RHO2_CALIBRATED)
    I_fluc_d = m05.calc_fb_intensity(q_data, p['A_fluc'], p['xi_fluc'], mu_fluc=p.get('mu_fluc', 4.0))
    c_final = I_prim_d + I_large_d + I_fluc_d + p['bkg']

    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.10)
    ax_fit = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_fit)

    ax_fit.loglog(q_data, y_data, 'o', color=SAXS_COLORS['data'], markersize=2.5, alpha=0.45, label='Composite scattering data (ND-i-005, 40 mM + CoCl$_2$)')
    ax_fit.loglog(q_dense, I_tot, color=SAXS_COLORS['total_fit'], linestyle='-', linewidth=2.2, label=f'Total unified fit ($\\chi^2_{{\\mathrm{{red}}}} = {p["redchi"]:.2f}$)')
    ax_fit.loglog(q_dense, I_prim, color=SAXS_COLORS['fractal_network'], linestyle='--', linewidth=1.8,
                  label=f'Spinodal fractal network ($\\langle R \\rangle = {p["Rmean"]:.2f}$ nm, $D_f = {p["D"]:.2f}$, $\\xi = {p["xi"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_large, color=SAXS_COLORS['droplets_porod'], linestyle='-.', linewidth=1.8,
                  label=f'Macroscopic droplets / Porod tail ($\\langle R \\rangle = {p["mean_radius_large"]:.0f}$ nm)')
    ax_fit.loglog(q_dense, I_fluc, color=SAXS_COLORS['pnc_fluctuations'], linestyle=':', linewidth=1.8,
                  label=f'Tanaka viscoelastic fluctuations ($\\xi_{{\\rm fluc}} = {p["xi_fluc"]:.1f}$ nm, $A_{{\\rm fluc}} = {p["A_fluc"]:.3f}$ cm$^{{-1}}$)')
    ax_fit.axhline(p['bkg'], color=SAXS_COLORS['background'], linestyle=':', linewidth=1.2,
                   label=f'Solvent background ({p["bkg"]:.2e} cm$^{{-1}}$)')

    ax_fit.set_ylabel('Intensity $I(q)$ [cm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_fit.set_title('ND-i-005: Multi-scale baseline scattering curve fit and component decomposition\n(40.0 mM CaCO$_3$, pH 11.3 with cobalt additive perturbation)', fontsize=11, fontweight='bold')
    ax_fit.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)
    ax_fit.legend(loc='upper right', frameon=True, fontsize=8.8, framealpha=0.92)
    ax_fit.set_xlim(0.0018, 2.65)
    ax_fit.set_ylim(bottom=1e-3, top=5e4)

    residuals = (y_data - c_final) / dy_data
    ax_res.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax_res.axhspan(-2, 2, color='#e2e8f0', alpha=0.5)
    ax_res.plot(q_data, residuals, 'o', color=SAXS_COLORS['residuals'], markersize=2.2, alpha=0.5)
    ax_res.set_xscale('log')
    ax_res.set_xlim(0.0018, 2.65)
    ax_res.set_xlabel('Scattering vector $q$ [nm$^{-1}$]', fontsize=12, fontweight='bold')
    ax_res.set_ylabel('Residuals [$\\sigma$]', fontsize=10, fontweight='bold')
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'composite_baseline_fit.pdf'), bbox_inches='tight')
    plt.close()
    print("ND-i-005 standardized baseline figure generated.")


def main():
    print("Standardizing composite baseline plots across all 6 reaction series...")
    plot_nd_i_003()
    plot_nd_i_004()
    plot_nd_i_002()
    plot_nd_i_009()
    plot_nd_i_001()
    plot_nd_i_005()
    print("All 6 composite baseline figures successfully standardized and saved.")


if __name__ == '__main__':
    main()
