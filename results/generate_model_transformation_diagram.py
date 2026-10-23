#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Model Transformation Diagram Generator (Vector SVG)
---------------------------------------------------
Generates a vector SVG diagram built directly around actual multi-scale
composite SAXS/WAXS datasets (q in [0.002, 2.5] nm^-1) and converged physical models.

Standardized project palette:
- Raw Data Points: #334155 (slate dark charcoal)
- Total Unified Model: #dc2626 (crimson red)
- Spinodal Fractal Network: #059669 (emerald green)
- Discrete ACC Nanospheres: #0284c7 (ocean blue)
- Secondary ACC Mode (10 mM): #0ea5e9 (sky blue)
- Giant Droplets / Porod Tail: #d97706 (warm amber)
- Incoherent Background: #64748b (slate grey)

Output: results/output/model_transformation_diagram.svg
"""

import os
import sys
import json
import importlib.util
import numpy as np
import pandas as pd

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, RESULTS_DIR)
from saxs_colors import SAXS_COLORS


def load_model_module(sdir, alias):
    mod_path = os.path.join(sdir, 'models', 'scattering_models.py')
    spec = importlib.util.spec_from_file_location(alias, mod_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def load_series_data(sdir, csv_name):
    path = os.path.join(sdir, 'data', csv_name)
    df = pd.read_csv(path)
    q = df['q'].values
    if '0' in df.columns:
        y = df['0'].values
    else:
        y = df['I'].values
    mask = (q >= 0.002) & (q <= 2.50) & (y > 0)
    return q[mask], y[mask]


def compute_svg_panel(q_data, y_data, model_curves, width=216, height=92, q_range=(0.002, 2.5), I_range=(1e-4, 2e5)):
    """
    Renders authentic log-log coordinate grid, experimental data points,
    and standardized physical model decomposition paths into SVG markup.
    """
    log_q_min, log_q_max = np.log10(q_range[0]), np.log10(q_range[1])
    log_I_min, log_I_max = np.log10(I_range[0]), np.log10(I_range[1])

    def to_svg(q_val, I_val):
        x = width * (np.log10(q_val) - log_q_min) / (log_q_max - log_q_min)
        y = height * (1.0 - (np.log10(I_val) - log_I_min) / (log_I_max - log_I_min))
        return np.clip(x, 0, width), np.clip(y, 0, height)

    elements = []

    # Container frame background
    elements.append(f'<rect width="{width}" height="{height}" rx="6" fill="#070d19" stroke="#1e293b" stroke-width="1"/>')

    # Grid lines at logarithmic decades
    for q_dec, label in [(0.01, '10⁻²'), (0.1, '10⁻¹'), (1.0, '10⁰')]:
        gx = width * (np.log10(q_dec) - log_q_min) / (log_q_max - log_q_min)
        elements.append(f'<line x1="{gx:.1f}" y1="8" x2="{gx:.1f}" y2="{height-8}" stroke="#1e293b" stroke-width="0.8" stroke-dasharray="2,2"/>')
        elements.append(f'<text x="{gx:.1f}" y="{height-2}" text-anchor="middle" fill="#475569" font-size="7">{label}</text>')

    for I_dec in [1e-2, 1e0, 1e2, 1e4]:
        gy = height * (1.0 - (np.log10(I_dec) - log_I_min) / (log_I_max - log_I_min))
        if 8 <= gy <= height - 8:
            elements.append(f'<line x1="10" y1="{gy:.1f}" x2="{width-10}" y2="{gy:.1f}" stroke="#1e293b" stroke-width="0.8" stroke-dasharray="2,2"/>')

    # Experimental data points (subsampled for clean rendering)
    step = max(1, len(q_data) // 35)
    for q_pt, y_pt in zip(q_data[::step], y_data[::step]):
        if y_pt > 0:
            sx, sy = to_svg(q_pt, y_pt)
            elements.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="1.3" fill="#64748b" opacity="0.6"/>')

    # Model curves with standardized colors & dash patterns
    for curve in model_curves:
        qc, yc, color, width_val, dash = curve['q'], curve['y'], curve['color'], curve['width'], curve.get('dash', None)
        coords = [to_svg(q_pt, y_pt) for q_pt, y_pt in zip(qc, yc) if y_pt > 0]
        if len(coords) > 1:
            d_str = 'M ' + ' L '.join([f'{x:.1f},{y:.1f}' for x, y in coords])
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ''
            elements.append(f'<path d="{d_str}" fill="none" stroke="{color}" stroke-width="{width_val}"{dash_attr}/>')

    return '\n      '.join(elements)


def build_curves_nd_i_003():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-003_10mM_pH11.3')
    m03 = load_model_module(sdir, 'svg_m03')
    with open(os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')) as f:
        p = json.load(f)
    q_data, y_data = load_series_data(sdir, 'dataNDi003_fitting_results_composite.csv')

    q_dense = np.geomspace(0.002, 2.50, 200)
    from lmfit import Parameters
    fit_params = Parameters()
    for k in ['Rmean', 'stdDev', 'phi', 'R2', 'stdDev2', 'phi2',
              'mean_radius_large', 'std_radius_large', 'total_volume_fraction_large', 'bkg']:
        fit_params.add(k, value=p[k])
    fit_params.add('rho2', value=1.0)

    y_tot = m03.calc_composite_model(fit_params, q_dense)
    y_acc1 = m03.calc_sphere_intensity(q_dense, p['Rmean'], p['stdDev'], p['phi'], bkg=0.0, rho2=1.0)
    y_drop = m03.calc_large_sphere_intensity(q_dense, p['mean_radius_large'], p['std_radius_large'], p['total_volume_fraction_large'], rho2=1.0)
    y_bkg = np.full_like(q_dense, p['bkg'])

    curves = [
        {'q': q_dense, 'y': y_drop, 'color': SAXS_COLORS['droplets_porod'], 'width': 1.5, 'dash': '3,2'},
        {'q': q_dense, 'y': y_acc1, 'color': SAXS_COLORS['acc_spheres'], 'width': 1.5, 'dash': '5,2,1,2'},
        {'q': q_dense, 'y': y_bkg, 'color': SAXS_COLORS['background'], 'width': 1.0, 'dash': '1,2'},
        {'q': q_dense, 'y': y_tot, 'color': SAXS_COLORS['total_fit'], 'width': 2.0},
    ]
    return q_data, y_data, curves


def build_curves_nd_i_004():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-004_15mM_pH11.3')
    m04 = load_model_module(sdir, 'svg_m04')
    with open(os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')) as f:
        p = json.load(f)
    q_data, y_data = load_series_data(sdir, 'dataNDi004_fitting_results_composite.csv')

    q_dense = np.geomspace(0.002, 2.50, 200)
    I_net = m04.calc_primary_fractal_intensity(q_dense, p['R_prim'], p['std_prim'], p['phi_net'], p['Dm'], p['xi'], rho2=m04.RHO2_CALIBRATED, integration_steps=80)
    I_acc = m04.calc_large_sphere_intensity(q_dense, p['R_acc'], p['std_acc'], p['phi_acc'], rho2=m04.RHO2_CALIBRATED, integration_steps=80)
    I_por = m04.calc_porod_intensity(q_dense, p['A_porod'], 4.0)
    I_bkg = np.full_like(q_dense, p['bkg'])
    I_tot = I_net + I_acc + I_por + I_bkg

    curves = [
        {'q': q_dense, 'y': I_por, 'color': SAXS_COLORS['droplets_porod'], 'width': 1.5, 'dash': '3,2'},
        {'q': q_dense, 'y': I_acc, 'color': SAXS_COLORS['acc_spheres'], 'width': 1.5, 'dash': '5,2,1,2'},
        {'q': q_dense, 'y': I_net, 'color': SAXS_COLORS['fractal_network'], 'width': 1.5, 'dash': '4,2'},
        {'q': q_dense, 'y': I_bkg, 'color': SAXS_COLORS['background'], 'width': 1.0, 'dash': '1,2'},
        {'q': q_dense, 'y': I_tot, 'color': SAXS_COLORS['total_fit'], 'width': 2.0},
    ]
    return q_data, y_data, curves


def build_curves_nd_i_002():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-002_20mM_pH11.3')
    m02 = load_model_module(sdir, 'svg_m02')
    with open(os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')) as f:
        p = json.load(f)
    q_data, y_data = load_series_data(sdir, 'dataNDi002_fitting_results_composite.csv')

    q_dense = np.geomspace(0.002, 2.50, 200)
    I_net = m02.calc_primary_fractal_intensity(q_dense, p['Rmean'], p['stdDev'], p['phi'], p['D'], p['xi'], rho2=m02.RHO2_CALIBRATED, integration_steps=100)
    I_por = m02.calc_porod_intensity(q_dense, p['A_porod'], 4.0)
    I_pnc = m02.calc_pnc_intensity(q_dense, p['R_pnc'], p['std_pnc'], p['phi_pnc'], rho2=m02.RHO2_CALIBRATED, integration_steps=80)
    I_bkg = np.full_like(q_dense, p['bkg'])
    I_tot = I_net + I_por + I_pnc + I_bkg

    curves = [
        {'q': q_dense, 'y': I_por, 'color': SAXS_COLORS['droplets_porod'], 'width': 1.5, 'dash': '3,2'},
        {'q': q_dense, 'y': I_net, 'color': SAXS_COLORS['fractal_network'], 'width': 1.5, 'dash': '4,2'},
        {'q': q_dense, 'y': I_pnc, 'color': SAXS_COLORS['pnc_fluctuations'], 'width': 1.5, 'dash': '2,2'},
        {'q': q_dense, 'y': I_bkg, 'color': SAXS_COLORS['background'], 'width': 1.0, 'dash': '1,2'},
        {'q': q_dense, 'y': I_tot, 'color': SAXS_COLORS['total_fit'], 'width': 2.0},
    ]
    return q_data, y_data, curves


def build_curves_nd_i_009():
    sdir = os.path.join(RESULTS_DIR, 'ND-i-009_40mM_pH10.0')
    m09 = load_model_module(sdir, 'svg_m09')
    with open(os.path.join(sdir, 'output', 'composite', 'baseline_parameters.json')) as f:
        p = json.load(f)
    q_data, y_data = load_series_data(sdir, 'dataNDi009_fitting_results_merged.csv')

    from lmfit import Parameters
    params = Parameters()
    for k, v in p.items():
        if not k.startswith('std_err_'):
            params.add(k, value=v)

    q_dense = np.geomspace(0.002, 2.50, 200)
    I_tot, I_prim, I_large, I_oz, I_bkg = m09.decompose_model_components(params, q_dense)

    curves = [
        {'q': q_dense, 'y': I_large, 'color': SAXS_COLORS['droplets_porod'], 'width': 1.5, 'dash': '3,2'},
        {'q': q_dense, 'y': I_prim, 'color': SAXS_COLORS['fractal_network'], 'width': 1.5, 'dash': '4,2'},
        {'q': q_dense, 'y': I_oz, 'color': SAXS_COLORS['pnc_fluctuations'], 'width': 1.5, 'dash': '2,2'},
        {'q': q_dense, 'y': I_bkg, 'color': SAXS_COLORS['background'], 'width': 1.0, 'dash': '1,2'},
        {'q': q_dense, 'y': I_tot, 'color': SAXS_COLORS['total_fit'], 'width': 2.0},
    ]
    return q_data, y_data, curves


def generate_diagram_svg():
    print("Generating authentic composite curve panels for vector diagram...")
    q03, y03, c03 = build_curves_nd_i_003()
    q04, y04, c04 = build_curves_nd_i_004()
    q02, y02, c02 = build_curves_nd_i_002()
    q09, y09, c09 = build_curves_nd_i_009()

    panel_10mM = compute_svg_panel(q03, y03, c03, width=216, height=92, I_range=(1e-4, 5e4))
    panel_15mM = compute_svg_panel(q04, y04, c04, width=216, height=92, I_range=(1e-4, 5e4))
    panel_20mM = compute_svg_panel(q02, y02, c02, width=216, height=92, I_range=(1e-4, 5e4))
    panel_40mM = compute_svg_panel(q09, y09, c09, width=386, height=92, I_range=(1e-4, 1e6))

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1340 870" width="100%" height="100%" style="background-color: #0b1120; border-radius: 12px; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;">
  <defs>
    <!-- Background Gradients -->
    <linearGradient id="bgGrad" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#0b1120"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    <linearGradient id="cardGrad1" x1="0%" y1="0%" x2="0%" y2="100%">
      <stop offset="0%" stop-color="#16223b"/>
      <stop offset="100%" stop-color="#0f172a"/>
    </linearGradient>
    
    <!-- Arrow Markers -->
    <marker id="arrowhead" markerWidth="9" markerHeight="6" refX="8" refY="3" orient="auto">
      <polygon points="0 0, 9 3, 0 6" fill="#38bdf8"/>
    </marker>
    <marker id="arrowGreen" markerWidth="9" markerHeight="6" refX="8" refY="3" orient="auto">
      <polygon points="0 0, 9 3, 0 6" fill="#4ade80"/>
    </marker>
    <marker id="arrowRed" markerWidth="9" markerHeight="6" refX="8" refY="3" orient="auto">
      <polygon points="0 0, 9 3, 0 6" fill="#f87171"/>
    </marker>
    <marker id="arrowPurple" markerWidth="9" markerHeight="6" refX="8" refY="3" orient="auto">
      <polygon points="0 0, 9 3, 0 6" fill="#a855f7"/>
    </marker>
  </defs>

  <!-- Canvas Background -->
  <rect width="100%" height="100%" fill="url(#bgGrad)" stroke="#1e293b" stroke-width="1.5" rx="12"/>

  <!-- ==================== HEADER ==================== -->
  <g transform="translate(45, 30)">
    <text x="0" y="0" fill="#ffffff" font-size="19" font-weight="700" letter-spacing="-0.01em">Model transformation across supersaturation: viscoelastic phase separation (VPS) continuum</text>
    <text x="0" y="20" fill="#94a3b8" font-size="12">Structural and mathematical evolution grounded in multi-scale composite SAXS datasets (q in [0.002, 2.5] nm⁻¹)</text>
  </g>

  <!-- ==================== STANDARDISED COMPONENT PALETTE LEGEND ==================== -->
  <g transform="translate(45, 62)">
    <rect x="0" y="0" width="1250" height="25" rx="6" fill="#070d19" stroke="#1e293b" stroke-width="1"/>
    
    <circle cx="16" cy="12.5" r="3.5" fill="#64748b" opacity="0.8"/>
    <text x="26" y="16" fill="#94a3b8" font-size="9" font-weight="600">Experimental scattering data</text>
    
    <line x1="185" y1="12.5" x2="210" y2="12.5" stroke="{SAXS_COLORS['total_fit']}" stroke-width="2.5"/>
    <text x="218" y="16" fill="#f87171" font-size="9" font-weight="600">Total unified fit</text>

    <line x1="330" y1="12.5" x2="355" y2="12.5" stroke="{SAXS_COLORS['fractal_network']}" stroke-width="2.0" stroke-dasharray="4,2"/>
    <text x="363" y="16" fill="#34d399" font-size="9" font-weight="600">Spinodal fractal network</text>

    <line x1="520" y1="12.5" x2="545" y2="12.5" stroke="{SAXS_COLORS['acc_spheres']}" stroke-width="2.0" stroke-dasharray="5,2,1,2"/>
    <text x="553" y="16" fill="#38bdf8" font-size="9" font-weight="600">Discrete ACC spheres</text>

    <line x1="700" y1="12.5" x2="725" y2="12.5" stroke="{SAXS_COLORS['pnc_fluctuations']}" stroke-width="2.0" stroke-dasharray="2,2"/>
    <text x="733" y="16" fill="#c084fc" font-size="9" font-weight="600">Tanaka PNC fluctuation mesh</text>

    <line x1="915" y1="12.5" x2="940" y2="12.5" stroke="{SAXS_COLORS['droplets_porod']}" stroke-width="2.0" stroke-dasharray="3,2"/>
    <text x="948" y="16" fill="#fbbf24" font-size="9" font-weight="600">Macroscopic droplets / Porod tail</text>

    <line x1="1135" y1="12.5" x2="1160" y2="12.5" stroke="{SAXS_COLORS['background']}" stroke-width="1.2" stroke-dasharray="1,2"/>
    <text x="1168" y="16" fill="#94a3b8" font-size="9" font-weight="600">Solvent background</text>
  </g>

  <!-- ==================== THERMODYNAMIC CONCENTRATION AXIS (BOTTOM) ==================== -->
  <g transform="translate(0, 735)">
    <line x1="50" y1="0" x2="1285" y2="0" stroke="#334155" stroke-width="3"/>
    <line x1="50" y1="0" x2="1280" y2="0" stroke="#38bdf8" stroke-width="2.5" marker-end="url(#arrowhead)"/>

    <text x="660" y="42" text-anchor="middle" fill="#38bdf8" font-size="12" font-weight="600" letter-spacing="0.02em">
      Initial reactant concentration [Ca²⁺]₀ (mM) → increasing thermodynamic supersaturation ln(S_ACC)
    </text>

    <!-- 10 mM -->
    <circle cx="165" cy="0" r="5" fill="#38bdf8"/>
    <text x="165" y="20" text-anchor="middle" fill="#f1f5f9" font-size="11.5" font-weight="600">10 mM (pure)</text>
    <text x="165" y="34" text-anchor="middle" fill="#94a3b8" font-size="10">ln(S) ≈ 4.1</text>

    <!-- 15 mM -->
    <circle cx="445" cy="0" r="5" fill="#2dd4bf"/>
    <text x="445" y="20" text-anchor="middle" fill="#f1f5f9" font-size="11.5" font-weight="600">15 mM (pure)</text>
    <text x="445" y="34" text-anchor="middle" fill="#94a3b8" font-size="10">ln(S) ≈ 4.5</text>

    <!-- 20 mM -->
    <circle cx="725" cy="0" r="5" fill="#4ade80"/>
    <text x="725" y="20" text-anchor="middle" fill="#f1f5f9" font-size="11.5" font-weight="600">20 mM (pure)</text>
    <text x="725" y="34" text-anchor="middle" fill="#94a3b8" font-size="10">ln(S) ≈ 4.8</text>

    <!-- 40 mM -->
    <circle cx="1060" cy="0" r="5" fill="#f87171"/>
    <text x="1060" y="20" text-anchor="middle" fill="#f1f5f9" font-size="11.5" font-weight="600">40 mM (primary pure CaCO₃ benchmark)</text>
    <text x="1060" y="34" text-anchor="middle" fill="#94a3b8" font-size="10">ln(S) ≈ 5.0 to 5.5 (ND-i-009 / ND-i-001)</text>
  </g>

  <!-- ==================== SPINODAL BOUNDARY DIVIDER ==================== -->
  <g transform="translate(305, 0)">
    <line x1="0" y1="96" x2="0" y2="735" stroke="#f59e0b" stroke-width="1.8" stroke-dasharray="6,4" opacity="0.85"/>
    <rect x="-115" y="96" width="230" height="24" rx="12" fill="#1e1b18" stroke="#f59e0b" stroke-width="1.2"/>
    <text x="0" y="112" text-anchor="middle" fill="#fbbf24" font-size="10.5" font-weight="600">Spinodal boundary C* ≈ 12 to 15 mM</text>
    <text x="-130" y="114" text-anchor="end" fill="#38bdf8" font-size="11" font-weight="600">Sub-spinodal nucleation</text>
    <text x="130" y="114" text-anchor="start" fill="#f87171" font-size="11" font-weight="600">Spinodal demixing &amp; VPS</text>
  </g>

  <!-- ==================== REGIME 1: SUB-SPINODAL (10 mM, ND-i-003) ==================== -->
  <g transform="translate(45, 130)">
    <rect width="240" height="525" rx="10" fill="url(#cardGrad1)" stroke="#38bdf8" stroke-width="1.6"/>
    
    <rect x="0" y="0" width="240" height="40" rx="10" fill="#0c4a6e" fill-opacity="0.4"/>
    <text x="14" y="24" fill="#ffffff" font-size="13" font-weight="600">10 mM (ND-i-003)</text>
    <rect x="135" y="10" width="95" height="18" rx="9" fill="rgba(56, 189, 248, 0.2)" stroke="#38bdf8" stroke-width="1"/>
    <text x="182" y="22" text-anchor="middle" fill="#38bdf8" font-size="8.5" font-weight="600">Discrete spheres</text>

    <text x="14" y="58" fill="#94a3b8" font-size="10" font-weight="600">Morphology and mechanisms</text>
    <text x="14" y="74" fill="#e2e8f0" font-size="11" font-weight="500">Discrete amorphous nanospheres</text>
    <text x="14" y="89" fill="#94a3b8" font-size="9.5">Isolated ACC spheres and sheared droplets</text>

    <!-- Mini SAXS Curve Frame (Authentic Composite Data) -->
    <g transform="translate(12, 102)">
      {panel_10mM}
      <text x="14" y="14" fill="{SAXS_COLORS['acc_spheres']}" font-size="8" font-weight="600">Bessel oscillations</text>
      <text x="130" y="84" fill="{SAXS_COLORS['droplets_porod']}" font-size="8">Droplet tail</text>
    </g>

    <!-- Mathematical Formula Box -->
    <g transform="translate(12, 206)">
      <rect width="216" height="84" rx="6" fill="#080e1e" stroke="#233554" stroke-width="1"/>
      <text x="10" y="18" fill="#38bdf8" font-size="10" font-weight="600" font-family="monospace">I_tot(q) = I_sph + I_drop + bkg</text>
      <text x="10" y="34" fill="#e2e8f0" font-size="9" font-family="monospace">I_sph = (Δρ)² ∫ n(R) V² P(q,R) dR</text>
      <text x="10" y="50" fill="#fbbf24" font-size="9" font-family="monospace">S(q) ≡ 1.0  [No fractal network]</text>
      <text x="10" y="66" fill="#94a3b8" font-size="8.5" font-family="monospace">I_drop ≈ 27 φ_lg / (8π R_lg q⁴)</text>
    </g>

    <!-- Quantitative Metrics & Physical Interpretation -->
    <g transform="translate(12, 302)">
      <text x="0" y="14" fill="#38bdf8" font-size="10" font-weight="600">Parameters and observations:</text>
      <text x="0" y="31" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">ACC radius:</tspan> ⟨R⟩: 60 to 107 nm (growth)</text>
      <text x="0" y="47" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Polydispersity:</tspan> σ ≈ 10 to 22 nm (20%)</text>
      <text x="0" y="63" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Precipitate φ:</tspan> 2.3e-5 to 4.3e-4 (86%)</text>
      <text x="0" y="79" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Induction period:</tspan> t_ind ≈ 1.41 s (slow)</text>
      <text x="0" y="95" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Droplets:</tspan> R_lg ≈ 2.17 μm, φ ≈ 5.1e-5</text>
      <text x="0" y="113" fill="#94a3b8" font-size="9" font-style="italic">Classical nucleation; no fractal network</text>
      <text x="0" y="130" fill="#38bdf8" font-size="8.5" font-weight="600">VPS: Sub-critical concentration</text>
    </g>

    <!-- Phase Proportions & Dominance Evolution -->
    <g transform="translate(12, 444)">
      <rect width="216" height="70" rx="6" fill="#080e1e" stroke="#1e293b" stroke-width="1"/>
      <text x="8" y="13" fill="#38bdf8" font-size="8.5" font-weight="600">Relative phase proportions &amp; dominance:</text>
      
      <!-- Segmented Bar (width 200, height 10, rx 2) -->
      <g transform="translate(8, 20)">
        <rect x="0" y="0" width="57" height="10" rx="2" fill="{SAXS_COLORS['acc_spheres']}"/>
        <rect x="58" y="0" width="142" height="10" rx="2" fill="{SAXS_COLORS['acc_spheres_mode2']}"/>
      </g>
      <text x="8" y="42" fill="#e2e8f0" font-size="8.5" font-weight="500">ACC spheres: <tspan fill="#38bdf8" font-weight="600">100%</tspan> | Fractal gel: <tspan fill="#64748b" font-weight="600">0%</tspan> | Tanaka PNC: <tspan fill="#64748b" font-weight="600">0%</tspan></text>
      <text x="8" y="56" fill="#94a3b8" font-size="8">Mode 1 (28% @ 21 nm) + Mode 2 (72% @ 60 nm); 0% gel network &amp; PNC mesh</text>
    </g>
  </g>

  <!-- Transition Arrow 1 -> 2 -->
  <g transform="translate(285, 375)">
    <line x1="5" y1="0" x2="35" y2="0" stroke="#38bdf8" stroke-width="2" marker-end="url(#arrowhead)"/>
    <text x="20" y="-8" text-anchor="middle" fill="#38bdf8" font-size="8" font-weight="600">Spinodal</text>
    <text x="20" y="14" text-anchor="middle" fill="#94a3b8" font-size="7.5">crossover</text>
  </g>

  <!-- ==================== REGIME 2: CROSSOVER BOUNDARY (15 mM, ND-i-004) ==================== -->
  <g transform="translate(325, 130)">
    <rect width="240" height="525" rx="10" fill="url(#cardGrad1)" stroke="#2dd4bf" stroke-width="1.6"/>
    
    <rect x="0" y="0" width="240" height="40" rx="10" fill="#134e4a" fill-opacity="0.4"/>
    <text x="14" y="24" fill="#ffffff" font-size="13" font-weight="600">15 mM (ND-i-004)</text>
    <rect x="135" y="10" width="95" height="18" rx="9" fill="rgba(45, 212, 191, 0.2)" stroke="#2dd4bf" stroke-width="1"/>
    <text x="182" y="22" text-anchor="middle" fill="#2dd4bf" font-size="8.5" font-weight="600">Crossover</text>

    <text x="14" y="58" fill="#94a3b8" font-size="10" font-weight="600">Morphology and mechanisms</text>
    <text x="14" y="74" fill="#e2e8f0" font-size="11" font-weight="500">Hybrid boundary state</text>
    <text x="14" y="89" fill="#94a3b8" font-size="9.5">Primary network nodes and ACC spheres</text>

    <!-- Mini SAXS Curve Frame (Authentic Composite Data) -->
    <g transform="translate(12, 102)">
      {panel_15mM}
      <text x="14" y="14" fill="{SAXS_COLORS['fractal_network']}" font-size="8" font-weight="600">Dual knee</text>
      <text x="110" y="84" fill="{SAXS_COLORS['acc_spheres']}" font-size="8">Nodes &amp; Spheres</text>
    </g>

    <!-- Mathematical Formula Box -->
    <g transform="translate(12, 206)">
      <rect width="216" height="84" rx="6" fill="#080e1e" stroke="#233554" stroke-width="1"/>
      <text x="10" y="18" fill="#2dd4bf" font-size="10" font-weight="600" font-family="monospace">I_tot = I_sph + I_net + I_drop</text>
      <text x="10" y="34" fill="#e2e8f0" font-size="9" font-family="monospace">I_net = P_prim(q) · S_net(q)</text>
      <text x="10" y="50" fill="#38bdf8" font-size="9" font-family="monospace">I_sph: ACC spheres (R ≈ 48 nm)</text>
      <text x="10" y="66" fill="#94a3b8" font-size="8.5" font-family="monospace">S_net: D_m = 2.85, ξ = 108.7 nm</text>
    </g>

    <!-- Quantitative Metrics -->
    <g transform="translate(12, 302)">
      <text x="0" y="14" fill="#2dd4bf" font-size="10" font-weight="600">Parameters and observations:</text>
      <text x="0" y="31" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">ACC spheres:</tspan> ⟨R_acc⟩: 15.9 to 43.8 nm</text>
      <text x="0" y="47" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Primary nodes:</tspan> R_prim ≈ 4.94 nm</text>
      <text x="0" y="63" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Correlation ξ:</tspan> 20.2 to 196.9 nm</text>
      <text x="0" y="79" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Fractal dimension:</tspan> D_m ≈ 2.85 to 2.33</text>
      <text x="0" y="95" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Induction period:</tspan> t_ind ≈ 0.35 s</text>
      <text x="0" y="113" fill="#94a3b8" font-size="9" font-style="italic">Coexistence of network and discrete spheres</text>
      <text x="0" y="130" fill="#2dd4bf" font-size="8.5" font-weight="600">VPS: Onset of dynamic asymmetry</text>
    </g>

    <!-- Phase Proportions & Dominance Evolution -->
    <g transform="translate(12, 444)">
      <rect width="216" height="70" rx="6" fill="#080e1e" stroke="#1e293b" stroke-width="1"/>
      <text x="8" y="13" fill="#2dd4bf" font-size="8.5" font-weight="600">Relative phase proportions &amp; dominance:</text>
      
      <!-- Segmented Bar (width 200, height 10, rx 2) -->
      <g transform="translate(8, 20)">
        <rect x="0" y="0" width="126" height="10" rx="2" fill="{SAXS_COLORS['acc_spheres']}"/>
        <rect x="127" y="0" width="73" height="10" rx="2" fill="{SAXS_COLORS['fractal_network']}"/>
      </g>
      <text x="8" y="42" fill="#e2e8f0" font-size="8.5" font-weight="500">Baseline: <tspan fill="#38bdf8" font-weight="600">63% ACC spheres</tspan> | <tspan fill="#34d399" font-weight="600">37% fractal gel</tspan></text>
      <text x="8" y="56" fill="#fbbf24" font-size="8">Dynamic evolution: 96% ACC (0.06 s) &rarr; 76% fractal gel (10.9 s)</text>
    </g>
  </g>

  <!-- Transition Arrow 2 -> 3 -->
  <g transform="translate(565, 375)">
    <line x1="5" y1="0" x2="35" y2="0" stroke="#4ade80" stroke-width="2" marker-end="url(#arrowGreen)"/>
    <text x="20" y="-8" text-anchor="middle" fill="#4ade80" font-size="8" font-weight="600">Percolation</text>
    <text x="20" y="14" text-anchor="middle" fill="#94a3b8" font-size="7.5">gelation</text>
  </g>

  <!-- ==================== REGIME 3: MASS FRACTAL GEL (20 mM, ND-i-002) ==================== -->
  <g transform="translate(605, 130)">
    <rect width="240" height="525" rx="10" fill="url(#cardGrad1)" stroke="#4ade80" stroke-width="1.6"/>
    
    <rect x="0" y="0" width="240" height="40" rx="10" fill="#14532d" fill-opacity="0.4"/>
    <text x="14" y="24" fill="#ffffff" font-size="13" font-weight="600">20 mM (ND-i-002)</text>
    <rect x="135" y="10" width="95" height="18" rx="9" fill="rgba(74, 222, 128, 0.2)" stroke="#4ade80" stroke-width="1"/>
    <text x="182" y="22" text-anchor="middle" fill="#4ade80" font-size="8.5" font-weight="600">Fractal gel</text>

    <text x="14" y="58" fill="#94a3b8" font-size="10" font-weight="600">Morphology and mechanisms</text>
    <text x="14" y="74" fill="#e2e8f0" font-size="11" font-weight="500">Percolated mass fractal network</text>
    <text x="14" y="89" fill="#94a3b8" font-size="9.5">Continuous gel mesh; spheres suppressed</text>

    <!-- Mini SAXS Curve Frame (Authentic Composite Data) -->
    <g transform="translate(12, 102)">
      {panel_20mM}
      <text x="14" y="14" fill="{SAXS_COLORS['fractal_network']}" font-size="8" font-weight="600">q⁻¹·⁸⁶ fractal</text>
      <text x="120" y="84" fill="{SAXS_COLORS['pnc_fluctuations']}" font-size="8">PNC shoulder</text>
    </g>

    <!-- Mathematical Formula Box -->
    <g transform="translate(12, 206)">
      <rect width="216" height="84" rx="6" fill="#080e1e" stroke="#233554" stroke-width="1"/>
      <text x="10" y="18" fill="#4ade80" font-size="10" font-weight="600" font-family="monospace">I_tot = P_prim · S_net + A_p/q⁴</text>
      <text x="10" y="34" fill="#e2e8f0" font-size="9" font-family="monospace">S(q, D_m, ξ) Teixeira model</text>
      <text x="10" y="50" fill="#f87171" font-size="9" font-family="monospace">I_sphere → 0  [Spheres vanish]</text>
      <text x="10" y="66" fill="#94a3b8" font-size="8.5" font-family="monospace">D_m = 1.86, ξ = 356 nm</text>
    </g>

    <!-- Quantitative Metrics -->
    <g transform="translate(12, 302)">
      <text x="0" y="14" fill="#4ade80" font-size="10" font-weight="600">Parameters and observations:</text>
      <text x="0" y="31" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Primary nodes:</tspan> R_prim ≈ 11.2 nm</text>
      <text x="0" y="47" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Fractal dimension:</tspan> D_m ≈ 1.86</text>
      <text x="0" y="63" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Correlation ξ:</tspan> 22.7 to 356 nm (spanning)</text>
      <text x="0" y="79" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Induction period:</tspan> t_ind &lt; 0.08 s (fast)</text>
      <text x="0" y="95" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Droplets:</tspan> R_lg ≈ 7.16 μm, φ ≈ 8.4e-5</text>
      <text x="0" y="113" fill="#94a3b8" font-size="9" font-style="italic">Spinodal gelation; discrete spheres absent</text>
      <text x="0" y="130" fill="#4ade80" font-size="8.5" font-weight="600">VPS: Transient elastic network</text>
    </g>

    <!-- Phase Proportions & Dominance Evolution -->
    <g transform="translate(12, 444)">
      <rect width="216" height="70" rx="6" fill="#080e1e" stroke="#1e293b" stroke-width="1"/>
      <text x="8" y="13" fill="#4ade80" font-size="8.5" font-weight="600">Relative phase proportions &amp; dominance:</text>
      
      <!-- Segmented Bar (width 200, height 10, rx 2) -->
      <g transform="translate(8, 20)">
        <rect x="0" y="0" width="152" height="10" rx="2" fill="{SAXS_COLORS['fractal_network']}"/>
        <rect x="153" y="0" width="47" height="10" rx="2" fill="{SAXS_COLORS['pnc_fluctuations']}"/>
      </g>
      <text x="8" y="42" fill="#e2e8f0" font-size="8.5" font-weight="500">Fractal gel: <tspan fill="#34d399" font-weight="600">76%</tspan> | Tanaka PNC mesh: <tspan fill="#c084fc" font-weight="600">24%</tspan> (0% spheres)</text>
      <text x="8" y="56" fill="#94a3b8" font-size="8">Transient viscoelastic PNC mesh (R_pnc &approx; 1.46 nm); discrete spheres 0%</text>
    </g>
  </g>

  <!-- Transition Arrow 3 -> 4 -->
  <g transform="translate(845, 375)">
    <line x1="5" y1="0" x2="35" y2="0" stroke="#f87171" stroke-width="2" marker-end="url(#arrowRed)"/>
    <text x="20" y="-8" text-anchor="middle" fill="#f87171" font-size="8" font-weight="600">Deep LLPS</text>
    <text x="20" y="14" text-anchor="middle" fill="#94a3b8" font-size="7.5">demixing</text>
  </g>

  <!-- ==================== REGIME 4: 40 mM PURE CaCO3 BENCHMARK (ND-i-009) ==================== -->
  <g transform="translate(885, 130)">
    <rect width="410" height="525" rx="10" fill="url(#cardGrad1)" stroke="#f87171" stroke-width="1.8"/>
    
    <rect x="0" y="0" width="410" height="40" rx="10" fill="#7f1d1d" fill-opacity="0.4"/>
    <text x="14" y="24" fill="#ffffff" font-size="13" font-weight="600">40 mM (ND-i-009, primary pure CaCO₃ benchmark)</text>
    <rect x="300" y="10" width="98" height="18" rx="9" fill="rgba(248, 113, 113, 0.2)" stroke="#f87171" stroke-width="1"/>
    <text x="349" y="22" text-anchor="middle" fill="#f87171" font-size="8.5" font-weight="600">VPS benchmark</text>

    <text x="14" y="58" fill="#94a3b8" font-size="10" font-weight="600">Morphology and mechanisms</text>
    <text x="14" y="74" fill="#e2e8f0" font-size="11" font-weight="500">Percolating spinodal network, coarsening droplets, and PNC mesh</text>
    <text x="14" y="89" fill="#94a3b8" font-size="9.5">Time-resolved percolation (351 frames, 15 ms/frame; compares to ND-i-001 at pH 11.3)</text>

    <!-- Mini SAXS Curve Frame (Authentic Composite Data) -->
    <g transform="translate(12, 102)">
      {panel_40mM}
      <text x="14" y="14" fill="{SAXS_COLORS['droplets_porod']}" font-size="8" font-weight="600">Droplet Porod tail</text>
      <text x="150" y="44" fill="{SAXS_COLORS['fractal_network']}" font-size="8" font-weight="600">q⁻¹·⁹² continuous gel</text>
      <text x="260" y="84" fill="{SAXS_COLORS['pnc_fluctuations']}" font-size="8" font-weight="600">PNC viscoelastic mesh</text>
    </g>

    <!-- Mathematical Formula Box -->
    <g transform="translate(12, 206)">
      <rect width="386" height="84" rx="6" fill="#080e1e" stroke="#233554" stroke-width="1"/>
      <text x="10" y="18" fill="#f87171" font-size="10" font-weight="600" font-family="monospace">I_tot(q) = P_prim · S_net + I_drop + I_fluc(q) + bkg</text>
      <text x="10" y="34" fill="#e2e8f0" font-size="9" font-family="monospace">S_net: continuous D_m(t) ∈ [1.881, 1.925], ξ(t): 46.7 nm → ∞</text>
      <text x="10" y="50" fill="#fbbf24" font-size="9" font-family="monospace">I_drop: phase inversion droplets ⟨R_lg⟩: 2.27 → 3.20 μm</text>
      <text x="10" y="66" fill="#a855f7" font-size="9" font-family="monospace">I_fluc = A_fluc / [1 + (q ξ_fluc)²/6]²  [PNC viscoelastic mesh]</text>
    </g>

    <!-- Quantitative Metrics & Physical Interpretation -->
    <g transform="translate(12, 302)">
      <text x="0" y="14" fill="#f87171" font-size="10" font-weight="600">Parameters and observations (Tanaka VPS verification):</text>
      <text x="0" y="31" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Time-resolved percolation:</tspan> finite ξ = 46.7 nm at frame 0 → continuous gel (ξ → ∞)</text>
      <text x="0" y="47" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Tanaka volume shrinking:</tspan> D_m(t) compacts from 1.881 to 1.925 (|dD_m| ≤ 0.004)</text>
      <text x="0" y="63" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">PNC viscoelastic mesh:</tspan> ξ_fluc: 2.72 to 2.84 nm; A_fluc: 0.142 to 0.197 cm⁻¹</text>
      <text x="0" y="79" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">Secondary droplets:</tspan> coarsening ⟨R_lg⟩: 2,274 to 3,196 nm (polydispersity 0.065)</text>
      <text x="0" y="95" fill="#cbd5e1" font-size="9.5">• <tspan font-weight="600">pH 11.3 pure counterpart (ND-i-001):</tspan> instant quench (D_m = 2.33, ξ = 735 nm, t_ind &lt; 0.06 s)</text>
      <text x="0" y="113" fill="#94a3b8" font-size="9" font-style="italic">Pure CaCO₃ benchmark resolving continuous pre-percolation and PNC mesh dynamics</text>
      <text x="0" y="130" fill="#f87171" font-size="8.5" font-weight="600">VPS: Stress relaxation and phase inversion</text>
    </g>

    <!-- Phase Proportions & Dominance Evolution -->
    <g transform="translate(12, 444)">
      <rect width="386" height="70" rx="6" fill="#080e1e" stroke="#1e293b" stroke-width="1"/>
      <text x="8" y="13" fill="#f87171" font-size="8.5" font-weight="600">Relative phase proportions &amp; dominance (pure benchmark):</text>
      
      <!-- Segmented Bar (width 370, height 10, rx 2) -->
      <g transform="translate(8, 20)">
        <rect x="0" y="0" width="270" height="10" rx="2" fill="{SAXS_COLORS['fractal_network']}"/>
        <rect x="271" y="0" width="95" height="10" rx="2" fill="{SAXS_COLORS['pnc_fluctuations']}"/>
        <rect x="367" y="0" width="3" height="10" rx="2" fill="{SAXS_COLORS['droplets_porod']}"/>
      </g>
      <text x="8" y="42" fill="#e2e8f0" font-size="8.5" font-weight="500">Fractal gel: <tspan fill="#34d399" font-weight="600">73%</tspan> | Tanaka PNC mesh: <tspan fill="#c084fc" font-weight="600">26%</tspan> | Droplets: <tspan fill="#fbbf24" font-weight="600">1%</tspan> (0% spheres)</text>
      <text x="8" y="56" fill="#94a3b8" font-size="8">Tanaka viscoelastic fluctuation mesh (ξ_fluc ≈ 2.7 nm) vitrifies into percolated gel nodes</text>
    </g>
  </g>

  <!-- ==================== ORTHOGONAL ADDITIVE BRANCH (Co2+, ND-i-005) ==================== -->
  <g transform="translate(885, 12)">
    <rect width="410" height="78" rx="8" fill="#1e1035" stroke="#a855f7" stroke-width="1.3"/>
    <text x="12" y="18" fill="#c084fc" font-size="11" font-weight="600">Chemical perturbation branch: Co²⁺ dopant (ND-i-005)</text>
    <rect x="290" y="8" width="110" height="16" rx="8" fill="rgba(168, 85, 247, 0.2)" stroke="#a855f7" stroke-width="0.8"/>
    <text x="345" y="19" text-anchor="middle" fill="#c084fc" font-size="8" font-weight="600">Additive modification</text>
    
    <text x="12" y="36" fill="#e2e8f0" font-size="9.5">• <tspan font-weight="600">Adsorption lowers interfacial tension γ:</tspan> enlarges nodes R_prim: 3.5 to 6.59 nm</text>
    <text x="12" y="52" fill="#e2e8f0" font-size="9.5">• <tspan font-weight="600">Extended network cutoff:</tspan> expands ξ: 735 to 1,190 nm; open fractal D_m: 2.33 to 2.10</text>
    <text x="12" y="68" fill="#94a3b8" font-size="9">• Only sample containing transition metal dopant; all other series are pure CaCO₃</text>
  </g>

  <!-- ==================== CALLOUT 1: UNIVERSAL SHEAR DROPLETS ==================== -->
  <g transform="translate(45, 668)">
    <rect width="580" height="52" rx="6" fill="#111827" stroke="#374151" stroke-width="1"/>
    <circle cx="20" cy="26" r="5" fill="#d97706"/>
    <text x="34" y="21" fill="#f1f5f9" font-size="10.5" font-weight="600">Stationary phase inversion droplets (universal across series)</text>
    <text x="34" y="38" fill="#94a3b8" font-size="9.5">⟨R_large⟩ ≈ 2.2 to 3.5 μm, φ_large ≈ (3.0 to 6.0) × 10⁻⁵ formed by mechanical strut rupture and capillary break-up</text>
  </g>

  <!-- ==================== CALLOUT 2: VISCOELASTIC MESH & PNCs ==================== -->
  <g transform="translate(635, 668)">
    <rect width="660" height="52" rx="6" fill="#141c2e" stroke="#9333ea" stroke-width="1.2"/>
    <circle cx="20" cy="26" r="5" fill="#9333ea"/>
    <text x="34" y="21" fill="#c084fc" font-size="10.5" font-weight="600">Viscoelastic mesh in the context of prenucleation clusters (PNCs)</text>
    <text x="34" y="38" fill="#cbd5e1" font-size="9.5">High-q fluctuation length ξ_fluc ≈ 1.5 to 3.0 nm reflects the transient mesh of entangled, hydrated PNC oligomers prior to solid ACC condensation</text>
  </g>

  <!-- ==================== FOOTER: ATTRIBUTION ==================== -->
  <g transform="translate(45, 840)">
    <text x="0" y="0" fill="#64748b" font-size="10">Data processing responsible: Tomasz Stawski, tomasz.stawski@bam.de | Viscoelastic phase separation (VPS) framework | Data: ESRF ID02</text>
  </g>
</svg>"""
    return svg


def main():
    svg_content = generate_diagram_svg()
    out_path = os.path.join(OUT_DIR, 'model_transformation_diagram.svg')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    print(f"Saved authentic, data-driven vector SVG diagram to: {out_path}")


if __name__ == '__main__':
    main()
