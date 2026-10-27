#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-002.
-----------------------------------------------------------------------
Generates a standalone, responsive HTML dashboard for reviewing:
1. Overview & Composite Anchor Fit (q in [0.002, 2.5] nm^-1)
2. Unified Cross-Scale Kinetic Trajectories (0.060 s to 10.877 s)
3. Short-Distance Time-Resolved Inspector (3.0 m, Frames 0 to 49, with live slider & +/- 1 sigma telemetry)
4. Long-Distance Time-Resolved Inspector (31.0 m, Frames 0 to 49, with live slider & +/- 1 sigma telemetry)
5. Comprehensive Parameter Tables (searchable & scrollable)
6. Physical Model Formulations & Equations rendered via MathJax.

Output:
  output/dashboard/index.html
"""

import os
import sys
import json
import pandas as pd

_curr = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(_curr, '..', '..')) if os.path.basename(_curr) in ['analysis', 'plotting', 'dashboard', 'reporting'] else os.path.abspath(os.path.join(_curr, '..'))
sys.path.insert(0, PROJECT_DIR)


def build_dashboard():
    dashboard_dir = os.path.join(PROJECT_DIR, 'output', 'dashboard')
    os.makedirs(dashboard_dir, exist_ok=True)

    # 1. Load Short Distance Parameters
    sd_csv = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    df_sd = pd.read_csv(sd_csv)
    short_telemetry = {}
    for _, r in df_sd.iterrows():
        f = int(r['frame_index'])
        t_sec = float(r['time_seconds'])
        short_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': f"{float(r['chi2_red']):.3f}",
            'Rmean': f"{float(r['Rmean']):.2f} ± {float(r.get('std_err_Rmean', 0.0)):.2f}",
            'stdDev': f"{float(r['stdDev']):.2f} ± {float(r.get('std_err_stdDev', 0.0)):.2f}",
            'D': f"{float(r['D']):.3f} ± {float(r.get('std_err_D', 0.0)):.3f}",
            'phi': f"{float(r['phi']) * 1e4:.2f} ± {float(r.get('std_err_phi', 0.0)) * 1e4:.2f}",
            'A_fluc': f"{float(r['A_fluc']):.4f} ± {float(r.get('std_err_A_fluc', 0.0)):.4f}",
            'xi_fluc': f"{float(r['xi_fluc']):.2f} ± {float(r.get('std_err_xi_fluc', 0.0)):.2f}",
            'bkg': f"{float(r['bkg']):.4f} ± {float(r.get('std_err_bkg', 0.0)):.4f}"
        }

    # 2. Load Long Distance Parameters
    ld_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')
    df_ld = pd.read_csv(ld_csv)
    long_telemetry = {}
    for _, r in df_ld.iterrows():
        f = int(r['frame_index'])
        t_sec = float(r['time_seconds'])
        long_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': f"{float(r['chi2_red']):.3f}",
            'mean_radius_large': f"{float(r['mean_radius_large']):.0f} ± {float(r.get('std_err_mean_radius_large', 0.0)):.0f}",
            'std_radius_large': f"{float(r['std_radius_large']):.0f} ± {float(r.get('std_err_std_radius_large', 0.0)):.0f}",
            'total_volume_fraction_large': f"{float(r['total_volume_fraction_large']) * 1e5:.2f} ± {float(r.get('std_err_total_volume_fraction_large', 0.0)) * 1e5:.2f}",
            'mean_radius_large_2': f"{float(r.get('mean_radius_large_2', 0.0)):.0f} ± {float(r.get('std_err_mean_radius_large_2', 0.0)):.0f}",
            'std_radius_large_2': f"{float(r.get('std_radius_large_2', 0.0)):.0f} ± {float(r.get('std_err_std_radius_large_2', 0.0)):.0f}",
            'total_volume_fraction_large_2': f"{float(r.get('total_volume_fraction_large_2', 0.0)) * 1e5:.2f} ± {float(r.get('std_err_total_volume_fraction_large_2', 0.0)) * 1e5:.2f}",
            'D': f"{float(r['D']):.3f} ± {float(r.get('std_err_D', 0.0)):.3f}",
            'xi': f"{float(r['xi']):.0f} ± {float(r.get('std_err_xi', 0.0)):.0f}",
            'phi': f"{float(r['phi']) * 1e4:.2f} ± {float(r.get('std_err_phi', 0.0)) * 1e4:.2f}"
        }

    # 3. Load Baseline Parameters
    bp_json = os.path.join(PROJECT_DIR, 'output', 'composite', 'baseline_parameters.json')
    with open(bp_json, 'r') as f:
        bp = json.load(f)

    # 4. Generate HTML content
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ND-i-002 Multiscale Kinetic Analysis Dashboard</title>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  <style>
    :root {{
      --bg-primary: #0f172a;
      --bg-secondary: #1e293b;
      --bg-card: #1e293b;
      --text-primary: #f8fafc;
      --text-secondary: #94a3b8;
      --accent: #38bdf8;
      --accent-hover: #0ea5e9;
      --border: #334155;
      --card-border: rgba(255, 255, 255, 0.08);
      --success: #10b981;
      --warning: #f59e0b;
      --danger: #ef4444;
    }}
    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background-color: var(--bg-primary);
      color: var(--text-primary);
      line-height: 1.5;
    }}
    header {{
      background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
      border-bottom: 1px solid var(--border);
      padding: 1.5rem 2rem;
    }}
    .header-content {{
      max-width: 1400px;
      margin: 0 auto;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 1rem;
    }}
    .badge {{
      display: inline-block;
      padding: 0.25rem 0.6rem;
      border-radius: 9999px;
      font-size: 0.75rem;
      font-weight: 600;
      background: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      border: 1px solid rgba(56, 189, 248, 0.3);
    }}
    nav.tabs {{
      background: var(--bg-secondary);
      border-bottom: 1px solid var(--border);
      position: sticky;
      top: 0;
      z-index: 100;
    }}
    .tabs-container {{
      max-width: 1400px;
      margin: 0 auto;
      display: flex;
      overflow-x: auto;
    }}
    .tab-btn {{
      padding: 0.9rem 1.4rem;
      background: transparent;
      border: none;
      color: var(--text-secondary);
      font-size: 0.9rem;
      font-weight: 500;
      cursor: pointer;
      border-bottom: 2px solid transparent;
      white-space: nowrap;
      transition: all 0.2s ease;
    }}
    .tab-btn:hover {{
      color: var(--text-primary);
      background: rgba(255, 255, 255, 0.02);
    }}
    .tab-btn.active {{
      color: var(--accent);
      border-bottom-color: var(--accent);
      font-weight: 600;
    }}
    main {{
      max-width: 1400px;
      margin: 2rem auto;
      padding: 0 1.5rem;
    }}
    .tab-content {{
      display: none;
    }}
    .tab-content.active {{
      display: block;
    }}
    .card {{
      background: var(--bg-card);
      border: 1px solid var(--card-border);
      border-radius: 12px;
      padding: 1.5rem;
      margin-bottom: 1.5rem;
      box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
    }}
    .card-title {{
      font-size: 1.15rem;
      font-weight: 600;
      margin-bottom: 1rem;
      color: var(--text-primary);
    }}
    .inspector-grid {{
      display: grid;
      grid-template-columns: 1fr 340px;
      gap: 1.5rem;
    }}
    @media (max-width: 1024px) {{
      .inspector-grid {{
        grid-template-columns: 1fr;
      }}
    }}
    .viewer-panel img {{
      width: 100%;
      height: auto;
      border-radius: 8px;
      background: #000;
      display: block;
      border: 1px solid var(--border);
    }}
    .slider-box {{
      margin: 1rem 0;
      background: rgba(0, 0, 0, 0.2);
      padding: 1rem;
      border-radius: 8px;
    }}
    .slider-box input[type="range"] {{
      width: 100%;
      height: 6px;
      border-radius: 3px;
      background: #334155;
      outline: none;
      margin-top: 0.5rem;
    }}
    .btn-group {{
      display: flex;
      gap: 0.5rem;
      margin-bottom: 1rem;
    }}
    .btn {{
      flex: 1;
      padding: 0.5rem 0.8rem;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: rgba(255, 255, 255, 0.05);
      color: var(--text-primary);
      font-size: 0.85rem;
      cursor: pointer;
      font-weight: 500;
    }}
    .btn:hover {{
      background: rgba(255, 255, 255, 0.1);
    }}
    .telemetry-table {{
      width: 100%;
      border-collapse: collapse;
      margin-bottom: 1rem;
      font-size: 0.85rem;
    }}
    .telemetry-table th, .telemetry-table td {{
      padding: 0.45rem 0.6rem;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    .telemetry-table th {{
      text-align: left;
      color: var(--text-secondary);
      font-weight: 500;
    }}
    .telemetry-table td {{
      text-align: right;
      font-weight: 600;
      font-family: monospace;
    }}
    .data-table-container {{
      max-height: 550px;
      overflow-y: auto;
      border: 1px solid var(--border);
      border-radius: 8px;
    }}
    table.data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
      text-align: left;
    }}
    table.data-table th, table.data-table td {{
      padding: 0.6rem 0.8rem;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    table.data-table th {{
      background: #111827;
      position: sticky;
      top: 0;
      color: var(--accent);
      font-weight: 600;
    }}
    table.data-table tr:hover {{
      background: rgba(255, 255, 255, 0.03);
    }}
    .equation-box {{
      background: rgba(0, 0, 0, 0.3);
      padding: 1.2rem;
      border-radius: 8px;
      margin: 1rem 0;
      border-left: 3px solid var(--accent);
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-content">
      <div>
        <h1 style="font-size: 1.4rem; font-weight: 700;">ND-i-002: Unified Multiscale SAXS Kinetics Suite</h1>
        <p style="font-size: 0.85rem; color: var(--text-secondary); margin-top: 0.25rem;">
          $20.0\text{{ mM }}\text{{CaCO}}_3$, pH 11.3, Pure System without Additive &bull; High-Speed Synchrotron SAXS (3.0 m & 31.0 m)
        </p>
      </div>
      <div style="display: flex; gap: 0.5rem; align-items: center;">
        <a href="../../../index.html" class="badge" style="text-decoration: none; cursor: pointer;">&larr; Master Dashboard</a>
        <span class="badge">t = 0.060 s &rarr; 10.88 s</span>
        <span class="badge" style="background: rgba(16, 185, 129, 0.15); color: var(--success); border-color: rgba(16, 185, 129, 0.3);">
          &plusmn;1&sigma; Uncertainties
        </span>
      </div>
    </div>
  </header>

  <nav class="tabs">
    <div class="tabs-container">
      <button class="tab-btn active" onclick="switchTab('overview')">1. Overview & Composite Curve</button>
      <button class="tab-btn" onclick="switchTab('unified_trajectories')">2. Cross-Scale Kinetics</button>
      <button class="tab-btn" onclick="switchTab('short_inspector')">3. Short Distance (3.0 m)</button>
      <button class="tab-btn" onclick="switchTab('long_inspector')">4. Long Distance (31.0 m)</button>
      <button class="tab-btn" onclick="switchTab('tables')">5. Parameter Tables</button>
      <button class="tab-btn" onclick="switchTab('theory')">6. Physical Formalism</button>
    </div>
  </nav>

  <main>
    <!-- TAB 1: OVERVIEW & BASELINE -->
    <div id="overview" class="tab-content active">
      <div class="card">
        <div class="card-title">Full-range composite baseline anchor curve ($q \\in [0.002, 2.5]\\text{{ nm}}^{{-1}}$)</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem; font-size: 0.9rem;">
          Stitched synchrotron composite curve spanning over 3 decades of momentum transfer. Decomposed into primary polydisperse
          building blocks, inter-particle mass fractal network, bimodal macroscopic dense liquid droplets, and sub-nanometer PNC clusters.
        </p>
        <img src="../composite/composite_baseline_fit.png" alt="Composite Baseline Fit" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>

      <div class="card">
        <div class="card-title">Converged structural baseline parameters (&plusmn;1&sigma;)</div>
        <table class="telemetry-table" style="font-size: 0.9rem;">
          <thead>
            <tr style="border-bottom: 2px solid var(--border);">
              <th style="color: var(--accent);">Parameter</th>
              <th style="color: var(--accent);">Symbol</th>
              <th style="color: var(--accent); text-align: right;">Value</th>
              <th style="color: var(--accent); text-align: right;">&plusmn;1&sigma; Error</th>
              <th style="color: var(--accent); text-align: left; padding-left: 1.5rem;">Physical Interpretation</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>Primary Building Block Radius</td><td>&lang;R<sub>prim</sub>&rang;</td>
              <td>{bp['Rmean']:.2f} nm</td><td>&plusmn; {bp.get('std_err_Rmean', 0.0):.2f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of primary building blocks</td>
            </tr>
            <tr>
              <td>Primary Core Dispersion</td><td>&sigma;<sub>prim</sub></td>
              <td>{bp['stdDev']:.2f} nm</td><td>&plusmn; {bp.get('std_err_stdDev', 0.0):.2f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Log-normal dispersion width</td>
            </tr>
            <tr>
              <td>Mass Fractal Dimension</td><td>D<sub>m</sub></td>
              <td>{bp['D']:.3f}</td><td>&plusmn; {bp.get('std_err_D', 0.0):.3f}</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Network space-filling scaling exponent</td>
            </tr>
            <tr>
              <td>Network Cutoff Length</td><td>&xi;</td>
              <td>{bp['xi']:.0f} nm</td><td>&plusmn; {bp.get('std_err_xi', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Cluster correlation cutoff length</td>
            </tr>
            <tr>
              <td>Network Volume Fraction</td><td>&phi;<sub>net</sub></td>
              <td>{bp['phi'] * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_phi', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of primary fractal network</td>
            </tr>
            <tr>
              <td>Primary Droplet Radius (Mode 1)</td><td>&lang;R<sub>1</sub>&rang;</td>
              <td>{bp.get('mean_radius_large', 0.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_mean_radius_large', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of primary dense liquid droplets</td>
            </tr>
            <tr>
              <td>Primary Droplet Dispersion (Mode 1)</td><td>&sigma;<sub>1</sub></td>
              <td>{bp.get('std_radius_large', 0.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_std_radius_large', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Log-normal droplet dispersion width (Mode 1)</td>
            </tr>
            <tr>
              <td>Primary Droplet Fraction (Mode 1)</td><td>&phi;<sub>1</sub></td>
              <td>{bp.get('total_volume_fraction_large', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_total_volume_fraction_large', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of primary macroscopic droplets</td>
            </tr>
            <tr>
              <td>Coalesced Giant Droplet Radius (Mode 2)</td><td>&lang;R<sub>2</sub>&rang;</td>
              <td>{bp.get('mean_radius_large_2', 0.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_mean_radius_large_2', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of shear-coalesced giant droplets</td>
            </tr>
            <tr>
              <td>Coalesced Droplet Dispersion (Mode 2)</td><td>&sigma;<sub>2</sub></td>
              <td>{bp.get('std_radius_large_2', 0.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_std_radius_large_2', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Log-normal droplet dispersion width (Mode 2)</td>
            </tr>
            <tr>
              <td>Coalesced Droplet Fraction (Mode 2)</td><td>&phi;<sub>2</sub></td>
              <td>{bp.get('total_volume_fraction_large_2', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_total_volume_fraction_large_2', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of shear-coalesced giant droplets</td>
            </tr>
            <tr>
              <td>Porod Interface Scale</td><td>A<sub>porod</sub></td>
              <td>{bp.get('A_porod', 0.0) * 1e8:.2f} &times; 10<sup>-8</sup></td><td>&plusmn; {bp.get('std_err_A_porod', 0.0) * 1e8:.2f} &times; 10<sup>-8</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Residual Porod surface interface scattering</td>
            </tr>
            <tr>
              <td>Fluctuation Correlation Length</td><td>&xi;<sub>fluc</sub></td>
              <td>{bp.get('xi_fluc', 7.84):.2f} nm</td><td>&plusmn; {bp.get('std_err_xi_fluc', 0.0):.2f} nm</td>
            </tr>
            <tr>
              <td>Fluctuation Amplitude</td><td>A<sub>fluc</sub></td>
              <td>{bp.get('A_fluc', 0.02):.4f} cm<sup>-1</sup></td><td>&plusmn; {bp.get('std_err_A_fluc', 0.0):.4f} cm<sup>-1</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Tanaka fluctuation mesh forward intensity</td>
            </tr>
            <tr>
              <td>Incoherent Background</td><td>bkg</td>
              <td>{bp['bkg']:.4e} cm<sup>-1</sup></td><td>&plusmn; {bp.get('std_err_bkg', 0.0):.4e}</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Flat background from solvent/incoherent scatter</td>
            </tr>
            <tr>
              <td>Goodness of Fit</td><td>&chi;<sup>2</sup><sub>red</sub></td>
              <td style="color: var(--success); font-weight: 700;">{bp['redchi']:.3f}</td><td>&mdash;</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Full composite curve reduced chi-squared</td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>

    <!-- TAB 2: UNIFIED TRAJECTORIES -->
    <div id="unified_trajectories" class="tab-content">
      <div class="card">
        <div class="card-title">Cross-scale kinetic trajectories synthesis ($t = 0.060\text{{ s}} \to 10.88\text{{ s}}$)</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem; font-size: 0.9rem;">
          Cross-scale synthesis combining 3.0 m (WAXS/MAXS) and 31.0 m (SAXS) measurements with shaded &plusmn;1&sigma; uncertainty envelopes.
        </p>
        <img src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Trajectories" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>
    </div>

    <!-- TAB 3: SHORT DISTANCE INSPECTOR -->
    <div id="short_inspector" class="tab-content">
      <div class="inspector-grid">
        <div class="viewer-panel card">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
            <div class="card-title" style="margin-bottom: 0;">Short-distance scattering pattern fit (3.0 m geometry)</div>
            <span class="badge" id="sd_frame_badge">Frame 00 (0.060 s)</span>
          </div>
          <img id="sd_img" src="../short_distance/fits/frame_000.png" alt="Short Distance Fit">
        </div>

        <div class="card">
          <div class="card-title">Frame controller & telemetry</div>
          
          <div class="slider-box">
            <div style="display: flex; justify-content: space-between; font-size: 0.85rem;">
              <span>Frame: <strong id="sd_frame_val">0</strong> / 49</span>
              <span>Time: <strong id="sd_time_val">0.060 s</strong></span>
            </div>
            <input type="range" id="sd_slider" min="0" max="49" value="0" oninput="updateShortFrame(this.value)">
          </div>

          <div class="btn-group">
            <button class="btn" onclick="stepShort(-1)">&larr; Previous</button>
            <button class="btn" id="sd_play_btn" onclick="toggleShortPlay()">Play &para;</button>
            <button class="btn" onclick="stepShort(1)">Next &rarr;</button>
          </div>

          <table class="telemetry-table">
            <tr><th>Goodness of Fit &chi;<sup>2</sup><sub>red</sub></th><td id="sd_chi2">-</td></tr>
            <tr><th>Primary Radius &lang;R<sub>prim</sub>&rang;</th><td id="sd_Rmean">-</td></tr>
            <tr><th>Primary Dispersion &sigma;<sub>prim</sub></th><td id="sd_stdDev">-</td></tr>
            <tr><th>Fractal Dimension D<sub>f</sub></th><td id="sd_D">-</td></tr>
            <tr><th>Primary Fraction &phi;<sub>prim</sub> &times; 10<sup>4</sup></th><td id="sd_phi">-</td></tr>
            <tr><th>Fluct. Amplitude A<sub>fluc</sub></th><td id="sd_A_fluc">-</td></tr>
            <tr><th>Fluct. Mesh &xi;<sub>fluc</sub></th><td id="sd_xi_fluc">-</td></tr>
            <tr><th>Flat Background</th><td id="sd_bkg">-</td></tr>
          </table>

          <div style="font-size: 0.75rem; color: var(--text-secondary); line-height: 1.4;">
            Point-to-point evolutionary model with viscoelastic concentration fluctuations (&mu; = 4.0).
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 4: LONG DISTANCE INSPECTOR -->
    <div id="long_inspector" class="tab-content">
      <div class="inspector-grid">
        <div class="viewer-panel card">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
            <div class="card-title" style="margin-bottom: 0;">Long-distance scattering pattern fit (31.0 m geometry)</div>
            <span class="badge" id="ld_frame_badge">Frame 00 (0.060 s)</span>
          </div>
          <img id="ld_img" src="../long_distance/fits/frame_000.png" alt="Long Distance Fit">
        </div>

        <div class="card">
          <div class="card-title">Frame controller & telemetry</div>
          
          <div class="slider-box">
            <div style="display: flex; justify-content: space-between; font-size: 0.85rem;">
              <span>Frame: <strong id="ld_frame_val">0</strong> / 49</span>
              <span>Time: <strong id="ld_time_val">0.060 s</strong></span>
            </div>
            <input type="range" id="ld_slider" min="0" max="49" value="0" oninput="updateLongFrame(this.value)">
          </div>

          <div class="btn-group">
            <button class="btn" onclick="stepLong(-1)">&larr; Previous</button>
            <button class="btn" id="ld_play_btn" onclick="toggleLongPlay()">Play &para;</button>
            <button class="btn" onclick="stepLong(1)">Next &rarr;</button>
          </div>

          <table class="telemetry-table">
            <tr><th>Goodness of Fit &chi;<sup>2</sup><sub>red</sub></th><td id="ld_chi2">-</td></tr>
            <tr><th>Primary Droplets &lang;R<sub>1</sub>&rang;</th><td id="ld_Rlarge">-</td></tr>
            <tr><th>Primary Fraction &phi;<sub>1</sub> &times; 10<sup>5</sup></th><td id="ld_philarge">-</td></tr>
            <tr><th>Coalesced Droplets &lang;R<sub>2</sub>&rang;</th><td id="ld_Rlarge2">-</td></tr>
            <tr><th>Coalesced Fraction &phi;<sub>2</sub> &times; 10<sup>5</sup></th><td id="ld_philarge2">-</td></tr>
            <tr><th>Fractal Dimension D<sub>f</sub></th><td id="ld_D">-</td></tr>
            <tr><th>Network Cutoff &xi;</th><td id="ld_xi">-</td></tr>
            <tr><th>Primary Fraction &phi;<sub>prim</sub> &times; 10<sup>4</sup></th><td id="ld_phi">-</td></tr>
          </table>

          <div style="font-size: 0.75rem; color: var(--text-secondary); line-height: 1.4;">
            Transferred primary building block core: &lang;R<sub>prim</sub>&rang; = {bp['Rmean']:.2f} nm, &sigma;<sub>prim</sub> = {bp['stdDev']:.2f} nm.
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 5: PARAMETER TABLES -->
    <div id="tables" class="tab-content">
      <div class="card">
        <div class="card-title">Short-distance fitted parameters (3.0 m, 50 frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>&lang;R<sub>prim</sub>&rang; [nm]</th>
                <th>&sigma;<sub>prim</sub> [nm]</th>
                <th>D<sub>f</sub></th>
                <th>&phi;<sub>prim</sub> &times; 10<sup>4</sup></th>
                <th>A<sub>fluc</sub> [cm<sup>-1</sup>]</th>
                <th>&xi;<sub>fluc</sub> [nm]</th>
                <th>Bkg [cm<sup>-1</sup>]</th>
              </tr>
            </thead>
            <tbody>
"""

    for _, r in df_sd.iterrows():
        f = int(r['frame_index'])
        html_content += f"""              <tr>
                <td>{f:02d}</td>
                <td>{float(r['time_seconds']):.3f}</td>
                <td>{float(r['chi2_red']):.3f}</td>
                <td>{float(r['Rmean']):.2f} &plusmn; {float(r.get('std_err_Rmean', 0.0)):.2f}</td>
                <td>{float(r['stdDev']):.2f} &plusmn; {float(r.get('std_err_stdDev', 0.0)):.2f}</td>
                <td>{float(r['D']):.3f} &plusmn; {float(r.get('std_err_D', 0.0)):.3f}</td>
                <td>{float(r['phi']) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi', 0.0)) * 1e4:.2f}</td>
                <td>{float(r['A_fluc']):.4f} &plusmn; {float(r.get('std_err_A_fluc', 0.0)):.4f}</td>
                <td>{float(r['xi_fluc']):.2f} &plusmn; {float(r.get('std_err_xi_fluc', 0.0)):.2f}</td>
                <td>{float(r['bkg']):.4f} &plusmn; {float(r.get('std_err_bkg', 0.0)):.4f}</td>
              </tr>
"""

    html_content += """            </tbody>
          </table>
        </div>
      </div>

      <div class="card" style="margin-top: 1.5rem;">
        <div class="card-title">Long-distance fitted parameters (31.0 m, 50 frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>&lang;R<sub>1</sub>&rang; [nm]</th>
                <th>&phi;<sub>1</sub> &times; 10<sup>5</sup></th>
                <th>&lang;R<sub>2</sub>&rang; [nm]</th>
                <th>&phi;<sub>2</sub> &times; 10<sup>5</sup></th>
                <th>D<sub>f</sub></th>
                <th>&xi; [nm]</th>
                <th>&phi;<sub>prim</sub> &times; 10<sup>4</sup></th>
              </tr>
            </thead>
            <tbody>
"""

    for _, r in df_ld.iterrows():
        f = int(r['frame_index'])
        html_content += f"""              <tr>
                <td>{f:02d}</td>
                <td>{float(r['time_seconds']):.3f}</td>
                <td>{float(r['chi2_red']):.3f}</td>
                <td>{float(r['mean_radius_large']):.0f} &plusmn; {float(r.get('std_err_mean_radius_large', 0.0)):.0f}</td>
                <td>{float(r['total_volume_fraction_large']) * 1e5:.2f}</td>
                <td>{float(r.get('mean_radius_large_2', 0.0)):.0f} &plusmn; {float(r.get('std_err_mean_radius_large_2', 0.0)):.0f}</td>
                <td>{float(r.get('total_volume_fraction_large_2', 0.0)) * 1e5:.2f}</td>
                <td>{float(r['D']):.3f} &plusmn; {float(r.get('std_err_D', 0.0)):.3f}</td>
                <td>{float(r['xi']):.0f} &plusmn; {float(r.get('std_err_xi', 0.0)):.0f}</td>
                <td>{float(r['phi']) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi', 0.0)) * 1e4:.2f}</td>
              </tr>
"""

    html_content += r"""            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 6: THEORY & FORMULATIONS -->
    <div id="theory" class="tab-content">
      <div class="card">
        <div class="card-title">Multiscale physical scattering formalism</div>
        <p style="color: var(--text-secondary); margin-bottom: 12px;">
          The total experimental scattering intensity \(I(q)\) is modelled across the entire momentum transfer range
          \(q \in [0.002, 2.5]\text{ nm}^{-1}\) as the coherent linear superposition of four structural contributions:
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{tot}}(q) = I_{\mathrm{prim}}(q) + I_{\mathrm{drop}}(q) + I_{\mathrm{fluc}}(q) + I_{\mathrm{bkg}}
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">1. Primary fractal network</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Polydisperse spherical building blocks following a log-normal size distribution \(n(R)\) with inter-particle
          Teixeira mass fractal structure factor \(S(q)\):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{prim}}(q) = (\Delta\rho)^2 \int_0^\infty n(R) V^2(R) P(q, R) S(q, R, D_f, \xi) \, dR
          \]
          \[
          S(q, R, D_f, \xi) = 1 + \frac{D_f \Gamma(D_f - 1)}{(q R)^{D_f}} \frac{\sin[(D_f - 1) \arctan(q \xi)]}{\left(1 + (q \xi)^{-2}\right)^{(D_f - 1)/2}}
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">2. Secondary bimodal macroscopic droplets</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Under spinodal liquid-liquid phase separation in capillary flow, droplets exhibit a bimodal size distribution
          comprising primary dense liquid droplets (\(R_1 \sim 1.46\ \mu\mathrm{m}\)) coexisting with shear-coalesced giant droplets (\(R_2 \sim 5.03\ \mu\mathrm{m}\)):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{drop}}(q) = I_{\mathrm{drop}, 1}(q) + I_{\mathrm{drop}, 2}(q)
          \]
          \[
          I_{\mathrm{drop}, i}(q) = (\Delta\rho)^2 \int_0^\infty n_i(R) V^2(R) P(q, R) \, dR \quad (i = 1, 2)
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">3. Viscoelastic network concentration fluctuations</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Dynamic viscoelastic network fluctuations (Fisher-Burford / Tanaka formulation with \(\mu_{\mathrm{fluc}} = 4.0\)) explaining the high-\(q\) excess:
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{fluc}}(q) = \frac{A_{\mathrm{fluc}}}{\left[1 + \frac{2}{3 \mu_{\mathrm{fluc}}} (q \xi_{\mathrm{fluc}})^2\right]^{\mu_{\mathrm{fluc}} / 2}}
          \]
        </div>
      </div>
    </div>
  </main>

  <script>
    const shortData = """ + json.dumps(short_telemetry) + """;
    const longData = """ + json.dumps(long_telemetry) + """;

    let sdTimer = null;
    let ldTimer = null;

    function switchTab(tabId) {
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
      event.target.classList.add('active');
      document.getElementById(tabId).classList.add('active');
      if (tabId === 'short_inspector') updateShortFrame(document.getElementById('sd_slider').value);
      if (tabId === 'long_inspector') updateLongFrame(document.getElementById('ld_slider').value);
    }

    function padFrame(f) {
      return String(f).padStart(3, '0');
    }

    function updateShortFrame(f) {
      f = parseInt(f);
      const d = shortData[f];
      if (!d) return;

      document.getElementById('sd_frame_val').innerText = f;
      document.getElementById('sd_time_val').innerText = d.time + ' s';
      document.getElementById('sd_frame_badge').innerText = 'Frame ' + String(f).padStart(2, '0') + ' (' + d.time + ' s)';
      document.getElementById('sd_img').src = '../short_distance/fits/frame_' + padFrame(f) + '.png';

      document.getElementById('sd_chi2').innerText = d.chi2;
      document.getElementById('sd_Rmean').innerText = d.Rmean + ' nm';
      document.getElementById('sd_stdDev').innerText = d.stdDev + ' nm';
      document.getElementById('sd_D').innerText = d.D;
      document.getElementById('sd_phi').innerText = d.phi;
      document.getElementById('sd_A_fluc').innerText = d.A_fluc;
      document.getElementById('sd_xi_fluc').innerText = d.xi_fluc + ' nm';
      document.getElementById('sd_bkg').innerText = d.bkg;
    }

    function stepShort(delta) {
      const s = document.getElementById('sd_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 49) val = 49;
      s.value = val;
      updateShortFrame(val);
    }

    function toggleShortPlay() {
      const btn = document.getElementById('sd_play_btn');
      if (sdTimer) {
        clearInterval(sdTimer);
        sdTimer = null;
        btn.innerText = 'Play ¶';
      } else {
        btn.innerText = 'Pause ❚❚';
        sdTimer = setInterval(() => {
          const s = document.getElementById('sd_slider');
          let val = (parseInt(s.value) + 1) % 50;
          s.value = val;
          updateShortFrame(val);
        }, 150);
      }
    }

    function updateLongFrame(f) {
      f = parseInt(f);
      const d = longData[f];
      if (!d) return;

      document.getElementById('ld_frame_val').innerText = f;
      document.getElementById('ld_time_val').innerText = d.time + ' s';
      document.getElementById('ld_frame_badge').innerText = 'Frame ' + String(f).padStart(2, '0') + ' (' + d.time + ' s)';
      document.getElementById('ld_img').src = '../long_distance/fits/frame_' + padFrame(f) + '.png';

      document.getElementById('ld_chi2').innerText = d.chi2;
      document.getElementById('ld_Rlarge').innerText = d.mean_radius_large + ' nm';
      document.getElementById('ld_philarge').innerText = d.total_volume_fraction_large;
      document.getElementById('ld_Rlarge2').innerText = (d.mean_radius_large_2 || '-') + ' nm';
      document.getElementById('ld_philarge2').innerText = d.total_volume_fraction_large_2 || '-';
      document.getElementById('ld_D').innerText = d.D;
      document.getElementById('ld_xi').innerText = d.xi + ' nm';
      document.getElementById('ld_phi').innerText = d.phi;
    }

    function stepLong(delta) {
      const s = document.getElementById('ld_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 49) val = 49;
      s.value = val;
      updateLongFrame(val);
    }

    function toggleLongPlay() {
      const btn = document.getElementById('ld_play_btn');
      if (ldTimer) {
        clearInterval(ldTimer);
        ldTimer = null;
        btn.innerText = 'Play ¶';
      } else {
        btn.innerText = 'Pause ❚❚';
        ldTimer = setInterval(() => {
          const s = document.getElementById('ld_slider');
          let val = (parseInt(s.value) + 1) % 50;
          s.value = val;
          updateLongFrame(val);
        }, 150);
      }
    }

    // Initialize telemetry
    updateShortFrame(0);
    updateLongFrame(0);
  </script>
</body>
</html>
"""

    dash_path = os.path.join(dashboard_dir, 'index.html')
    with open(dash_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Master interactive dashboard compiled: {dash_path}")
    print("=" * 80)


if __name__ == '__main__':
    build_dashboard()
