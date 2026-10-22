#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-005.
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

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
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
  <title>ND-i-005 Multiscale Kinetic Analysis Dashboard</title>
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
      padding-bottom: 60px;
    }}
    header {{
      background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
      border-bottom: 1px solid var(--border);
      padding: 24px 32px;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 16px;
    }}
    .header-title h1 {{
      font-size: 1.6rem;
      font-weight: 700;
      color: #ffffff;
      letter-spacing: -0.02em;
    }}
    .header-title p {{
      font-size: 0.95rem;
      color: var(--text-secondary);
      margin-top: 4px;
    }}
    .badge {{
      display: inline-block;
      padding: 4px 10px;
      font-size: 0.8rem;
      font-weight: 600;
      border-radius: 9999px;
      background-color: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      border: 1px solid rgba(56, 189, 248, 0.3);
    }}
    .badge-success {{
      background-color: rgba(16, 185, 129, 0.15);
      color: var(--success);
      border: 1px solid rgba(16, 185, 129, 0.3);
    }}
    nav {{
      background-color: var(--bg-secondary);
      border-bottom: 1px solid var(--border);
      padding: 0 32px;
      display: flex;
      gap: 8px;
      overflow-x: auto;
    }}
    .tab-btn {{
      background: none;
      border: none;
      color: var(--text-secondary);
      padding: 14px 18px;
      font-size: 0.95rem;
      font-weight: 600;
      cursor: pointer;
      border-bottom: 2px solid transparent;
      transition: all 0.2s ease;
      white-space: nowrap;
    }}
    .tab-btn:hover {{
      color: var(--text-primary);
    }}
    .tab-btn.active {{
      color: var(--accent);
      border-bottom: 2px solid var(--accent);
    }}
    main {{
      max-width: 1440px;
      margin: 0 auto;
      padding: 32px;
    }}
    .tab-content {{
      display: none;
    }}
    .tab-content.active {{
      display: block;
      animation: fadeIn 0.3s ease;
    }}
    @keyframes fadeIn {{
      from {{ opacity: 0; transform: translateY(4px); }}
      to {{ opacity: 1; transform: translateY(0); }}
    }}
    .card {{
      background-color: var(--bg-card);
      border: 1px solid var(--card-border);
      border-radius: 12px;
      padding: 24px;
      margin-bottom: 24px;
      box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
    }}
    .card-title {{
      font-size: 1.2rem;
      font-weight: 600;
      margin-bottom: 16px;
      display: flex;
      align-items: center;
      gap: 10px;
      color: #ffffff;
    }}
    .grid-2 {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 24px;
    }}
    @media (max-width: 1024px) {{
      .grid-2 {{
        grid-template-columns: 1fr;
      }}
    }}
    .inspector-layout {{
      display: grid;
      grid-template-columns: 1fr 380px;
      gap: 24px;
    }}
    @media (max-width: 1100px) {{
      .inspector-layout {{
        grid-template-columns: 1fr;
      }}
    }}
    .viewer-box {{
      background-color: #0b1120;
      border: 1px solid var(--border);
      border-radius: 8px;
      overflow: hidden;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      padding: 12px;
    }}
    .viewer-img {{
      max-width: 100%;
      height: auto;
      border-radius: 4px;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.4);
    }}
    .controls-panel {{
      background-color: var(--bg-secondary);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 20px;
      display: flex;
      flex-direction: column;
      gap: 18px;
    }}
    .slider-row {{
      display: flex;
      flex-direction: column;
      gap: 8px;
    }}
    .slider-labels {{
      display: flex;
      justify-content: space-between;
      font-size: 0.85rem;
      color: var(--text-secondary);
    }}
    input[type=range] {{
      width: 100%;
      accent-color: var(--accent);
      cursor: pointer;
    }}
    .btn-group {{
      display: flex;
      gap: 8px;
    }}
    .btn {{
      flex: 1;
      padding: 8px 12px;
      border-radius: 6px;
      font-weight: 600;
      font-size: 0.85rem;
      cursor: pointer;
      border: 1px solid var(--border);
      background-color: #273549;
      color: var(--text-primary);
      transition: all 0.15s ease;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 6px;
    }}
    .btn:hover {{
      background-color: var(--accent);
      color: #0f172a;
      border-color: var(--accent);
    }}
    .telemetry-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.85rem;
    }}
    .telemetry-table th, .telemetry-table td {{
      padding: 8px 10px;
      text-align: left;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    .telemetry-table th {{
      color: var(--text-secondary);
      font-weight: 500;
      width: 50%;
    }}
    .telemetry-table td {{
      color: #ffffff;
      font-weight: 600;
      text-align: right;
      font-variant-numeric: tabular-nums;
    }}
    .data-table-container {{
      max-height: 520px;
      overflow-y: auto;
      border: 1px solid var(--border);
      border-radius: 8px;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
    }}
    .data-table th {{
      position: sticky;
      top: 0;
      background-color: #243044;
      color: var(--text-primary);
      font-weight: 600;
      padding: 10px 12px;
      text-align: right;
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
    }}
    .data-table th:first-child {{
      text-align: left;
    }}
    .data-table td {{
      padding: 7px 12px;
      text-align: right;
      border-bottom: 1px solid rgba(255, 255, 255, 0.04);
      color: var(--text-secondary);
      font-variant-numeric: tabular-nums;
    }}
    .data-table td:first-child {{
      text-align: left;
      color: var(--text-primary);
      font-weight: 600;
    }}
    .data-table tr:hover td {{
      background-color: rgba(56, 189, 248, 0.05);
      color: #ffffff;
    }}
    .equation-box {{
      background-color: #0b1120;
      border-left: 3px solid var(--accent);
      border-radius: 4px;
      padding: 16px 20px;
      margin: 16px 0;
      overflow-x: auto;
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-title">
      <h1>ND-i-005 Multiscale Kinetic Analysis Dashboard</h1>
      <p>40 mM CaCO<sub>3</sub>, pH 11.3 with Cobalt Additive &bull; High-Speed Time-Resolved SAXS (0.060 s to 10.877 s)</p>
    </div>
    <div style="display: flex; gap: 10px; align-items: center;">
      <a href="../../../index.html" class="badge" style="text-decoration: none; cursor: pointer;">&larr; Master Dashboard</a>
      <span class="badge badge-success">50 Time Frames</span>
      <span class="badge">q: 0.002 &ndash; 2.5 nm<sup>-1</sup></span>
      <span class="badge">Uncertainties: &plusmn;1&sigma;</span>
    </div>
  </header>

  <nav>
    <button class="tab-btn active" onclick="switchTab('overview')">1. Overview & Composite Curve</button>
    <button class="tab-btn" onclick="switchTab('cross_scale')">2. Cross-Scale Kinetics</button>
    <button class="tab-btn" onclick="switchTab('short_inspector')">3. Short Distance (3.0 m)</button>
    <button class="tab-btn" onclick="switchTab('long_inspector')">4. Long Distance (31.0 m)</button>
    <button class="tab-btn" onclick="switchTab('tables')">5. Parameter Tables</button>
    <button class="tab-btn" onclick="switchTab('theory')">6. Physical Formalism</button>
  </nav>

  <main>
    <!-- TAB 1: OVERVIEW & BASELINE -->
    <div id="overview" class="tab-content active">
      <div class="grid-2">
        <div class="card">
          <div class="card-title">Full Multiscale Composite Baseline Curve</div>
          <p style="color: var(--text-secondary); margin-bottom: 16px; font-size: 0.9rem;">
            Stitched wide-angle scattering curve combining the final reaction stages across both cameras (q: 0.002 to 2.5 nm<sup>-1</sup>).
            Resolves primary amorphous building blocks, the mass fractal network, secondary droplets, and sub-nanometer clusters.
          </p>
          <div class="viewer-box">
            <img class="viewer-img" src="../composite/composite_baseline_fit.png" alt="Composite Baseline Fit">
          </div>
        </div>

        <div class="card">
          <div class="card-title">Converged Baseline Structural Parameters (&plusmn;1&sigma;)</div>
          <table class="telemetry-table" style="margin-top: 12px;">
            <tr><th>Primary Core Radius &lang;R<sub>prim</sub>&rang;</th><td>{bp['Rmean']:.2f} &plusmn; {bp.get('std_err_Rmean', 0.0):.2f} nm</td></tr>
            <tr><th>Primary Core Dispersion &sigma;<sub>prim</sub></th><td>{bp['stdDev']:.2f} &plusmn; {bp.get('std_err_stdDev', 0.0):.2f} nm</td></tr>
            <tr><th>Mass Fractal Dimension D<sub>f</sub></th><td>{bp['D']:.3f} &plusmn; {bp.get('std_err_D', 0.0):.3f}</td></tr>
            <tr><th>Network Correlation Cutoff &xi;</th><td>{bp['xi']:.0f} &plusmn; {bp.get('std_err_xi', 0.0):.0f} nm</td></tr>
            <tr><th>Secondary Droplet Radius &lang;R<sub>large</sub>&rang;</th><td>{bp['mean_radius_large']:.0f} &plusmn; {bp.get('std_err_mean_radius_large', 0.0):.0f} nm</td></tr>
            <tr><th>Secondary Droplet Dispersion &sigma;<sub>large</sub></th><td>{bp['std_radius_large']:.0f} &plusmn; {bp.get('std_err_std_radius_large', 0.0):.0f} nm</td></tr>
            <tr><th>Droplet Volume Fraction &phi;<sub>large</sub></th><td>{bp['total_volume_fraction_large'] * 1e5:.2f} &plusmn; {bp.get('std_err_total_volume_fraction_large', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td></tr>
            <tr><th>Sub-Nanometer Cluster Radius R<sub>pnc</sub></th><td>{bp['R_pnc']:.2f} &plusmn; {bp.get('std_err_R_pnc', 0.0):.2f} nm</td></tr>
            <tr><th>Cluster Population Fraction &phi;<sub>pnc</sub></th><td>{bp['phi_pnc'] * 1e4:.2f} &plusmn; {bp.get('std_err_phi_pnc', 0.0) * 1e4:.2f} &times; 10<sup>-4</sup></td></tr>
            <tr><th>Incoherent Flat Background</th><td>{bp['bkg']:.6f} &plusmn; {bp.get('std_err_bkg', 0.0):.6f} cm<sup>-1</sup></td></tr>
            <tr><th>Overall Baseline Fit Quality &chi;<sup>2</sup><sub>red</sub></th><td><span class="badge badge-success">{bp['redchi']:.3f}</span></td></tr>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 2: CROSS-SCALE KINETICS -->
    <div id="cross_scale" class="tab-content">
      <div class="card">
        <div class="card-title">Unified Multi-Scale Kinetic Synthesis (Real Time 0.060 s to 10.877 s)</div>
        <p style="color: var(--text-secondary); margin-bottom: 16px; font-size: 0.9rem;">
          Cross-scale synchronization demonstrating consistent structural length scale hierarchy across 4 decades,
          rigorous agreement between short-distance and long-distance fractal dimensions (D<sub>f</sub> &plusmn; 1&sigma;), and
          coupled volume fraction evolution.
        </p>
        <div class="viewer-box">
          <img class="viewer-img" src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Cross-Scale Trajectories">
        </div>
      </div>
    </div>

    <!-- TAB 3: SHORT DISTANCE INSPECTOR -->
    <div id="short_inspector" class="tab-content">
      <div class="inspector-layout">
        <div class="card" style="margin-bottom: 0;">
          <div class="card-title">
            <span>Short-Distance Frame Fit (3.0 m Camera)</span>
            <span id="sd_frame_badge" class="badge">Frame 00 (0.060 s)</span>
          </div>
          <div class="viewer-box">
            <img id="sd_img" class="viewer-img" src="../short_distance/fits/frame_000.png" alt="Short Distance Fit">
          </div>
        </div>

        <div class="controls-panel">
          <div class="slider-row">
            <div class="slider-labels">
              <span>Time Frame: <b id="sd_frame_val" style="color: var(--accent);">0</b> / 49</span>
              <span>Reaction Time: <b id="sd_time_val" style="color: #ffffff;">0.060 s</b></span>
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
            <tr><th>Primary Core &lang;R<sub>prim</sub>&rang;</th><td id="sd_Rmean">-</td></tr>
            <tr><th>Core Dispersion &sigma;<sub>prim</sub></th><td id="sd_stdDev">-</td></tr>
            <tr><th>Fractal Dimension D<sub>f</sub></th><td id="sd_D">-</td></tr>
            <tr><th>Primary Fraction &phi;<sub>prim</sub> &times; 10<sup>4</sup></th><td id="sd_phi">-</td></tr>
            <tr><th>Viscoelastic Fluct. Ampl. A<sub>fluc</sub></th><td id="sd_A_fluc">-</td></tr>
            <tr><th>Fluctuation Mesh &xi;<sub>fluc</sub></th><td id="sd_xi_fluc">-</td></tr>
            <tr><th>Incoherent Bkg [cm<sup>-1</sup>]</th><td id="sd_bkg">-</td></tr>
          </table>

          <div style="font-size: 0.75rem; color: var(--text-secondary); line-height: 1.4;">
            Use slider or left/right keyboard arrow keys to navigate frames.
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 4: LONG DISTANCE INSPECTOR -->
    <div id="long_inspector" class="tab-content">
      <div class="inspector-layout">
        <div class="card" style="margin-bottom: 0;">
          <div class="card-title">
            <span>Long-Distance Frame Fit (31.0 m Camera)</span>
            <span id="ld_frame_badge" class="badge">Frame 00 (0.060 s)</span>
          </div>
          <div class="viewer-box">
            <img id="ld_img" class="viewer-img" src="../long_distance/fits/frame_000.png" alt="Long Distance Fit">
          </div>
        </div>

        <div class="controls-panel">
          <div class="slider-row">
            <div class="slider-labels">
              <span>Time Frame: <b id="ld_frame_val" style="color: var(--accent);">0</b> / 49</span>
              <span>Reaction Time: <b id="ld_time_val" style="color: #ffffff;">0.060 s</b></span>
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
            <tr><th>Droplet Radius &lang;R<sub>large</sub>&rang;</th><td id="ld_Rlarge">-</td></tr>
            <tr><th>Droplet Dispersion &sigma;<sub>large</sub></th><td id="ld_stdLarge">-</td></tr>
            <tr><th>Droplet Fraction &phi;<sub>large</sub> &times; 10<sup>5</sup></th><td id="ld_philarge">-</td></tr>
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
        <div class="card-title">Short-Distance Fitted Parameters (3.0 m, 50 Frames)</div>
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
        t_sec = float(r['time_seconds'])
        html_content += f"""              <tr>
                <td>Frame {f:02d}</td>
                <td>{t_sec:.3f}</td>
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

      <div class="card">
        <div class="card-title">Long-Distance Fitted Parameters (31.0 m, 50 Frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>&lang;R<sub>large</sub>&rang; [nm]</th>
                <th>&sigma;<sub>large</sub> [nm]</th>
                <th>&phi;<sub>large</sub> &times; 10<sup>5</sup></th>
                <th>D<sub>f</sub></th>
                <th>&xi; [nm]</th>
                <th>&phi;<sub>prim</sub> &times; 10<sup>4</sup></th>
              </tr>
            </thead>
            <tbody>
"""

    for _, r in df_ld.iterrows():
        f = int(r['frame_index'])
        t_sec = float(r['time_seconds'])
        html_content += f"""              <tr>
                <td>Frame {f:02d}</td>
                <td>{t_sec:.3f}</td>
                <td>{float(r['chi2_red']):.3f}</td>
                <td>{float(r['mean_radius_large']):.0f} &plusmn; {float(r.get('std_err_mean_radius_large', 0.0)):.0f}</td>
                <td>{float(r['std_radius_large']):.0f} &plusmn; {float(r.get('std_err_std_radius_large', 0.0)):.0f}</td>
                <td>{float(r['total_volume_fraction_large']) * 1e5:.2f} &plusmn; {float(r.get('std_err_total_volume_fraction_large', 0.0)) * 1e5:.2f}</td>
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
        <div class="card-title">Multiscale Physical Scattering Formalism</div>
        <p style="color: var(--text-secondary); margin-bottom: 12px;">
          The total experimental scattering intensity \(I(q)\) is modelled across the entire momentum transfer range
          \(q \in [0.002, 2.5]\text{ nm}^{-1}\) as the coherent linear superposition of four structural contributions:
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{tot}}(q) = I_{\mathrm{prim}}(q) + I_{\mathrm{large}}(q) + I_{\mathrm{fluc}}(q) + I_{\mathrm{bkg}}
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">1. Primary Fractal Network</h3>
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

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">2. Secondary Macroscopic Droplets</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Dilute secondary dense liquid phase droplets with log-normal polydispersity and constant relative dispersion
          \(\sigma_{\mathrm{large}} / \langle R_{\mathrm{large}} \rangle = 0.143\):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{large}}(q) = (\Delta\rho)^2 \int_0^\infty n_{\mathrm{large}}(R) V^2(R) P(q, R) \, dR
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">3. Viscoelastic Network Concentration Fluctuations</h3>
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
      document.getElementById('ld_stdLarge').innerText = d.std_radius_large + ' nm';
      document.getElementById('ld_philarge').innerText = d.total_volume_fraction_large;
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

    // Keyboard navigation
    window.addEventListener('keydown', (e) => {
      const activeTab = document.querySelector('.tab-content.active').id;
      if (activeTab === 'short_inspector') {
        if (e.key === 'ArrowLeft') stepShort(-1);
        if (e.key === 'ArrowRight') stepShort(1);
      } else if (activeTab === 'long_inspector') {
        if (e.key === 'ArrowLeft') stepLong(-1);
        if (e.key === 'ArrowRight') stepLong(1);
      }
    });

    // Initialize frame 0
    updateShortFrame(0);
    updateLongFrame(0);
  </script>
</body>
</html>
"""

    dashboard_file = os.path.join(dashboard_dir, 'index.html')
    with open(dashboard_file, 'w') as f_out:
        f_out.write(html_content)
    print(f"Saved master interactive dashboard: {dashboard_file}")


if __name__ == '__main__':
    build_dashboard()
