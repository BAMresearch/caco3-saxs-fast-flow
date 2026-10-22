#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-006.
-----------------------------------------------------------------------
Generates a standalone, responsive HTML dashboard for reviewing:
1. Overview & Composite Anchor Fit (q in [0.002, 2.5] nm^-1)
2. Unified Cross-Scale Kinetic Trajectories (0.060 s to 10.877 s)
3. Short-Distance Time-Resolved Inspector (3.0 m, Frames 0 to 49, with live slider & +/- 1 sigma telemetry)
4. Long-Distance Time-Resolved Inspector (31.0 m, Frames 0 to 49, with live slider & +/- 1 sigma telemetry)
5. Comprehensive Parameter Tables (searchable & scrollable)
6. Physical Model Formulations & Equations rendered via MathJax.

System: 40.0 mM CaCO3, pH 9.75 (Spinodal Mass Fractal Network & Large Droplets, No ACC Particles)

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
            'Dm': f"{float(r['Dm']):.3f} ± {float(r.get('std_err_Dm', 0.0)):.3f}",
            'phi_net': f"{float(r.get('phi_net', 0.0)) * 1e4:.2f} ± {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}",
            'A_porod': f"{float(r.get('A_porod', 0.0)) * 1e6:.2f} ± {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}",
            'R_prim': f"{float(r.get('R_prim', 1.50)):.2f}",
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
            'Dm': f"{float(r['Dm']):.3f} ± {float(r.get('std_err_Dm', 0.0)):.3f}",
            'xi': f"{float(r['xi']):.0f} ± {float(r.get('std_err_xi', 0.0)):.0f}",
            'phi_net': f"{float(r['phi_net']) * 1e4:.2f} ± {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}",
            'R_large': f"{float(r.get('mean_radius_large', 3113.0)):.0f} ± {float(r.get('std_err_mean_radius_large', 0.0)):.0f}",
            'phi_large': f"{float(r.get('total_volume_fraction_large', 0.0)) * 1e4:.2f} ± {float(r.get('std_err_total_volume_fraction_large', 0.0)) * 1e4:.2f}",
            'A_porod': f"{float(r.get('A_porod', 0.0)) * 1e6:.2f} ± {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}",
            'bkg': f"{float(r['bkg']):.4f} ± {float(r.get('std_err_bkg', 0.0)):.4f}"
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
  <title>ND-i-006 Kinetic SAXS Suite (Fractal Network & Large Spheres)</title>
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
      background: var(--bg-primary);
      color: var(--text-primary);
      line-height: 1.5;
      padding-bottom: 3rem;
    }}
    header {{
      background: var(--bg-secondary);
      border-bottom: 1px solid var(--border);
      padding: 1.25rem 2rem;
      position: sticky;
      top: 0;
      z-index: 100;
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
    .title-group h1 {{
      font-size: 1.4rem;
      font-weight: 700;
      color: var(--text-primary);
    }}
    .title-group p {{
      font-size: 0.85rem;
      color: var(--text-secondary);
    }}
    .badge {{
      display: inline-block;
      padding: 0.25rem 0.6rem;
      font-size: 0.75rem;
      font-weight: 600;
      border-radius: 9999px;
      background: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      border: 1px solid rgba(56, 189, 248, 0.3);
    }}
    .nav-tabs {{
      max-width: 1400px;
      margin: 1.5rem auto 1rem auto;
      padding: 0 1.5rem;
      display: flex;
      gap: 0.5rem;
      border-bottom: 1px solid var(--border);
      overflow-x: auto;
    }}
    .tab-btn {{
      background: transparent;
      border: none;
      color: var(--text-secondary);
      font-size: 0.9rem;
      font-weight: 600;
      padding: 0.75rem 1.25rem;
      cursor: pointer;
      border-bottom: 2px solid transparent;
      transition: all 0.2s;
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
      max-width: 1400px;
      margin: 0 auto;
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
      border-radius: 8px;
      padding: 1.5rem;
      margin-bottom: 1.5rem;
      box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
    }}
    .card-title {{
      font-size: 1.1rem;
      font-weight: 700;
      margin-bottom: 1rem;
      color: var(--text-primary);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    .inspector-grid {{
      display: grid;
      grid-template-columns: 2fr 1fr;
      gap: 1.5rem;
    }}
    @media (max-width: 900px) {{
      .inspector-grid {{
        grid-template-columns: 1fr;
      }}
    }}
    .viewer-panel {{
      background: rgba(0,0,0,0.25);
      border-radius: 8px;
      padding: 1rem;
      border: 1px solid var(--border);
    }}
    .viewer-panel img {{
      width: 100%;
      height: auto;
      border-radius: 6px;
      display: block;
      background: #fff;
    }}
    .slider-box {{
      margin: 1rem 0;
      display: flex;
      flex-direction: column;
      gap: 0.5rem;
    }}
    .slider-box input[type=range] {{
      width: 100%;
      accent-color: var(--accent);
      cursor: pointer;
    }}
    .btn-group {{
      display: flex;
      gap: 0.5rem;
      justify-content: center;
    }}
    .btn {{
      background: var(--bg-secondary);
      border: 1px solid var(--border);
      color: var(--text-primary);
      padding: 0.4rem 0.9rem;
      border-radius: 4px;
      font-size: 0.85rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.2s;
    }}
    .btn:hover {{
      background: var(--border);
    }}
    .telemetry-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.85rem;
    }}
    .telemetry-table th, .telemetry-table td {{
      padding: 0.5rem 0.75rem;
      text-align: left;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    .telemetry-table th {{
      color: var(--text-secondary);
      font-weight: 500;
      width: 50%;
    }}
    .telemetry-table td {{
      font-weight: 600;
      color: var(--accent);
      text-align: right;
    }}
    .data-table-container {{
      max-height: 480px;
      overflow-y: auto;
      border: 1px solid var(--border);
      border-radius: 6px;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
    }}
    .data-table th, .data-table td {{
      padding: 0.5rem 0.75rem;
      text-align: right;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    .data-table th {{
      background: var(--bg-secondary);
      color: var(--text-secondary);
      position: sticky;
      top: 0;
      font-weight: 600;
    }}
    .data-table tr:hover {{
      background: rgba(255, 255, 255, 0.03);
    }}
    .equation-box {{
      background: rgba(0, 0, 0, 0.3);
      padding: 1rem;
      border-radius: 6px;
      border-left: 3px solid var(--accent);
      margin: 1rem 0;
      overflow-x: auto;
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-content">
      <div class="title-group">
        <h1>ND-i-006 Kinetic SAXS Suite</h1>
        <p>40.0 mM CaCO<sub>3</sub>, pH 9.75 &mdash; Spinodal Mass Fractal Network &amp; Large Droplets (No ACC Particles)</p>
      </div>
      <div style="display: flex; gap: 0.5rem; align-items: center;">
        <a href="../../../index.html" class="badge" style="text-decoration: none; cursor: pointer;">&larr; Master Dashboard</a>
        <span class="badge" style="background: rgba(16, 185, 129, 0.2); color: #10b981;">No ACC Nanospheres</span>
        <span class="badge">50 Time-Resolved Frames</span>
      </div>
    </div>
  </header>

  <nav class="nav-tabs">
    <button class="tab-btn active" onclick="switchTab('overview')">1. Overview & Composite Curve</button>
    <button class="tab-btn" onclick="switchTab('unified_trajectories')">2. Cross-Scale Kinetics</button>
    <button class="tab-btn" onclick="switchTab('short_inspector')">3. Short Distance (3.0 m)</button>
    <button class="tab-btn" onclick="switchTab('long_inspector')">4. Long Distance (31.0 m)</button>
    <button class="tab-btn" onclick="switchTab('tables')">5. Parameter Tables</button>
    <button class="tab-btn" onclick="switchTab('theory')">6. Physical Formalism</button>
  </nav>

  <main>
    <!-- TAB 1: OVERVIEW & BASELINE -->
    <div id="overview" class="tab-content active">
      <div class="card">
        <div class="card-title">Experimental Series Overview & Scientific Rationale</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem;">
          At 40.0 mM CaCO<sub>3</sub> and pH 9.75, precipitation proceeds via <strong>liquid-liquid spinodal decomposition</strong>
          forming an interconnected percolated mass fractal network along with secondary macroscopic large droplets/spheres.
          Experimental SAXS profiles confirm the <strong>absence of discrete ACC nanospheres</strong> (no intermediate ~30–50 nm particle knee).
          The complete multi-scale composite scattering curve spans four decades in momentum transfer (<em>q</em> &isin; [0.002, 2.5] nm<sup>-1</sup>).
        </p>
        <img src="../composite/composite_baseline_fit.png" alt="Composite Baseline Fit" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>

      <div class="card">
        <div class="card-title">Converged Structural Baseline Parameters (&plusmn;1&sigma;)</div>
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
              <td>Primary Building Block Radius</td><td>R<sub>prim</sub></td>
              <td>{bp['R_prim']:.2f} nm</td><td>&plusmn; {bp.get('std_err_R_prim', 0.0):.2f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of amorphous building blocks</td>
            </tr>
            <tr>
              <td>Mass Fractal Dimension</td><td>D<sub>m</sub></td>
              <td>{bp['Dm']:.3f}</td><td>&plusmn; {bp.get('std_err_Dm', 0.0):.3f}</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Spinodal network space-filling scaling exponent</td>
            </tr>
            <tr>
              <td>Network Cutoff Length</td><td>&xi;</td>
              <td>{bp['xi']:.0f} nm</td><td>&plusmn; {bp.get('std_err_xi', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Percolation cluster correlation cutoff length</td>
            </tr>
            <tr>
              <td>Network Volume Fraction</td><td>&phi;<sub>net</sub></td>
              <td>{bp['phi_net'] * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_phi_net', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of the mesoscale percolated network</td>
            </tr>
            <tr>
              <td>Large Sphere / Droplet Radius</td><td>&lang;R<sub>large</sub>&rang;</td>
              <td>{bp.get('mean_radius_large', 3113.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_mean_radius_large', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of macroscopic phase-separated droplets</td>
            </tr>
            <tr>
              <td>Large Sphere Polydispersity</td><td>&sigma;<sub>large</sub></td>
              <td>{bp.get('std_radius_large', 267.0):.0f} nm</td><td>&plusmn; {bp.get('std_err_std_radius_large', 0.0):.0f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Polydispersity of secondary droplets</td>
            </tr>
            <tr>
              <td>Large Sphere Volume Fraction</td><td>&phi;<sub>large</sub></td>
              <td>{bp.get('total_volume_fraction_large', 1.1e-4) * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_total_volume_fraction_large', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of macroscopic dense liquid phase</td>
            </tr>
            <tr>
              <td>Porod Surface Constant</td><td>A<sub>porod</sub></td>
              <td>{bp['A_porod'] * 1e7:.2f} &times; 10<sup>-7</sup></td><td>&plusmn; {bp.get('std_err_A_porod', 0.0) * 1e7:.2f} &times; 10<sup>-7</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Macroscopic coarsened domain surface scattering amplitude</td>
            </tr>
            <tr>
              <td>Incoherent Background</td><td>bkg</td>
              <td>{bp['bkg']:.4f} cm<sup>-1</sup></td><td>&plusmn; {bp.get('std_err_bkg', 0.0):.4f} cm<sup>-1</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">High-q flat incoherent background</td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>

    <!-- TAB 2: UNIFIED TRAJECTORIES -->
    <div id="unified_trajectories" class="tab-content">
      <div class="card">
        <div class="card-title">Cross-Scale Kinetic Trajectories (0.060 s to 10.877 s)</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem;">
          Synthesised evolution across both camera geometries (Short Distance 3.0 m & Long Distance 31.0 m) illustrating
          spinodal network percolation and macroscopic droplet evolution (without intermediate ACC nanospheres).
        </p>
        <img src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Cross Scale Trajectories" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>
    </div>

    <!-- TAB 3: SHORT DISTANCE INSPECTOR -->
    <div id="short_inspector" class="tab-content">
      <div class="card">
        <div class="card-title">Short-Distance Time-Resolved Inspector (3.0 m, q &isin; [0.022, 2.54] nm<sup>-1</sup>)</div>
        <div class="inspector-grid">
          <div class="viewer-panel">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem;">
              <span id="sd_frame_badge" class="badge" style="font-size: 0.85rem;">Frame 00 (0.060 s)</span>
              <span style="font-size: 0.85rem; color: var(--text-secondary);">50 frames total</span>
            </div>
            <img id="sd_img" src="../short_distance/fits/frame_00.png" alt="Short Frame Fit">
            <div class="slider-box">
              <input type="range" id="sd_slider" min="0" max="49" value="0" oninput="updateShortFrame(this.value)">
              <div style="display: flex; justify-content: space-between; font-size: 0.75rem; color: var(--text-secondary);">
                <span>Frame 0 (0.060 s)</span>
                <span>Frame 25 (0.435 s)</span>
                <span>Frame 49 (10.877 s)</span>
              </div>
            </div>
            <div class="btn-group">
              <button class="btn" onclick="stepShort(-1)">&larr; Prev</button>
              <button class="btn" id="sd_play_btn" onclick="toggleShortPlay()">Play &#9658;</button>
              <button class="btn" onclick="stepShort(1)">Next &rarr;</button>
            </div>
          </div>
          <div class="telemetry-panel">
            <div style="background: rgba(0,0,0,0.2); padding: 1rem; border-radius: 6px; border: 1px solid var(--border);">
              <div style="font-weight: 600; font-size: 0.95rem; margin-bottom: 0.75rem; color: var(--accent);">Frame Parameters</div>
              <table class="telemetry-table">
                <tr><th>Frame Index</th><td id="sd_frame_val">0</td></tr>
                <tr><th>Reaction Time</th><td id="sd_time_val">0.060 s</td></tr>
                <tr><th>Reduced &chi;<sup>2</sup></th><td id="sd_chi2">-</td></tr>
                <tr><th>Fractal Dim D<sub>m</sub></th><td id="sd_Dm">-</td></tr>
                <tr><th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th><td id="sd_phi_net">-</td></tr>
                <tr><th>A<sub>porod</sub> &times; 10<sup>6</sup></th><td id="sd_A_porod">-</td></tr>
                <tr><th>R<sub>prim</sub> (fixed)</th><td id="sd_R_prim">1.50 nm</td></tr>
                <tr><th>Background</th><td id="sd_bkg">-</td></tr>
              </table>
            </div>
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 4: LONG DISTANCE INSPECTOR -->
    <div id="long_inspector" class="tab-content">
      <div class="card">
        <div class="card-title">Long-Distance Time-Resolved Inspector (31.0 m, q &isin; [0.002, 0.22] nm<sup>-1</sup>)</div>
        <div class="inspector-grid">
          <div class="viewer-panel">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem;">
              <span id="ld_frame_badge" class="badge" style="font-size: 0.85rem;">Frame 00 (0.060 s)</span>
              <span style="font-size: 0.85rem; color: var(--text-secondary);">50 frames total</span>
            </div>
            <img id="ld_img" src="../long_distance/fits/frame_00.png" alt="Long Frame Fit">
            <div class="slider-box">
              <input type="range" id="ld_slider" min="0" max="49" value="0" oninput="updateLongFrame(this.value)">
              <div style="display: flex; justify-content: space-between; font-size: 0.75rem; color: var(--text-secondary);">
                <span>Frame 0 (0.060 s)</span>
                <span>Frame 25 (0.435 s)</span>
                <span>Frame 49 (10.877 s)</span>
              </div>
            </div>
            <div class="btn-group">
              <button class="btn" onclick="stepLong(-1)">&larr; Prev</button>
              <button class="btn" id="ld_play_btn" onclick="toggleLongPlay()">Play &#9658;</button>
              <button class="btn" onclick="stepLong(1)">Next &rarr;</button>
            </div>
          </div>
          <div class="telemetry-panel">
            <div style="background: rgba(0,0,0,0.2); padding: 1rem; border-radius: 6px; border: 1px solid var(--border);">
              <div style="font-weight: 600; font-size: 0.95rem; margin-bottom: 0.75rem; color: var(--accent);">Frame Parameters</div>
              <table class="telemetry-table">
                <tr><th>Frame Index</th><td id="ld_frame_val">0</td></tr>
                <tr><th>Reaction Time</th><td id="ld_time_val">0.060 s</td></tr>
                <tr><th>Reduced &chi;<sup>2</sup></th><td id="ld_chi2">-</td></tr>
                <tr><th>Fractal Dim D<sub>m</sub></th><td id="ld_Dm">-</td></tr>
                <tr><th>Cutoff &xi;</th><td id="ld_xi">-</td></tr>
                <tr><th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th><td id="ld_phi_net">-</td></tr>
                <tr><th>&lang;R<sub>large</sub>&rang;</th><td id="ld_R_large">-</td></tr>
                <tr><th>&phi;<sub>large</sub> &times; 10<sup>4</sup></th><td id="ld_phi_large">-</td></tr>
                <tr><th>A<sub>porod</sub> &times; 10<sup>6</sup></th><td id="ld_A_porod">-</td></tr>
              </table>
            </div>
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 5: PARAMETER TABLES -->
    <div id="tables" class="tab-content">
      <div class="card">
        <div class="card-title">Short Distance Parameters Table (3.0 m, 50 frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>D<sub>m</sub></th>
                <th>&plusmn;1&sigma;</th>
                <th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th>
                <th>A<sub>porod</sub> &times; 10<sup>6</sup></th>
                <th>bkg [cm<sup>-1</sup>]</th>
              </tr>
            </thead>
            <tbody>
              {"".join(f"<tr><td>{int(r['frame_index'])}</td><td>{float(r['time_seconds']):.3f}</td><td>{float(r['chi2_red']):.3f}</td><td>{float(r['Dm']):.3f}</td><td>{float(r.get('std_err_Dm', 0.0)):.3f}</td><td>{float(r.get('phi_net', 0.0))*1e4:.2f}</td><td>{float(r.get('A_porod', 0.0))*1e6:.2f}</td><td>{float(r['bkg']):.4f}</td></tr>" for _, r in df_sd.iterrows())}
            </tbody>
          </table>
        </div>
      </div>

      <div class="card">
        <div class="card-title">Long Distance Parameters Table (31.0 m, 50 frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>D<sub>m</sub></th>
                <th>&plusmn;1&sigma;</th>
                <th>&xi; [nm]</th>
                <th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th>
                <th>R<sub>large</sub> [nm]</th>
                <th>&phi;<sub>large</sub> &times; 10<sup>4</sup></th>
                <th>A<sub>porod</sub> &times; 10<sup>6</sup></th>
              </tr>
            </thead>
            <tbody>
              {"".join(f"<tr><td>{int(r['frame_index'])}</td><td>{float(r['time_seconds']):.3f}</td><td>{float(r['chi2_red']):.3f}</td><td>{float(r['Dm']):.3f}</td><td>{float(r.get('std_err_Dm', 0.0)):.3f}</td><td>{float(r['xi']):.0f}</td><td>{float(r.get('phi_net', 0.0))*1e4:.2f}</td><td>{float(r.get('mean_radius_large', 3113.0)):.0f}</td><td>{float(r.get('total_volume_fraction_large', 0.0))*1e4:.2f}</td><td>{float(r.get('A_porod', 0.0))*1e6:.2f}</td></tr>" for _, r in df_ld.iterrows())}
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 6: PHYSICAL FORMALISM -->
    <div id="theory" class="tab-content">
      <div class="card">
        <div class="card-title">Physical Model Formulations: Spinodal Network & Secondary Large Spheres</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem;">
          At 40.0 mM CaCO<sub>3</sub> (pH 9.75), precipitation is captured by the coherent linear superposition of the
          mesoscale Teixeira mass-fractal network, secondary macroscopic dense liquid droplets/spheres, and asymptotic Porod surface scattering:
        </p>

        <div class="equation-box">
          $$I(q) = I_{{\\rm net}}(q; R_{{\\rm prim}}, \\sigma_{{\\rm prim}}, \\phi_{{\\rm net}}, D_m, \\xi) + I_{{\\rm large}}(q; R_{{\\rm large}}, \\sigma_{{\\rm large}}, \\phi_{{\\rm large}}) + \\frac{{A_{{\\rm porod}}}}{{q^4}} + \\text{{bkg}}$$
        </div>

        <h4 style="color: #fff; margin: 1.5rem 0 0.5rem 0;">1. Teixeira Mass Fractal Network</h4>
        <p style="color: var(--text-secondary);">
          Captures the self-similar spatial correlation of the spinodal network structure:
        </p>
        <div class="equation-box">
          $$I_{{\\rm net}}(q) = \\phi_{{\\rm net}} V_{{\\rm prim}} \\Delta\\rho^2 P_{{\\rm prim}}(q) S(q)$$
          $$S(q) = 1 + \\frac{{1}}{{(q R_{{\\rm prim}})^{{D_m}}}} \\frac{{D_m \\Gamma(D_m - 1)}}{{[1 + 1/(q\\xi)^2]^{{(D_m - 1)/2}}}} \\sin\\left[(D_m - 1) \\arctan(q\\xi)\\right]$$
        </div>

        <h4 style="color: #fff; margin: 1.5rem 0 0.5rem 0;">2. Secondary Macroscopic Droplets / Large Spheres</h4>
        <p style="color: var(--text-secondary);">
          Captures the low-q upturn and dense phase liquid droplets (R &approx; 3 &mu;m):
        </p>
        <div class="equation-box">
          $$I_{{\\rm large}}(q) = \\phi_{{\\rm large}} \\Delta\\rho^2 \\frac{{\\int_0^\\infty V(R)^2 f(R) P(q, R)\\, dR}}{{\\int_0^\\infty V(R) f(R)\\, dR}}$$
        </div>
      </div>
    </div>
  </main>

  <script>
    const shortData = {json.dumps(short_telemetry)};
    const longData = {json.dumps(long_telemetry)};

    let sdTimer = null;
    let ldTimer = null;

    function switchTab(tabId) {{
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));

      event.target.classList.add('active');
      document.getElementById(tabId).classList.add('active');
    }}

    function padFrame(f) {{
      return String(f).padStart(2, '0');
    }}

    function updateShortFrame(f) {{
      f = parseInt(f);
      const d = shortData[f];
      if (!d) return;

      document.getElementById('sd_frame_val').innerText = f;
      document.getElementById('sd_time_val').innerText = d.time + ' s';
      document.getElementById('sd_frame_badge').innerText = 'Frame ' + String(f).padStart(2, '0') + ' (' + d.time + ' s)';
      document.getElementById('sd_img').src = '../short_distance/fits/frame_' + padFrame(f) + '.png';

      document.getElementById('sd_chi2').innerText = d.chi2;
      document.getElementById('sd_Dm').innerText = d.Dm;
      document.getElementById('sd_phi_net').innerText = d.phi_net;
      document.getElementById('sd_A_porod').innerText = d.A_porod;
      document.getElementById('sd_R_prim').innerText = d.R_prim + ' nm';
      document.getElementById('sd_bkg').innerText = d.bkg;
    }}

    function stepShort(delta) {{
      const s = document.getElementById('sd_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 49) val = 49;
      s.value = val;
      updateShortFrame(val);
    }}

    function toggleShortPlay() {{
      const btn = document.getElementById('sd_play_btn');
      if (sdTimer) {{
        clearInterval(sdTimer);
        sdTimer = null;
        btn.innerText = 'Play ►';
      }} else {{
        btn.innerText = 'Pause ❚❚';
        sdTimer = setInterval(() => {{
          const s = document.getElementById('sd_slider');
          let val = (parseInt(s.value) + 1) % 50;
          s.value = val;
          updateShortFrame(val);
        }}, 150);
      }}
    }}

    function updateLongFrame(f) {{
      f = parseInt(f);
      const d = longData[f];
      if (!d) return;

      document.getElementById('ld_frame_val').innerText = f;
      document.getElementById('ld_time_val').innerText = d.time + ' s';
      document.getElementById('ld_frame_badge').innerText = 'Frame ' + String(f).padStart(2, '0') + ' (' + d.time + ' s)';
      document.getElementById('ld_img').src = '../long_distance/fits/frame_' + padFrame(f) + '.png';

      document.getElementById('ld_chi2').innerText = d.chi2;
      document.getElementById('ld_Dm').innerText = d.Dm;
      document.getElementById('ld_xi').innerText = d.xi + ' nm';
      document.getElementById('ld_phi_net').innerText = d.phi_net;
      document.getElementById('ld_R_large').innerText = d.R_large + ' nm';
      document.getElementById('ld_phi_large').innerText = d.phi_large;
      document.getElementById('ld_A_porod').innerText = d.A_porod;
    }}

    function stepLong(delta) {{
      const s = document.getElementById('ld_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 49) val = 49;
      s.value = val;
      updateLongFrame(val);
    }}

    function toggleLongPlay() {{
      const btn = document.getElementById('ld_play_btn');
      if (ldTimer) {{
        clearInterval(ldTimer);
        ldTimer = null;
        btn.innerText = 'Play ►';
      }} else {{
        btn.innerText = 'Pause ❚❚';
        ldTimer = setInterval(() => {{
          const s = document.getElementById('ld_slider');
          let val = (parseInt(s.value) + 1) % 50;
          s.value = val;
          updateLongFrame(val);
        }}, 150);
      }}
    }}

    updateShortFrame(0);
    updateLongFrame(0);
  </script>
</body>
</html>
"""

    dashboard_path = os.path.join(dashboard_dir, 'index.html')
    with open(dashboard_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Generated unified dashboard: {dashboard_path} ({len(html_content):,} bytes)")


if __name__ == '__main__':
    build_dashboard()
