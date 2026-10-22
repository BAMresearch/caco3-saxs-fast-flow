#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-004.
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
            'R_acc': f"{float(r['R_acc']):.2f} ± {float(r.get('std_err_R_acc', 0.0)):.2f}",
            'std_acc': f"{float(r['std_acc']):.2f} ± {float(r.get('std_err_std_acc', 0.0)):.2f}",
            'phi_acc': f"{float(r['phi_acc']) * 1e4:.2f} ± {float(r.get('std_err_phi_acc', 0.0)) * 1e4:.2f}",
            'phi_net': f"{float(r.get('phi_net', 0.0)) * 1e4:.2f} ± {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}",
            'A_porod': f"{float(r.get('A_porod', 0.0)) * 1e6:.2f} ± {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}",
            'R_prim': f"{float(r.get('R_prim', 3.5)):.2f}",
            'bkg': f"{float(r['bkg']):.4f} ± {float(r.get('std_err_bkg', 0.0)):.4f}"
        }

    # 2. Load Long Distance Parameters
    ld_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')
    df_ld = pd.read_csv(ld_csv)
    long_telemetry = {}
    for _, r in df_ld.iterrows():
        f = int(r['frame_index'])
        t_sec = float(r['time_seconds'])
        p_net = float(r['phi_net'])
        p_acc = float(r['phi_acc'])
        tot = p_net + p_acc
        c_net = (p_net / tot * 100.0) if tot > 0 else 0.0
        c_acc = (p_acc / tot * 100.0) if tot > 0 else 0.0
        dom = "Fractal Network" if c_net >= 50.0 else "ACC Nanospheres"
        long_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': f"{float(r['chi2_red']):.3f}",
            'Dm': f"{float(r['Dm']):.3f} ± {float(r.get('std_err_Dm', 0.0)):.3f}",
            'xi': f"{float(r['xi']):.0f} ± {float(r.get('std_err_xi', 0.0)):.0f}",
            'phi_net': f"{p_net * 1e4:.2f} ± {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}",
            'A_porod': f"{float(r['A_porod']) * 1e6:.2f} ± {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}",
            'R_acc_transferred': f"{float(r.get('R_acc_transferred', 0.0)):.2f}",
            'phi_acc': f"{p_acc * 1e4:.2f} ± {float(r.get('std_err_phi_acc', 0.0)) * 1e4:.2f}",
            'chi_net': f"{c_net:.1f}%",
            'chi_acc': f"{c_acc:.1f}%",
            'dominant_phase': dom,
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
  <title>ND-i-004 Multiscale Kinetic Analysis Dashboard</title>
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
      background: linear-gradient(180deg, #1e293b 0%, #0f172a 100%);
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
    .title-group h1 {{
      font-size: 1.5rem;
      font-weight: 700;
      letter-spacing: -0.025em;
      color: #fff;
    }}
    .title-group p {{
      color: var(--text-secondary);
      font-size: 0.875rem;
      margin-top: 0.25rem;
    }}
    .nav-tabs {{
      display: flex;
      gap: 0.5rem;
      border-bottom: 1px solid var(--border);
      background-color: var(--bg-secondary);
      padding: 0 2rem;
      overflow-x: auto;
    }}
    .tab-btn {{
      background: none;
      border: none;
      color: var(--text-secondary);
      padding: 0.75rem 1rem;
      font-size: 0.875rem;
      font-weight: 600;
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
      border-bottom-color: var(--accent);
    }}
    main {{
      max-width: 1400px;
      margin: 0 auto;
      padding: 2rem;
    }}
    .tab-content {{
      display: none;
    }}
    .tab-content.active {{
      display: block;
    }}
    .card {{
      background-color: var(--bg-card);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1.5rem;
      margin-bottom: 1.5rem;
    }}
    .card-title {{
      font-size: 1.125rem;
      font-weight: 600;
      margin-bottom: 1rem;
      color: #fff;
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
      display: flex;
      flex-direction: column;
    }}
    .viewer-panel img {{
      width: 100%;
      height: auto;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: #000;
    }}
    .slider-box {{
      margin: 1.5rem 0;
      display: flex;
      flex-direction: column;
      gap: 0.5rem;
    }}
    input[type=range] {{
      width: 100%;
      accent-color: var(--accent);
    }}
    .btn-group {{
      display: flex;
      gap: 0.5rem;
      margin-bottom: 1.5rem;
    }}
    .btn {{
      background-color: var(--border);
      color: #fff;
      border: none;
      padding: 0.5rem 1rem;
      border-radius: 4px;
      font-size: 0.85rem;
      cursor: pointer;
      font-weight: 600;
      transition: background 0.2s;
    }}
    .btn:hover {{
      background-color: var(--accent);
    }}
    .telemetry-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.85rem;
    }}
    .telemetry-table tr {{
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}
    .telemetry-table th {{
      text-align: left;
      padding: 0.5rem 0;
      color: var(--text-secondary);
      font-weight: 500;
    }}
    .telemetry-table td {{
      text-align: right;
      padding: 0.5rem 0;
      font-weight: 600;
      color: #fff;
    }}
    .badge {{
      background-color: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      padding: 0.25rem 0.5rem;
      border-radius: 4px;
      font-size: 0.75rem;
      font-weight: 600;
    }}
    .data-table-container {{
      max-height: 500px;
      overflow-y: auto;
      border: 1px solid var(--border);
      border-radius: 6px;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.825rem;
    }}
    .data-table th, .data-table td {{
      padding: 0.6rem 0.8rem;
      text-align: right;
      border-bottom: 1px solid var(--border);
    }}
    .data-table th {{
      background-color: #0b1120;
      position: sticky;
      top: 0;
      color: var(--accent);
      font-weight: 600;
    }}
    .data-table tr:hover {{
      background-color: rgba(255, 255, 255, 0.02);
    }}
    .data-table td:first-child, .data-table th:first-child {{
      text-align: left;
    }}
    .equation-box {{
      background-color: #0b1120;
      border-left: 3px solid var(--accent);
      padding: 1rem;
      margin: 1rem 0;
      border-radius: 0 4px 4px 0;
      overflow-x: auto;
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-content">
      <div class="title-group">
        <h1>ND-i-004 Kinetic SAXS Suite</h1>
        <p>15.0 mM CaCO<sub>3</sub>, pH 11.3 &mdash; Hybrid Spinodal Decomposition & ACC Nanosphere Precipitation</p>
      </div>
      <div style="display: flex; gap: 0.5rem; align-items: center;">
        <a href="../../../index.html" class="badge" style="text-decoration: none; cursor: pointer;">&larr; Master Dashboard</a>
        <span class="badge" style="background: rgba(16, 185, 129, 0.2); color: #10b981;">Validated Hybrid Model</span>
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
          At 15.0 mM initial supersaturation and pH 11.3, calcium carbonate precipitation occurs directly at the physical boundary
          between <strong>mesoscale spinodal gel network formation</strong> and <strong>discrete amorphous calcium carbonate (ACC) nanosphere precipitation</strong>.
          The complete multi-scale composite scattering curve spans four decades in momentum transfer (<em>q</em> &isin; [0.002, 2.5] nm<sup>-1</sup>), capturing both length scale regimes.
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
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Building block radius of the spinodal network</td>
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
              <td>ACC Nanosphere Radius</td><td>&lang;R<sub>ACC</sub>&rang;</td>
              <td>{bp['R_acc']:.1f} nm</td><td>&plusmn; {bp.get('std_err_R_acc', 0.0):.1f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Mean radius of amorphous calcium carbonate nanospheres</td>
            </tr>
            <tr>
              <td>ACC Polydispersity</td><td>&sigma;<sub>ACC</sub></td>
              <td>{bp['std_acc']:.1f} nm</td><td>&plusmn; {bp.get('std_err_std_acc', 0.0):.1f} nm</td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Standard deviation of Lognormal size distribution</td>
            </tr>
            <tr>
              <td>ACC Volume Fraction</td><td>&phi;<sub>ACC</sub></td>
              <td>{bp['phi_acc'] * 1e5:.2f} &times; 10<sup>-5</sup></td><td>&plusmn; {bp.get('std_err_phi_acc', 0.0) * 1e5:.2f} &times; 10<sup>-5</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Volume fraction of discrete ACC nanospheres</td>
            </tr>
            <tr>
              <td>Porod Surface Amplitude</td><td>A<sub>porod</sub></td>
              <td>{bp['A_porod'] * 1e6:.2f} &times; 10<sup>-6</sup></td><td>&plusmn; {bp.get('std_err_A_porod', 0.0) * 1e6:.2f} &times; 10<sup>-6</sup></td>
              <td style="text-align: left; padding-left: 1.5rem; color: var(--text-secondary);">Macroscopic coarsening domain surface scattering</td>
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
      <div class="card" style="border-left: 4px solid var(--accent); margin-bottom: 1.5rem;">
        <div class="card-title" style="color: var(--accent);">Relative Component Proportions and Dynamic Phase Inversion (Crossover Kinetics)</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem; font-size: 0.95rem; line-height: 1.6;">
          At 15.0 mM initial supersaturation, time-resolved measurements quantify the dynamic phase inversion between 
          <strong style="color: #38bdf8;">discrete ACC spheres</strong> (&chi;<sub>ACC</sub>) and the 
          <strong style="color: #34d399;">mesoscale spinodal fractal network</strong> (&chi;<sub>net</sub>):
        </p>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1rem; margin-bottom: 1.25rem;">
          <div style="background: rgba(2, 132, 199, 0.1); border: 1px solid rgba(2, 132, 199, 0.3); border-radius: 8px; padding: 1rem;">
            <div style="color: #38bdf8; font-weight: 700; font-size: 0.85rem; text-transform: uppercase;">1. Initial Stage (t = 0.060 s)</div>
            <div style="font-size: 1.25rem; font-weight: 700; color: #fff; margin: 0.25rem 0;">95.8% ACC Spheres</div>
            <p style="color: var(--text-secondary); font-size: 0.85rem;">Initial nucleation of discrete ACC spheres (&lang;R<sub>ACC</sub>&rang; &approx; 15.9 nm) comprises the primary precipitate fraction.</p>
          </div>
          <div style="background: rgba(245, 158, 11, 0.1); border: 1px solid rgba(245, 158, 11, 0.3); border-radius: 8px; padding: 1rem;">
            <div style="color: #fbbf24; font-weight: 700; font-size: 0.85rem; text-transform: uppercase;">2. Dynamic Crossover (t &approx; 0.15 to 0.25 s)</div>
            <div style="font-size: 1.25rem; font-weight: 700; color: #fff; margin: 0.25rem 0;">50 / 50 Parity Inversion</div>
            <p style="color: var(--text-secondary); font-size: 0.85rem;">Mesoscale spinodal network formation develops across the solution volume; sphere radius growth arrests at &sim;44 nm.</p>
          </div>
          <div style="background: rgba(5, 150, 105, 0.1); border: 1px solid rgba(5, 150, 105, 0.3); border-radius: 8px; padding: 1rem;">
            <div style="color: #34d399; font-weight: 700; font-size: 0.85rem; text-transform: uppercase;">3. Terminal Regime (t &gt; 1.0 s to 10.9 s)</div>
            <div style="font-size: 1.25rem; font-weight: 700; color: #fff; margin: 0.25rem 0;">76.3% Fractal Network</div>
            <p style="color: var(--text-secondary); font-size: 0.85rem;">Percolated spinodal gel network (D<sub>m</sub> &approx; 2.33, &xi; &approx; 197 nm) comprises the predominant precipitate fraction; discrete ACC spheres remain as an arrested secondary population.</p>
          </div>
        </div>
        <img src="../unified_trajectories/relative_phase_evolution.png" alt="Relative Phase Evolution" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>

      <div class="card">
        <div class="card-title">Multiscale Kinetic Trajectories Synthesis ($t = 0.060\text{{ s}} \to 10.88\text{{ s}}$)</div>
        <p style="color: var(--text-secondary); margin-bottom: 1rem; font-size: 0.9rem;">
          Cross-scale synthesis combining 3.0 m (WAXS/MAXS) and 31.0 m (SAXS) measurements with shaded &plusmn;1&sigma; uncertainty envelopes and standardised component palette.
        </p>
        <img src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Trajectories" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
      </div>
    </div>

    <!-- TAB 3: SHORT DISTANCE INSPECTOR -->
    <div id="short_inspector" class="tab-content">
      <div class="inspector-grid">
        <div class="viewer-panel card">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
            <div class="card-title" style="margin-bottom: 0;">Short-Distance Scattering Pattern Fit (3.0 m Geometry)</div>
            <span class="badge" id="sd_frame_badge">Frame 00 (0.060 s)</span>
          </div>
          <img id="sd_img" src="../short_distance/fits/frame_00.png" alt="Short Distance Fit">
        </div>

        <div class="card">
          <div class="card-title">Frame Controller & Telemetry</div>
          
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
            <tr><th>ACC Radius &lang;R<sub>ACC</sub>&rang;</th><td id="sd_R_acc">-</td></tr>
            <tr><th>ACC Polydispersity &sigma;<sub>ACC</sub></th><td id="sd_std_acc">-</td></tr>
            <tr><th>ACC Fraction &phi;<sub>ACC</sub> &times; 10<sup>4</sup></th><td id="sd_phi_acc">-</td></tr>
            <tr><th>Network Fraction &phi;<sub>net</sub> &times; 10<sup>4</sup></th><td id="sd_phi_net">-</td></tr>
            <tr><th>Building Block R<sub>prim</sub></th><td id="sd_R_prim">-</td></tr>
            <tr><th>Porod Amplitude A<sub>porod</sub> &times; 10<sup>6</sup></th><td id="sd_A_porod">-</td></tr>
            <tr><th>Flat Background</th><td id="sd_bkg">-</td></tr>
          </table>

          <div style="font-size: 0.75rem; color: var(--text-secondary); line-height: 1.4; margin-top: 1rem;">
            Physically unified model: Polydisperse ACC nanospheres (Lognormal) + Spinodal network primary building blocks (3.5 nm) + Porod surface tail.
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 4: LONG DISTANCE INSPECTOR -->
    <div id="long_inspector" class="tab-content">
      <div class="inspector-grid">
        <div class="viewer-panel card">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
            <div class="card-title" style="margin-bottom: 0;">Long-Distance Scattering Pattern Fit (31.0 m Geometry)</div>
            <span class="badge" id="ld_frame_badge">Frame 00 (0.060 s)</span>
          </div>
          <img id="ld_img" src="../long_distance/fits/frame_00.png" alt="Long Distance Fit">
        </div>

        <div class="card">
          <div class="card-title">Frame Controller & Telemetry</div>
          
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
            <tr><th>Fractal Dimension D<sub>m</sub></th><td id="ld_Dm">-</td></tr>
            <tr><th>Network Cutoff &xi;</th><td id="ld_xi">-</td></tr>
            <tr><th>Network Fraction &phi;<sub>net</sub> &times; 10<sup>4</sup></th><td id="ld_phi_net">-</td></tr>
            <tr><th>Porod Amplitude A<sub>porod</sub> &times; 10<sup>6</sup></th><td id="ld_A_porod">-</td></tr>
            <tr><th>Transferred &lang;R<sub>ACC</sub>&rang;</th><td id="ld_R_acc">-</td></tr>
            <tr><th>ACC Fraction &phi;<sub>ACC</sub> &times; 10<sup>4</sup></th><td id="ld_phi_acc">-</td></tr>
            <tr><th>Relative Network Share &chi;<sub>net</sub></th><td id="ld_chi_net" style="color: #34d399; font-weight: 600;">-</td></tr>
            <tr><th>Relative ACC Share &chi;<sub>ACC</sub></th><td id="ld_chi_acc" style="color: #38bdf8; font-weight: 600;">-</td></tr>
            <tr><th>Dominant Phase</th><td id="ld_dominant" style="font-weight: 700;">-</td></tr>
          </table>

          <div style="font-size: 0.75rem; color: var(--text-secondary); line-height: 1.4; margin-top: 1rem;">
            Transferred ACC nanosphere population: seamlessly coupled from Short Distance.
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
                <th>&lang;R<sub>ACC</sub>&rang; [nm]</th>
                <th>&sigma;<sub>ACC</sub> [nm]</th>
                <th>&phi;<sub>ACC</sub> &times; 10<sup>4</sup></th>
                <th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th>
                <th>A<sub>porod</sub> &times; 10<sup>6</sup></th>
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
                <td>{float(r['R_acc']):.2f} &plusmn; {float(r.get('std_err_R_acc', 0.0)):.2f}</td>
                <td>{float(r['std_acc']):.2f} &plusmn; {float(r.get('std_err_std_acc', 0.0)):.2f}</td>
                <td>{float(r['phi_acc']) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi_acc', 0.0)) * 1e4:.2f}</td>
                <td>{float(r.get('phi_net', 0.0)) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}</td>
                <td>{float(r.get('A_porod', 0.0)) * 1e6:.2f} &plusmn; {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}</td>
                <td>{float(r['bkg']):.4f} &plusmn; {float(r.get('std_err_bkg', 0.0)):.4f}</td>
              </tr>
"""

    html_content += """            </tbody>
          </table>
        </div>
      </div>

      <div class="card" style="margin-top: 1.5rem;">
        <div class="card-title">Long-Distance Fitted Parameters (31.0 m, 50 Frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>D<sub>m</sub></th>
                <th>&xi; [nm]</th>
                <th>&phi;<sub>net</sub> &times; 10<sup>4</sup></th>
                <th>A<sub>porod</sub> &times; 10<sup>6</sup></th>
                <th>&lang;R<sub>ACC</sub>&rang; [nm]</th>
                <th>&phi;<sub>ACC</sub> &times; 10<sup>4</sup></th>
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
                <td>{float(r['Dm']):.3f} &plusmn; {float(r.get('std_err_Dm', 0.0)):.3f}</td>
                <td>{float(r['xi']):.0f} &plusmn; {float(r.get('std_err_xi', 0.0)):.0f}</td>
                <td>{float(r['phi_net']) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi_net', 0.0)) * 1e4:.2f}</td>
                <td>{float(r['A_porod']) * 1e6:.2f} &plusmn; {float(r.get('std_err_A_porod', 0.0)) * 1e6:.2f}</td>
                <td>{float(r.get('R_acc_transferred', 0.0)):.2f}</td>
                <td>{float(r['phi_acc']) * 1e4:.2f} &plusmn; {float(r.get('std_err_phi_acc', 0.0)) * 1e4:.2f}</td>
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
        <div class="card-title">Hybrid Multiscale Physical Scattering Formalism</div>
        <p style="color: var(--text-secondary); margin-bottom: 12px;">
          At 15.0 mM CaCO<sub>3</sub> (pH 11.3), phase separation occurs directly at the crossover phase boundary between
          spinodal decomposition / gel network formation and discrete amorphous calcium carbonate (ACC) precipitation.
          The total scattering intensity \(I_{\mathrm{tot}}(q)\) is modelled consistently across both camera geometries
          and the composite dataset \(q \in [0.002, 2.5]\text{ nm}^{-1}\) as:
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{tot}}(q) = I_{\mathrm{fractal}}(q; R_{\mathrm{prim}}, \sigma_{\mathrm{prim}}, D_m, \xi, \phi_{\mathrm{net}}) + I_{\mathrm{ACC}}(q; \langle R_{\mathrm{ACC}} \rangle, \sigma_{\mathrm{ACC}}, \phi_{\mathrm{ACC}}) + I_{\mathrm{porod}}(q; A_{\mathrm{porod}}) + I_{\mathrm{bkg}}
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">1. Spinodal Mass Fractal Network & Primary Building Blocks</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Primary building blocks (\(R_{\mathrm{prim}} = 3.5\text{ nm}\), \(\sigma_{\mathrm{prim}} = 1.5\text{ nm}\)) interconnected via a 3D percolated mass fractal network described by Teixeira's structure factor \(S(q)\). At short distance (\(q > 0.1\text{ nm}^{-1}\)), \(S(q) \to 1\), capturing the building block form factor \(P(q)\):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{fractal}}(q) = (\Delta\rho)^2 \int_0^\infty n(R) V^2(R) P(q, R) S(q, R, D_m, \xi) \, dR
          \]
          \[
          S(q, R, D_m, \xi) = 1 + \frac{D_m \Gamma(D_m - 1)}{(q R)^{D_m}} \frac{\sin[(D_m - 1) \arctan(q \xi)]}{\left(1 + (q \xi)^{-2}\right)^{(D_m - 1)/2}}
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">2. Discrete Polydisperse ACC Nanospheres</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Amorphous calcium carbonate nanospheres growing via nucleation & growth with narrow polydispersity (\(\sim 17\%\)), producing form-factor Bessel oscillations:
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{ACC}}(q) = (\Delta\rho)^2 \int_0^\infty n_{\mathrm{ACC}}(R) V^2(R) P(q, R) \, dR
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">3. Asymptotic Porod Surface Scattering</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Macroscopic domain surface scattering from large phase-separated coarsened regions (\(q R \gg 1\)):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{porod}}(q) = \frac{A_{\mathrm{porod}}}{q^4}
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
      return String(f).padStart(2, '0');
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
      document.getElementById('sd_R_acc').innerText = d.R_acc + ' nm';
      document.getElementById('sd_std_acc').innerText = d.std_acc + ' nm';
      document.getElementById('sd_phi_acc').innerText = d.phi_acc;
      document.getElementById('sd_phi_net').innerText = d.phi_net;
      document.getElementById('sd_R_prim').innerText = d.R_prim + ' nm';
      document.getElementById('sd_A_porod').innerText = d.A_porod;
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
      document.getElementById('ld_Dm').innerText = d.Dm;
      document.getElementById('ld_xi').innerText = d.xi + ' nm';
      document.getElementById('ld_phi_net').innerText = d.phi_net;
      document.getElementById('ld_A_porod').innerText = d.A_porod;
      document.getElementById('ld_R_acc').innerText = d.R_acc_transferred + ' nm';
      document.getElementById('ld_phi_acc').innerText = d.phi_acc;
      document.getElementById('ld_chi_net').innerText = d.chi_net;
      document.getElementById('ld_chi_acc').innerText = d.chi_acc;
      const domEl = document.getElementById('ld_dominant');
      domEl.innerText = d.dominant_phase;
      domEl.style.color = d.dominant_phase === 'Fractal Network' ? '#34d399' : '#38bdf8';
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

    dashboard_path = os.path.join(dashboard_dir, 'index.html')
    with open(dashboard_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Generated unified dashboard: {dashboard_path} ({len(html_content):,} bytes)")


if __name__ == '__main__':
    build_dashboard()
