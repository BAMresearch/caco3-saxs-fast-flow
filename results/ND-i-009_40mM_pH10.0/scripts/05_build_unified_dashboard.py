#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de
"""
Script 05: Build Unified Multi-Scale Interactive Dashboard for ND-i-009.
-----------------------------------------------------------------------
Generates a standalone, responsive HTML dashboard for reviewing:
1. Overview & Composite Curve (q in [0.002, 2.5] nm^-1)
2. Unified Cross-Scale Kinetic Trajectories (0.045 s to 5.295 s)
3. Short-Distance Time-Resolved Inspector (3.0 m, Frames 0 to 350, with live slider & +/- 1 sigma telemetry)
4. Long-Distance Time-Resolved Inspector (31.0 m, Frames 0 to 349, with live slider & +/- 1 sigma telemetry)
5. Comprehensive Parameter Tables (searchable & scrollable)
6. Multiscale Physical Scattering Formalism & Equations rendered via MathJax.

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

    # 1. Load Short Distance Parameters (Flowing Viscoelastic VPS)
    p_visco_csv = os.path.join(PROJECT_DIR, 'output', 'short_distance', 'parameters.csv')
    df_visco = pd.read_csv(p_visco_csv)
    short_telemetry = {}
    for _, r in df_visco.iterrows():
        f = int(r['frame'])
        t_sec = float(r.get('time_seconds', 0.045 + 0.015 * f))

        r_mean = float(r['Rmean'])
        err_r = float(r.get('std_err_Rmean', 0.0))
        r_str = f"{r_mean:.2f} ± {err_r:.2f}" if err_r > 0 else f"{r_mean:.2f}"

        d_val = float(r['D'])
        err_d = float(r.get('std_err_D', 0.0))
        d_str = f"{d_val:.4f} ± {err_d:.4f}" if err_d > 0 else f"{d_val:.4f}"

        xi_fluc = float(r.get('xi_fluc', 2.80))
        err_xi_fluc = float(r.get('std_err_xi_fluc', 0.0))
        xi_fluc_str = f"{xi_fluc:.2f} ± {err_xi_fluc:.2f}" if err_xi_fluc > 0 else f"{xi_fluc:.2f}"

        a_fluc = float(r.get('A_fluc', 0.0))
        err_a_fluc = float(r.get('std_err_A_fluc', 0.0))
        a_fluc_str = f"{a_fluc:.3f} ± {err_a_fluc:.3f}" if err_a_fluc > 0 else f"{a_fluc:.3f}"

        ratio = float(r.get('ratio_A_xi_83', 0.0))
        err_ratio = float(r.get('std_err_ratio', 0.0))
        ratio_str = f"{ratio:.5f} ± {err_ratio:.5f}" if err_ratio > 0 else f"{ratio:.5f}"

        phi = float(r['phi'])
        err_phi = float(r.get('std_err_phi', 0.0))
        phi_str = f"{phi:.2e} ± {err_phi:.2e}" if err_phi > 0 else f"{phi:.2e}"

        short_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': round(float(r['redchi']), 3),
            'hq_chi2': round(float(r.get('hq_redchi', 0.0)), 3),
            'Rmean': r_str,
            'D': d_str,
            'phi': phi_str,
            'xi': round(float(r.get('xi', 1e7)), 1) if float(r.get('xi', 1e7)) < 1e5 else '∞',
            'xi_fluc': xi_fluc_str,
            'A_fluc': a_fluc_str,
            'ratio': ratio_str
        }

    # 2. Load Long Distance Parameters (Option C Forward)
    ld_fwd_csv = os.path.join(PROJECT_DIR, 'output', 'long_distance', 'parameters.csv')
    df_fwd = pd.read_csv(ld_fwd_csv)
    long_telemetry = {}
    for _, r_f in df_fwd.iterrows():
        f = int(r_f['frame'])
        t_sec = float(r_f.get('time_seconds', 0.045 + 0.015 * f))
        has_l = bool(r_f.get('has_large_spheres', False)) and float(r_f.get('total_volume_fraction_large', 0.0)) > 0

        d_val = float(r_f['D'])
        err_d = float(r_f.get('std_err_D', 0.0))
        d_str = f"{d_val:.3f} ± {err_d:.3f}" if err_d > 0 else f"{d_val:.3f}"

        phi = float(r_f['phi'])
        err_phi = float(r_f.get('std_err_phi', 0.0))
        phi_str = f"{phi:.2e} ± {err_phi:.2e}" if err_phi > 0 else f"{phi:.2e}"

        r_large = float(r_f.get('mean_radius_large', 0.0))
        err_r_large = float(r_f.get('std_err_mean_radius_large', 0.0))
        r_large_str = f"{r_large:.0f} ± {err_r_large:.0f}" if err_r_large > 0 and has_l else f"{r_large:.0f}"

        phi_large = float(r_f.get('total_volume_fraction_large', 0.0))
        err_phi_large = float(r_f.get('std_err_total_volume_fraction_large', 0.0))
        phi_large_str = f"{phi_large:.2e} ± {err_phi_large:.2e}" if err_phi_large > 0 and has_l else f"{phi_large:.2e}"

        long_telemetry[f] = {
            'time': f"{t_sec:.3f}",
            'chi2': round(float(r_f['redchi']), 3),
            'D': d_str,
            'phi': phi_str,
            'xi': '∞',
            'has_large': has_l,
            'R_large': r_large_str if has_l else "0",
            'std_large': round(float(r_f.get('std_radius_large', 0.0)), 0) if has_l else 0,
            'phi_large': phi_large_str if has_l else '0.00e+00'
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
  <title>ND-i-009 Multiscale Kinetic Analysis Dashboard</title>
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
      gap: 6px;
    }}
    .btn {{
      flex: 1;
      padding: 8px 10px;
      border-radius: 6px;
      font-weight: 600;
      font-size: 0.82rem;
      cursor: pointer;
      border: 1px solid var(--border);
      background-color: #273549;
      color: var(--text-primary);
      transition: all 0.15s ease;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 4px;
    }}
    .btn:hover {{
      background-color: var(--accent);
      color: #0f172a;
      border-color: var(--accent);
    }}
    .btn-play {{
      background-color: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      border-color: rgba(56, 189, 248, 0.4);
    }}
    .btn-pause {{
      background-color: rgba(239, 68, 68, 0.2);
      color: var(--danger);
      border-color: rgba(239, 68, 68, 0.4);
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
      width: 52%;
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
    .hero-img {{
      width: 100%;
      border-radius: 8px;
      border: 1px solid var(--border);
      margin-top: 12px;
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-title">
      <h1>ND-i-009 Multiscale Kinetic Analysis Dashboard</h1>
      <p>CaCO<sub>3</sub> Viscoelastic Phase Separation & Droplet Growth &bull; High-Speed Time-Resolved SAXS (0.045 s to 5.295 s)</p>
    </div>
    <div style="display: flex; gap: 10px; align-items: center; flex-wrap: wrap;">
      <a href="../../../index.html" class="badge" style="text-decoration: none; cursor: pointer;">&larr; Master Dashboard</a>
      <span class="badge badge-success">351 Time Frames</span>
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
    <!-- TAB 1: OVERVIEW & COMPOSITE CURVE -->
    <div id="overview" class="tab-content active">
      <div class="grid-2">
        <div class="card">
          <div class="card-title">Full Multiscale Composite SAXS Scattering Curve</div>
          <p style="color: var(--text-secondary); margin-bottom: 16px; font-size: 0.9rem;">
            Stitched wide-angle composite scattering curve combining the final reaction stages across both detector distances (q: 0.002 to 2.5 nm<sup>-1</sup>).
            Resolves primary amorphous building blocks, the mass fractal network, secondary macroscopic droplets, and dynamic concentration fluctuations.
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
            <tr><th>Network Correlation Cutoff &xi;</th><td>{bp['xi']:.1f} &plusmn; {bp.get('std_err_xi', 0.0):.1f} nm</td></tr>
            <tr><th>Secondary Droplet Radius &lang;R<sub>large</sub>&rang;</th><td>{bp['mean_radius_large']:.0f} &plusmn; {bp.get('std_err_mean_radius_large', 0.0):.0f} nm</td></tr>
            <tr><th>Secondary Droplet Dispersion &sigma;<sub>large</sub></th><td>{bp['std_radius_large']:.0f} &plusmn; {bp.get('std_err_std_radius_large', 0.0):.0f} nm</td></tr>
            <tr><th>Droplet Volume Fraction &phi;<sub>large</sub></th><td>{bp['total_volume_fraction_large'] * 1e3:.2f} &plusmn; {bp.get('std_err_total_volume_fraction_large', 0.0) * 1e3:.2f} &times; 10<sup>-3</sup></td></tr>
            <tr><th>OZ Fluctuation Amplitude A<sub>OZ</sub></th><td>{bp['A_OZ']:.1f} &plusmn; {bp.get('std_err_A_OZ', 0.0):.1f}</td></tr>
            <tr><th>OZ Fluctuation Length &xi;<sub>OZ</sub></th><td>{bp['xi_OZ']:.1f} &plusmn; {bp.get('std_err_xi_OZ', 0.0):.1f} nm</td></tr>
            <tr><th>Fractal Volume Fraction &phi;<sub>prim</sub></th><td>{bp['phi'] * 1e3:.2f} &plusmn; {bp.get('std_err_phi', 0.0) * 1e3:.2f} &times; 10<sup>-3</sup></td></tr>
            <tr><th>Overall Baseline Fit Quality &chi;<sup>2</sup><sub>red</sub></th><td><span class="badge badge-success">{bp['redchi']:.3f}</span></td></tr>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 2: CROSS-SCALE KINETICS -->
    <div id="cross_scale" class="tab-content">
      <div class="card">
        <div class="card-title">Unified Multi-Scale Kinetic Synthesis (0.045 s to 5.295 s)</div>
        <p style="color: var(--text-secondary); margin-bottom: 16px; font-size: 0.9rem;">
          Cross-scale synchronization demonstrating consistent structural length scale hierarchy across 4 decades,
          rigorous agreement between short-distance and long-distance fractal dimensions (D<sub>f</sub> &plusmn; 1&sigma;), and
          smooth continuous growth of macroscopic droplets.
        </p>
        <div class="viewer-box">
          <img class="viewer-img" src="../unified_trajectories/unified_cross_scale_trajectories.png" alt="Unified Cross-Scale Trajectories">
        </div>
      </div>

      <div class="grid-2">
        <div class="card">
          <div class="card-title">Short-Distance Trajectories (3.0 m Camera)</div>
          <img class="hero-img" src="../short_distance/viscoelastic/trajectories.png" alt="Short Distance Trajectories">
        </div>
        <div class="card">
          <div class="card-title">Long-Distance Trajectories (31.0 m Camera)</div>
          <img class="hero-img" src="../long_distance/option_c_forward/trajectories.png" alt="Long Distance Trajectories">
        </div>
      </div>
    </div>

    <!-- TAB 3: SHORT DISTANCE INSPECTOR -->
    <div id="short_inspector" class="tab-content">
      <div class="inspector-layout">
        <div class="card" style="margin-bottom: 0;">
          <div class="card-title">
            <span>Short-Distance Frame Fit (3.0 m Camera)</span>
            <span id="sd_frame_badge" class="badge">Frame 000 (0.045 s)</span>
          </div>
          <div class="viewer-box">
            <img id="sd_img" class="viewer-img" src="../short_distance/fits/fit_frame_000.png" alt="Short Distance Fit">
          </div>
        </div>

        <div class="controls-panel">
          <div class="slider-row">
            <div class="slider-labels">
              <span>Time Frame: <b id="sd_frame_val" style="color: var(--accent);">0</b> / 350</span>
              <span>Reaction Time: <b id="sd_time_val" style="color: #ffffff;">0.045 s</b></span>
            </div>
            <input type="range" id="sd_slider" min="0" max="350" value="0" oninput="updateShortFrame(this.value)">
          </div>

          <div class="btn-group">
            <button class="btn" onclick="stepShort(-10)">&minus;10</button>
            <button class="btn" onclick="stepShort(-1)">&minus;1</button>
            <button class="btn btn-play" id="sd_play_btn" onclick="toggleShortPlay()">Play</button>
            <button class="btn" onclick="stepShort(1)">+1</button>
            <button class="btn" onclick="stepShort(10)">+10</button>
          </div>

          <table class="telemetry-table">
            <tr><th>Goodness of Fit &chi;<sup>2</sup><sub>red</sub></th><td id="sd_chi2">-</td></tr>
            <tr><th>High-q &chi;<sup>2</sup><sub>red</sub></th><td id="sd_hq_chi2">-</td></tr>
            <tr><th>Primary Core &lang;R<sub>prim</sub>&rang;</th><td id="sd_Rmean">-</td></tr>
            <tr><th>Fractal Dimension D<sub>f</sub></th><td id="sd_D">-</td></tr>
            <tr><th>Network Cutoff &xi;</th><td id="sd_xi">-</td></tr>
            <tr><th>Correlation Length &xi;<sub>fluc</sub></th><td id="sd_xi_fluc">-</td></tr>
            <tr><th>Fluctuation Ampl. A<sub>fluc</sub></th><td id="sd_A_fluc">-</td></tr>
            <tr><th>Scaling Ratio A/&xi;<sup>8/3</sup></th><td id="sd_ratio">-</td></tr>
            <tr><th>Volume Fraction &phi;</th><td id="sd_phi">-</td></tr>
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
            <span id="ld_frame_badge" class="badge">Frame 000 (0.045 s)</span>
          </div>
          <div class="viewer-box">
            <img id="ld_img" class="viewer-img" src="../long_distance/fits/fit_frame_000.png" alt="Long Distance Fit">
          </div>
        </div>

        <div class="controls-panel">
          <div class="slider-row">
            <div class="slider-labels">
              <span>Time Frame: <b id="ld_frame_val" style="color: var(--accent);">0</b> / 349</span>
              <span>Reaction Time: <b id="ld_time_val" style="color: #ffffff;">0.045 s</b></span>
            </div>
            <input type="range" id="ld_slider" min="0" max="349" value="0" oninput="updateLongFrame(this.value)">
          </div>

          <div class="btn-group">
            <button class="btn" onclick="stepLong(-10)">&minus;10</button>
            <button class="btn" onclick="stepLong(-1)">&minus;1</button>
            <button class="btn btn-play" id="ld_play_btn" onclick="toggleLongPlay()">Play</button>
            <button class="btn" onclick="stepLong(1)">+1</button>
            <button class="btn" onclick="stepLong(10)">+10</button>
          </div>

          <table class="telemetry-table">
            <tr><th>Goodness of Fit &chi;<sup>2</sup><sub>red</sub></th><td id="ld_chi2">-</td></tr>
            <tr><th>Fractal Dimension D<sub>f</sub></th><td id="ld_D">-</td></tr>
            <tr><th>Network Cutoff &xi;</th><td id="ld_xi">-</td></tr>
            <tr><th>Has Large Droplets</th><td id="ld_haslarge">-</td></tr>
            <tr><th>Droplet Radius &lang;R<sub>large</sub>&rang;</th><td id="ld_rlarge">-</td></tr>
            <tr><th>Droplet Dispersion &sigma;<sub>large</sub></th><td id="ld_stdlarge">-</td></tr>
            <tr><th>Droplet Fraction &phi;<sub>large</sub></th><td id="ld_philarge">-</td></tr>
            <tr><th>Fractal Fraction &phi;</th><td id="ld_phi">-</td></tr>
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
        <div class="card-title">Short-Distance Fitted Parameters (3.0 m, 351 Frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>&lang;R<sub>prim</sub>&rang; [nm]</th>
                <th>D<sub>f</sub></th>
                <th>&xi; [nm]</th>
                <th>&xi;<sub>fluc</sub> [nm]</th>
                <th>A<sub>fluc</sub></th>
                <th>A/&xi;<sup>8/3</sup></th>
                <th>&phi;</th>
              </tr>
            </thead>
            <tbody>
"""

    for _, r in df_visco.iterrows():
        f = int(r['frame'])
        t_sec = float(r.get('time_seconds', 0.045 + 0.015 * f))
        r_mean = float(r['Rmean'])
        err_r = float(r.get('std_err_Rmean', 0.0))
        d_val = float(r['D'])
        err_d = float(r.get('std_err_D', 0.0))
        xi_val = float(r.get('xi', 1e7))
        xi_str = f"{xi_val:.1f}" if xi_val < 1e5 else "∞"
        xi_fluc = float(r.get('xi_fluc', 2.80))
        err_xi_fluc = float(r.get('std_err_xi_fluc', 0.0))
        a_fluc = float(r.get('A_fluc', 0.0))
        err_a_fluc = float(r.get('std_err_A_fluc', 0.0))
        ratio = float(r.get('ratio_A_xi_83', 0.0))
        phi = float(r['phi'])

        html_content += f"""              <tr>
                <td>Frame {f:03d}</td>
                <td>{t_sec:.3f}</td>
                <td>{float(r['redchi']):.3f}</td>
                <td>{r_mean:.2f} &plusmn; {err_r:.2f}</td>
                <td>{d_val:.4f} &plusmn; {err_d:.4f}</td>
                <td>{xi_str}</td>
                <td>{xi_fluc:.2f} &plusmn; {err_xi_fluc:.2f}</td>
                <td>{a_fluc:.3f} &plusmn; {err_a_fluc:.3f}</td>
                <td>{ratio:.5f}</td>
                <td>{phi:.2e}</td>
              </tr>
"""

    html_content += """            </tbody>
          </table>
        </div>
      </div>

      <div class="card">
        <div class="card-title">Long-Distance Fitted Parameters (31.0 m, 350 Frames)</div>
        <div class="data-table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Frame</th>
                <th>Time [s]</th>
                <th>&chi;<sup>2</sup><sub>red</sub></th>
                <th>D<sub>f</sub></th>
                <th>&xi; [nm]</th>
                <th>Has Large</th>
                <th>&lang;R<sub>large</sub>&rang; [nm]</th>
                <th>&sigma;<sub>large</sub> [nm]</th>
                <th>&phi;<sub>large</sub></th>
                <th>&phi;</th>
              </tr>
            </thead>
            <tbody>
"""

    for _, r_f in df_fwd.iterrows():
        f = int(r_f['frame'])
        t_sec = float(r_f.get('time_seconds', 0.045 + 0.015 * f))
        has_l = bool(r_f.get('has_large_spheres', False)) and float(r_f.get('total_volume_fraction_large', 0.0)) > 0
        d_val = float(r_f['D'])
        err_d = float(r_f.get('std_err_D', 0.0))
        r_large = float(r_f.get('mean_radius_large', 0.0))
        err_r_large = float(r_f.get('std_err_mean_radius_large', 0.0))
        r_large_str = f"{r_large:.0f} &plusmn; {err_r_large:.0f}" if err_r_large > 0 and has_l else f"{r_large:.0f}"
        std_large = round(float(r_f.get('std_radius_large', 0.0)), 0) if has_l else 0
        phi_large = float(r_f.get('total_volume_fraction_large', 0.0))
        phi_l_str = f"{phi_large:.2e}" if has_l else "0.00e+00"
        phi = float(r_f['phi'])

        html_content += f"""              <tr>
                <td>Frame {f:03d}</td>
                <td>{t_sec:.3f}</td>
                <td>{float(r_f['redchi']):.3f}</td>
                <td>{d_val:.3f} &plusmn; {err_d:.3f}</td>
                <td>&infin;</td>
                <td>{"True" if has_l else "False"}</td>
                <td>{r_large_str if has_l else "0"}</td>
                <td>{std_large if has_l else 0}</td>
                <td>{phi_l_str}</td>
                <td>{phi:.2e}</td>
              </tr>
"""

    html_content += r"""            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 6: PHYSICAL FORMALISM -->
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
          \(\sigma_{\mathrm{large}} / \langle R_{\mathrm{large}} \rangle = 0.065\):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{large}}(q) = (\Delta\rho)^2 \int_0^\infty n_{\mathrm{large}}(R) V^2(R) P(q, R) \, dR
          \]
        </div>

        <h3 style="font-size: 1.05rem; margin-top: 20px; color: #ffffff;">3. Viscoelastic Network Concentration Fluctuations</h3>
        <p style="color: var(--text-secondary); margin-top: 6px;">
          Dynamic viscoelastic network fluctuations explaining the high-\(q\) excess (Ornstein-Zernike / Tanaka dynamic network):
        </p>
        <div class="equation-box">
          \[
          I_{\mathrm{fluc}}(q) = \frac{A_{\mathrm{fluc}}(t)}{\left[1 + \frac{1}{6} (q \xi_{\mathrm{fluc}}(t))^2\right]^2}
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
      if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise();
      }
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
      document.getElementById('sd_frame_badge').innerText = 'Frame ' + padFrame(f) + ' (' + d.time + ' s)';
      document.getElementById('sd_img').src = '../short_distance/fits/fit_frame_' + padFrame(f) + '.png';

      document.getElementById('sd_chi2').innerText = d.chi2;
      document.getElementById('sd_hq_chi2').innerText = d.hq_chi2;
      document.getElementById('sd_Rmean').innerText = d.Rmean + ' nm';
      document.getElementById('sd_D').innerText = d.D;
      document.getElementById('sd_xi').innerText = d.xi === '∞' ? '∞' : (d.xi + ' nm');
      document.getElementById('sd_xi_fluc').innerText = d.xi_fluc + ' nm';
      document.getElementById('sd_A_fluc').innerText = d.A_fluc;
      document.getElementById('sd_ratio').innerText = d.ratio;
      document.getElementById('sd_phi').innerText = d.phi;
    }

    function stepShort(delta) {
      const s = document.getElementById('sd_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 350) val = 350;
      s.value = val;
      updateShortFrame(val);
    }

    function toggleShortPlay() {
      const btn = document.getElementById('sd_play_btn');
      if (sdTimer) {
        clearInterval(sdTimer);
        sdTimer = null;
        btn.innerText = 'Play';
        btn.classList.remove('btn-pause');
      } else {
        btn.innerText = 'Pause';
        btn.classList.add('btn-pause');
        sdTimer = setInterval(() => {
          const s = document.getElementById('sd_slider');
          let val = (parseInt(s.value) + 1);
          if (val > 350) val = 0;
          s.value = val;
          updateShortFrame(val);
        }, 120);
      }
    }

    function updateLongFrame(f) {
      f = parseInt(f);
      const d = longData[f];
      if (!d) return;

      document.getElementById('ld_frame_val').innerText = f;
      document.getElementById('ld_time_val').innerText = d.time + ' s';
      document.getElementById('ld_frame_badge').innerText = 'Frame ' + padFrame(f) + ' (' + d.time + ' s)';
      document.getElementById('ld_img').src = '../long_distance/fits/fit_frame_' + padFrame(f) + '.png';

      document.getElementById('ld_chi2').innerText = d.chi2;
      document.getElementById('ld_D').innerText = d.D;
      document.getElementById('ld_xi').innerText = d.xi;
      document.getElementById('ld_haslarge').innerText = d.has_large ? 'True' : 'False';
      document.getElementById('ld_rlarge').innerText = (d.R_large || 0) + ' nm';
      document.getElementById('ld_stdlarge').innerText = (d.std_large || 0) + ' nm (6.5%)';
      document.getElementById('ld_philarge').innerText = d.phi_large || '0.00e+00';
      document.getElementById('ld_phi').innerText = d.phi || '---';
    }

    function stepLong(delta) {
      const s = document.getElementById('ld_slider');
      let val = parseInt(s.value) + delta;
      if (val < 0) val = 0;
      if (val > 349) val = 349;
      s.value = val;
      updateLongFrame(val);
    }

    function toggleLongPlay() {
      const btn = document.getElementById('ld_play_btn');
      if (ldTimer) {
        clearInterval(ldTimer);
        ldTimer = null;
        btn.innerText = 'Play';
        btn.classList.remove('btn-pause');
      } else {
        btn.innerText = 'Pause';
        btn.classList.add('btn-pause');
        ldTimer = setInterval(() => {
          const s = document.getElementById('ld_slider');
          let val = (parseInt(s.value) + 1);
          if (val > 349) val = 0;
          s.value = val;
          updateLongFrame(val);
        }, 120);
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
    with open(dashboard_file, 'w', encoding='utf-8') as f_out:
        f_out.write(html_content)
    print(f"Saved master interactive dashboard: {dashboard_file}")


if __name__ == '__main__':
    build_dashboard()
