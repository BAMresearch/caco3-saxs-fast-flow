#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master Multi-Scale SAXS Synthesis Dashboard Generator
-----------------------------------------------------
Author: Tomasz Stawski, tomasz.stawski@bam.de

Generates results/index.html and results/dashboard/index.html providing
centralised navigation, cross-series parameter comparison, vector model
transformation diagram, Viscoelastic Phase Separation (VPS) verification,
Primary benchmark: 40 mM pure CaCO3 (ND-i-009).
"""

import os
import json
import shutil

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))

SERIES_METADATA = [
    {
        "id": "ND-i-003",
        "dir": "ND-i-003_10mM_pH11.3",
        "conc_mM": 10.0,
        "pH": 11.3,
        "additive": "None (pure CaCO₃)",
        "frames": 50,
        "time_range": "0.060 to 10.877 s",
        "regime": "Sub-spinodal nucleation",
        "dominant_morphology": "Discrete amorphous nanospheres + sheared droplets",
        "key_mechanisms": [
            "Below spinodal threshold: no mass fractal network (S(q) = 1)",
            "Classical nucleation of spherical ACC particles (<R_acc> = 138.7 nm, sigma = 21.0 nm)",
            "Intermediate-to-high q dominated by Porod q^-4 decay and form-factor oscillations",
            "Flat solvent background at 3.0 m camera geometry without colloidal excess"
        ]
    },
    {
        "id": "ND-i-004",
        "dir": "ND-i-004_15mM_pH11.3",
        "conc_mM": 15.0,
        "pH": 11.3,
        "additive": "None (pure CaCO₃)",
        "frames": 50,
        "time_range": "0.060 to 10.877 s",
        "regime": "Crossover phase boundary",
        "dominant_morphology": "Coexisting primary network nodes and discrete ACC nanospheres",
        "key_mechanisms": [
            "Boundary regime: simultaneous emergence of network and ACC particles",
            "Primary building blocks (R_prim = 4.94 nm) form network (Dm = 2.85, xi = 108.7 nm)",
            "Polydisperse ACC nanospheres (<R_acc> = 47.7 nm, phi_acc = 1.94e-4)",
            "Asymptotic Porod surface scattering from coarsening domains (A_porod = 1.34e-8)"
        ]
    },
    {
        "id": "ND-i-002",
        "dir": "ND-i-002_20mM_pH11.3",
        "conc_mM": 20.0,
        "pH": 11.3,
        "additive": "None (pure CaCO₃)",
        "frames": 50,
        "time_range": "0.060 to 10.877 s",
        "regime": "Moderate spinodal gelation",
        "dominant_morphology": "Percolated mass fractal network + Porod interface",
        "key_mechanisms": [
            "Spinodal liquid-liquid phase separation into mass fractal network",
            "Primary spherical building blocks (R_prim = 11.2 nm, sigma = 4.4 nm)",
            "Network scaling: D_f = 1.86, correlation cutoff xi = 356 nm",
            "Asymptotic Porod surface tail (A_porod = 5.28e-8) without discrete ACC form factor"
        ]
    },
    {
        "id": "ND-i-009",
        "dir": "ND-i-009_40mM_pH10.0",
        "conc_mM": 40.0,
        "pH": 10.0,
        "additive": "None (pure CaCO₃ primary case)",
        "frames": 351,
        "time_range": "0.045 to 5.295 s (15 ms/frame)",
        "regime": "Primary pure 40 mM benchmark",
        "dominant_morphology": "Continuous evolutionary fractal network + coarsening droplets + PNC mesh",
        "key_mechanisms": [
            "Primary pure CaCO₃ benchmark with 351 time-resolved frames (15 ms/frame)",
            "Time-resolved pre-percolation (xi = 46.7 nm -> inf) smoothly bridging into infinite gel",
            "Continuous point-to-point drift of fractal dimension D_f(t) in [1.881, 1.925]",
            "Viscoelastic fluctuation mesh (xi_fluc = 2.72 -> 2.84 nm) reflecting entangled PNC oligomers",
            "Secondary liquid droplets coarsening continuously (<R_large> = 2.27 -> 3.20 um)"
        ]
    },
    {
        "id": "ND-i-001",
        "dir": "ND-i-001_40mM_pH11.3",
        "conc_mM": 40.0,
        "pH": 11.3,
        "additive": "None (pure CaCO₃, pH 11.3)",
        "frames": 50,
        "time_range": "0.060 to 10.877 s",
        "regime": "Rapid quench spinodal demixing",
        "dominant_morphology": "Bicontinuous mass fractal gel + macroscopic emulsion droplets",
        "key_mechanisms": [
            "Pure CaCO₃ counterpart at pH 11.3 with instantaneous barrierless quench (t_ind < 0.06 s)",
            "Primary network units (R_prim = 3.50 nm, D_f = 2.33, xi = 735 nm)",
            "Capillary-pinched macroscopic droplets (<R_large> = 3.01 um, phi_large = 4.98e-5)",
            "Internal viscoelastic fluctuation mesh (xi_fluc = 1.34 nm)"
        ]
    },
    {
        "id": "ND-i-005",
        "dir": "ND-i-005_40mM_pH11.3_CoCl2",
        "conc_mM": 40.0,
        "pH": 11.3,
        "additive": "CoCl₂ (additive modification)",
        "frames": 50,
        "time_range": "0.060 to 10.877 s",
        "regime": "Additive perturbation branch",
        "dominant_morphology": "Extended fractal network + macroscopic droplets + elevated fluctuations",
        "key_mechanisms": [
            "Chemical dopant perturbation branch (only sample containing additive; all others are pure CaCO₃)",
            "Co²⁺ adsorption lowers interfacial energy gamma, enlarging primary nodes (R_prim = 6.59 nm)",
            "Extended correlation cutoff (xi = 1,190 nm) and open network geometry (D_f = 2.10)",
            "Macroscopic droplet population (<R_large> = 3.53 um, phi_large = 2.92e-5)"
        ]
    }
]


def load_all_parameters():
    data = {}
    for meta in SERIES_METADATA:
        json_path = os.path.join(RESULTS_DIR, meta["dir"], "output", "composite", "baseline_parameters.json")
        if os.path.exists(json_path):
            with open(json_path, "r", encoding="utf-8") as f:
                data[meta["id"]] = json.load(f)
        else:
            data[meta["id"]] = {}
    return data


def load_model_transformation_svg():
    svg_path = os.path.join(RESULTS_DIR, "output", "model_transformation_diagram.svg")
    if os.path.exists(svg_path):
        with open(svg_path, "r", encoding="utf-8") as f:
            return f.read()
    return "<p style='color: var(--danger);'>SVG diagram not found. Please run generate_model_transformation_diagram.py.</p>"


def load_trajectories_json():
    json_path = os.path.join(RESULTS_DIR, "output", "cross_series_trajectories.json")
    if os.path.exists(json_path):
        with open(json_path, "r", encoding="utf-8") as f:
            return f.read()
    return "{}"


def build_master_dashboard():
    params = load_all_parameters()
    svg_diagram = load_model_transformation_svg()
    trajectories_json_raw = load_trajectories_json()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CaCO3 multi-scale SAXS synthesis master dashboard</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  <style>
    :root {{
      --bg-primary: #0b0f19;
      --bg-secondary: #131b2e;
      --bg-card: #17223b;
      --bg-card-hover: #1e2c4a;
      --text-primary: #f1f5f9;
      --text-secondary: #94a3b8;
      --text-muted: #64748b;
      --accent: #38bdf8;
      --accent-hover: #0284c7;
      --border: #233354;
      --border-light: rgba(255, 255, 255, 0.08);
      --success: #10b981;
      --warning: #f59e0b;
      --danger: #ef4444;
      --code-bg: #070a12;
      --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      --font-mono: "JetBrains Mono", Consolas, Monaco, "Courier New", monospace;
    }}
    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}
    body {{
      font-family: var(--font-sans);
      background-color: var(--bg-primary);
      color: var(--text-primary);
      line-height: 1.55;
      padding-bottom: 80px;
    }}
    header {{
      background: linear-gradient(180deg, #131b2e 0%, #0b0f19 100%);
      border-bottom: 1px solid var(--border);
      padding: 32px 40px;
    }}
    .header-container {{
      max-width: 1480px;
      margin: 0 auto;
      display: flex;
      justify-content: space-between;
      align-items: flex-start;
      flex-wrap: wrap;
      gap: 20px;
    }}
    .header-title h1 {{
      font-size: 1.65rem;
      font-weight: 700;
      color: #ffffff;
      letter-spacing: -0.01em;
      margin-bottom: 6px;
    }}
    .header-title p {{
      font-size: 0.92rem;
      color: var(--text-secondary);
      max-width: 980px;
    }}
    .badge-bar {{
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
    }}
    .badge {{
      display: inline-block;
      padding: 4px 10px;
      font-size: 0.78rem;
      font-weight: 600;
      border-radius: 9999px;
      background-color: rgba(56, 189, 248, 0.12);
      color: var(--accent);
      border: 1px solid rgba(56, 189, 248, 0.28);
      font-family: var(--font-mono);
    }}
    .badge-success {{
      background-color: rgba(16, 185, 129, 0.12);
      color: var(--success);
      border-color: rgba(16, 185, 129, 0.28);
    }}
    nav.master-nav {{
      background-color: var(--bg-secondary);
      border-bottom: 1px solid var(--border);
      position: sticky;
      top: 0;
      z-index: 100;
      padding: 0 40px;
    }}
    .nav-container {{
      max-width: 1480px;
      margin: 0 auto;
      display: flex;
      gap: 8px;
      overflow-x: auto;
    }}
    .nav-btn {{
      background: none;
      border: none;
      color: var(--text-secondary);
      padding: 14px 20px;
      font-size: 0.88rem;
      font-weight: 600;
      cursor: pointer;
      border-bottom: 2px solid transparent;
      transition: all 0.15s ease;
      white-space: nowrap;
    }}
    .nav-btn:hover {{
      color: var(--text-primary);
    }}
    .nav-btn.active {{
      color: var(--accent);
      border-bottom-color: var(--accent);
    }}
    main {{
      max-width: 1480px;
      margin: 0 auto;
      padding: 32px 40px;
    }}
    .section-content {{
      display: none;
    }}
    .section-content.active {{
      display: block;
      animation: fadeIn 0.25s ease-in-out;
    }}
    @keyframes fadeIn {{
      from {{ opacity: 0; transform: translateY(4px); }}
      to {{ opacity: 1; transform: translateY(0); }}
    }}
    .card {{
      background-color: var(--bg-card);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 24px;
      margin-bottom: 24px;
    }}
    .card-title {{
      font-size: 1.15rem;
      font-weight: 600;
      color: #ffffff;
      margin-bottom: 14px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 10px;
    }}
    .grid-series {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(460px, 1fr));
      gap: 20px;
      margin-bottom: 24px;
    }}
    @media (max-width: 980px) {{
      .grid-series {{
        grid-template-columns: 1fr;
      }}
    }}
    .series-card {{
      background-color: var(--bg-card);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 20px;
      display: flex;
      flex-direction: column;
      justify-content: space-between;
      transition: border-color 0.15s ease, background-color 0.15s ease;
    }}
    .series-card:hover {{
      border-color: var(--accent);
      background-color: var(--bg-card-hover);
    }}
    .series-header {{
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      margin-bottom: 10px;
      border-bottom: 1px solid var(--border-light);
      padding-bottom: 10px;
    }}
    .series-id {{
      font-size: 1.15rem;
      font-weight: 700;
      color: #ffffff;
      font-family: var(--font-mono);
    }}
    .series-conditions {{
      font-size: 0.86rem;
      color: var(--accent);
      font-weight: 600;
    }}
    .series-regime {{
      font-size: 0.8rem;
      letter-spacing: 0.02em;
      color: var(--warning);
      margin-bottom: 10px;
      font-weight: 600;
    }}
    .series-mechanisms {{
      list-style-type: none;
      font-size: 0.84rem;
      color: var(--text-secondary);
      margin-bottom: 16px;
    }}
    .series-mechanisms li {{
      margin-bottom: 6px;
      padding-left: 16px;
      position: relative;
    }}
    .series-mechanisms li::before {{
      content: "•";
      color: var(--accent);
      position: absolute;
      left: 0;
      font-weight: bold;
    }}
    .series-metrics {{
      background-color: var(--code-bg);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 10px 14px;
      font-family: var(--font-mono);
      font-size: 0.8rem;
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 6px 12px;
      margin-bottom: 16px;
    }}
    .metric-item {{
      display: flex;
      justify-content: space-between;
    }}
    .metric-label {{
      color: var(--text-muted);
    }}
    .metric-value {{
      color: #ffffff;
      font-weight: 600;
    }}
    .btn-dashboard {{
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
      padding: 10px 16px;
      border-radius: 6px;
      background-color: rgba(56, 189, 248, 0.15);
      color: var(--accent);
      border: 1px solid rgba(56, 189, 248, 0.35);
      font-weight: 600;
      font-size: 0.86rem;
      text-decoration: none;
      transition: all 0.15s ease;
      cursor: pointer;
    }}
    .btn-dashboard:hover {{
      background-color: var(--accent);
      color: #0b0f19;
      border-color: var(--accent);
    }}
    .table-container {{
      overflow-x: auto;
      border: 1px solid var(--border);
      border-radius: 8px;
      margin-top: 12px;
    }}
    table.data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.84rem;
    }}
    table.data-table th {{
      background-color: #19253d;
      color: var(--text-primary);
      padding: 10px 12px;
      text-align: left;
      font-weight: 600;
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
    }}
    table.data-table th.num, table.data-table td.num {{
      text-align: right;
      font-family: var(--font-mono);
    }}
    table.data-table td {{
      padding: 9px 12px;
      border-bottom: 1px solid var(--border-light);
      color: var(--text-secondary);
      white-space: nowrap;
    }}
    table.data-table tr:hover td {{
      background-color: rgba(56, 189, 248, 0.04);
      color: #ffffff;
    }}
    .text-block {{
      color: var(--text-secondary);
      font-size: 0.92rem;
      line-height: 1.65;
      margin-bottom: 16px;
    }}
    .text-block strong {{
      color: #ffffff;
    }}
    .text-block h3 {{
      color: #ffffff;
      font-size: 1.02rem;
      font-weight: 600;
      margin: 22px 0 10px 0;
      border-left: 3px solid var(--accent);
      padding-left: 10px;
    }}
    .equation-box {{
      background-color: var(--code-bg);
      border-left: 3px solid var(--accent);
      border-radius: 4px;
      padding: 16px 20px;
      margin: 16px 0;
      overflow-x: auto;
      font-size: 0.95rem;
    }}
    .gallery-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(460px, 1fr));
      gap: 20px;
    }}
    .gallery-item {{
      background-color: #0b1120;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px;
      display: flex;
      flex-direction: column;
    }}
    .gallery-caption {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 10px;
      font-size: 0.86rem;
      font-weight: 600;
      color: #ffffff;
    }}
    .gallery-img {{
      width: 100%;
      height: auto;
      border-radius: 4px;
      border: 1px solid var(--border-light);
    }}
    /* Trajectory Explorer Styling */
    .controls-bar {{
      display: flex;
      flex-wrap: wrap;
      gap: 16px;
      align-items: center;
      justify-content: space-between;
      background-color: #111a2e;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px 18px;
      margin-bottom: 16px;
    }}
    .control-group {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}
    .control-label {{
      font-size: 0.8rem;
      font-weight: 600;
      color: var(--text-muted);
      letter-spacing: 0.02em;
    }}
    .toggle-btn {{
      background-color: #17223b;
      color: var(--text-secondary);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 6px 12px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.15s ease;
    }}
    .toggle-btn:hover {{
      border-color: var(--accent);
      color: #ffffff;
    }}
    .toggle-btn.active {{
      background-color: var(--accent);
      color: #0b0f19;
      border-color: var(--accent);
    }}
    .series-checkbox-group {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: center;
    }}
    .series-checkbox-item {{
      display: inline-flex;
      align-items: center;
      gap: 6px;
      font-size: 0.82rem;
      color: var(--text-secondary);
      cursor: pointer;
      user-select: none;
    }}
    .series-checkbox-item input {{
      cursor: pointer;
      accent-color: var(--accent);
    }}
    .color-dot {{
      width: 10px;
      height: 10px;
      border-radius: 50%;
      display: inline-block;
    }}
    .chart-wrapper {{
      position: relative;
      width: 100%;
      height: 480px;
      background-color: #070a12;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 14px;
      margin-bottom: 20px;
    }}
    .svg-container {{
      width: 100%;
      overflow-x: auto;
      background-color: #0b1120;
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 12px;
      margin: 18px 0;
    }}
    .svg-container svg {{
      display: block;
      width: 100%;
      min-width: 1120px;
      height: auto;
    }}
    .pnc-box {{
      background-color: rgba(37, 99, 235, 0.08);
      border: 1px solid rgba(56, 189, 248, 0.25);
      border-left: 4px solid var(--accent);
      border-radius: 6px;
      padding: 16px 20px;
      margin: 18px 0;
    }}
    .pnc-box h4 {{
      color: var(--accent);
      font-size: 0.95rem;
      font-weight: 600;
      margin-bottom: 8px;
    }}
  </style>
</head>
<body>

  <header>
    <div class="header-container">
      <div class="header-title">
        <h1>CaCO<sub>3</sub> multi-scale SAXS synthesis master dashboard</h1>
        <p>
          Technical cross-series synthesis of time-resolved synchrotron SAXS data (ESRF ID02, 3.0 m and 31.0 m camera geometries,
          q: 0.002 to 2.5 nm<sup>-1</sup>). Evaluates non-classical nucleation, spinodal liquid-liquid phase separation,
          and macroscopic emulsion coarsening across concentration, pH, and additive variations.
        </p>
        <div style="font-size: 0.85rem; color: #94a3b8; margin-top: 6px;">
          Author: Tomasz Stawski (tomasz.stawski@bam.de), Federal Institute for Materials Research and Testing (BAM)
        </div>
      </div>
      <div class="badge-bar">
        <span class="badge badge-success">6 reaction series</span>
        <span class="badge">10 – 40 mM CaCO<sub>3</sub></span>
        <span class="badge">Pure CaCO<sub>3</sub> primary benchmark</span>
        <span class="badge">&plusmn;1&sigma; covariance</span>
      </div>
    </div>
  </header>

  <nav class="master-nav">
    <div class="nav-container">
      <button class="nav-btn active" onclick="switchSection('hub')">1. Dataset hub and navigation</button>
      <button class="nav-btn" onclick="switchSection('evolution')">2. Physicochemical evolution and model continuum</button>
      <button class="nav-btn" onclick="switchSection('trajectories')">3. Cross-series trajectory synthesis and kinetic comparison</button>
      <button class="nav-btn" onclick="switchSection('parameters')">4. Comparative parameter matrix</button>
      <button class="nav-btn" onclick="switchSection('gallery')">5. Composite scattering curves</button>
      <button class="nav-btn" onclick="switchSection('formalism')">6. Physical scattering formalism</button>
    </div>
  </nav>

  <main>
    <!-- SECTION 1: HUB & DIRECT ACCESS CARDS -->
    <div id="sec-hub" class="section-content active">
      <div class="card">
        <div class="card-title">
          <span>Reaction series overview and direct dashboard links</span>
          <span style="font-size: 0.82rem; color: var(--text-muted); font-weight: normal;">Select a series to open its independent 6-tab time-resolved dashboard</span>
        </div>
        <p class="text-block" style="margin-bottom: 20px;">
          Each series contains an independent 6-tab dashboard with interactive frame-by-frame inspectors, live parameter telemetry,
          uncertainty envelopes (&plusmn;1&sigma;), and full data tables. The 40 mM pure CaCO₃ dataset (ND-i-009) serves as the primary benchmark.
        </p>

        <div class="grid-series">
"""

    for meta in SERIES_METADATA:
        sid = meta["id"]
        p = params.get(sid, {})
        chi2 = p.get("redchi") or p.get("chi2_red", 0.0)

        r_prim = p.get("Rmean") or p.get("R_prim", 0.0)
        df_val = p.get("D") or p.get("Dm", 0.0)
        xi_val = p.get("xi", 0.0)
        r_large = p.get("mean_radius_large", 0.0)
        r_acc = p.get("R_acc", 0.0)

        xi_str = f"{xi_val:.0f} nm" if 0 < xi_val < 1e5 else ("inf" if xi_val >= 1e5 else "N/A")
        r_prim_str = f"{r_prim:.2f} nm" if r_prim > 0 else "N/A"
        df_str = f"{df_val:.3f}" if df_val > 0 else "N/A"
        drop_str = f"{r_large:.0f} nm" if r_large > 0 else (f"{r_acc:.1f} nm" if r_acc > 0 else "N/A")
        drop_label = "R_acc:" if r_acc > 0 and r_large == 0 else "R_large:"

        html += f"""          <div class="series-card">
            <div>
              <div class="series-header">
                <span class="series-id">{sid}</span>
                <span class="series-conditions">{meta['conc_mM']:.0f} mM | {meta['pH']} | {meta['additive']}</span>
              </div>
              <div class="series-regime">{meta['regime']}</div>
              <ul class="series-mechanisms">
"""
        for m in meta["key_mechanisms"]:
            html += f"                <li>{m}</li>\n"

        html += f"""              </ul>
              <div class="series-metrics">
                <div class="metric-item"><span class="metric-label">Chi2_red:</span><span class="metric-value" style="color: var(--success);">{chi2:.3f}</span></div>
                <div class="metric-item"><span class="metric-label">R_prim:</span><span class="metric-value">{r_prim_str}</span></div>
                <div class="metric-item"><span class="metric-label">D_fractal:</span><span class="metric-value">{df_str}</span></div>
                <div class="metric-item"><span class="metric-label">Cutoff xi:</span><span class="metric-value">{xi_str}</span></div>
                <div class="metric-item"><span class="metric-label">{drop_label}</span><span class="metric-value">{drop_str}</span></div>
                <div class="metric-item"><span class="metric-label">Frames:</span><span class="metric-value">{meta['frames']}</span></div>
              </div>
            </div>
            <a class="btn-dashboard" href="./{meta['dir']}/output/dashboard/index.html">
              Open {sid} dashboard &rarr;
            </a>
          </div>
"""

    html += f"""        </div>

        <!-- DATASET EXCLUSIONS NOTE -->
        <div style="margin-top: 1.5rem; padding: 1rem 1.25rem; background: rgba(239, 68, 68, 0.08); border: 1px solid rgba(239, 68, 68, 0.25); border-radius: 8px;">
          <div style="font-weight: 600; color: #ef4444; font-size: 0.95rem; margin-bottom: 0.35rem;">
            Excluded Experimental Series (Corrupted / Broken / Unused Datasets)
          </div>
          <div style="color: var(--text-secondary); font-size: 0.85rem; line-height: 1.5;">
            <strong>ND-i-006</strong> (40.0 mM, pH 9.75), <strong>ND-i-008</strong> (40.0 mM, pH 9.5), and <strong>ND-i-007</strong> (40.0 mM, pH 9.0) are strictly excluded from comparative and kinetic analysis suites.
          </div>
        </div>
      </div>
    </div>

    <!-- SECTION 2: PHYSICOCHEMICAL EVOLUTION & MODEL CONTINUUM -->
    <div id="sec-evolution" class="section-content">
      <div class="card">
        <div class="card-title">
          <span>Model transformation across supersaturation and phase space</span>
          <span class="badge badge-success">Data-driven Vector SVG</span>
        </div>
        
        <div class="text-block">
          The reaction series systematically explore calcium carbonate precipitation across reactant concentration, 
          solution pH, and chemical additives. All samples except ND-i-005 are pure CaCO₃. The primary 40 mM pure CaCO₃ benchmark 
          is <strong>ND-i-009</strong>, which resolves continuous time-resolved percolation over 351 frames (15 ms/frame).
          The vector diagram below is built directly around the full-range composite experimental scattering curves (\\(q \\in [0.002, 2.5]\\text{{ nm}}^{{-1}}\\)),
          tracking how physical morphologies and mathematical model components evolve across the spinodal phase boundary with unique, standardised component colours.
        </div>

        <!-- STANDARDISED PALETTE BANNER -->
        <div style="margin: 0.75rem 0 1.25rem 0; padding: 0.75rem 1rem; background: #070d19; border: 1px solid #1e293b; border-radius: 8px; display: flex; flex-wrap: wrap; gap: 1.25rem; font-size: 0.82rem; align-items: center;">
          <span style="font-weight: 700; color: #cbd5e1; text-transform: uppercase; font-size: 0.75rem; letter-spacing: 0.05em;">Standardised Model Components:</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #94a3b8;"><span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#64748b;"></span> Experimental Data</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #f87171;"><span style="display:inline-block; width:14px; height:3px; background:#dc2626;"></span> Total Unified Fit</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #34d399;"><span style="display:inline-block; width:14px; height:3px; background:#059669;"></span> Spinodal Fractal Network</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #38bdf8;"><span style="display:inline-block; width:14px; height:3px; background:#0284c7;"></span> Discrete ACC Spheres</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #c084fc;"><span style="display:inline-block; width:14px; height:3px; background:#9333ea;"></span> Tanaka PNC Fluctuation Mesh</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #fbbf24;"><span style="display:inline-block; width:14px; height:3px; background:#d97706;"></span> Macroscopic Droplets / Porod Tail</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #94a3b8;"><span style="display:inline-block; width:14px; height:2px; background:#64748b;"></span> Solvent Background</span>
        </div>

        <!-- EMBEDDED INLINE SVG DIAGRAM -->
        <div class="svg-container">
          {svg_diagram}
        </div>

        <!-- DEDICATED PNC AND VISCOELASTIC MESH DISCUSSION -->
        <div class="pnc-box">
          <h4>Role of the viscoelastic mesh in the context of prenucleation clusters (PNCs)</h4>
          <p class="text-block" style="margin-bottom: 8px;">
            In the non-classical nucleation framework, hydrated calcium and carbonate ions form dynamic coordination oligomers or 
            prenucleation clusters (PNCs) with characteristic dimensions \\(\\sim 0.6 - 2.0\\text{{ nm}}\\). In dilute sub-spinodal solutions, 
            these species remain in dynamic equilibrium with free ions without forming small-angle scattering excess.
          </p>
          <p class="text-block" style="margin-bottom: 8px;">
            Upon crossing the spinodal boundary into liquid-liquid phase separation (LLPS), the dense hydrated liquid precursor phase 
            rapidly concentrates the ions (local concentration \\(\\sim 5 - 10\\text{{ M}}\\)). Under these concentrated conditions, 
            PNCs do not persist as isolated globular particles with sharp boundaries; instead, their hydration shells overlap and 
            interpenetrate, forming a <strong>transient, dynamically entangled network or viscoelastic mesh</strong>.
          </p>
          <p class="text-block" style="margin-bottom: 0;">
            The dynamic concentration fluctuation component in the SAXS model:
            \\[
            I_{{\\rm fluc}}(q) = \\frac{{A_{{\\rm fluc}}}}{{\\left[ 1 + \\frac{{1}}{{6}} (q \\xi_{{\\rm fluc}})^2 \\right]^2}}
            \\]
            directly measures the <strong>mesh size (mean coordination distance between entangled PNC nodes)</strong> 
            \\(\\xi_{{\\rm fluc}} \\approx 1.5 - 3.0\\text{{ nm}}\\). The Tanaka decay exponent \\(\\mu = 4.0\\) reflects continuous 
            concentration gradients without sharp dielectric phase boundaries. As the system ages and expels water, this viscoelastic mesh 
            structurally arrests and vitrifies into the primary amorphous calcium carbonate building blocks 
            (\\(\\langle R_{{\\rm prim}} \\rangle \\approx 2.8 - 6.5\\text{{ nm}}\\)) that constitute the mesoscale mass fractal gel.
          </p>
        </div>

        <div class="text-block">
          <h3>Concentration axis in pure CaCO₃ (10 &rarr; 15 &rarr; 20 &rarr; 40 mM)</h3>
          <p>
            Varying reactant concentration in pure CaCO₃ establishes a sharp boundary at <strong>C* &approx; 12 – 15 mM</strong>:
          </p>
          <ul style="margin: 8px 0 12px 20px;">
            <li>
              <strong>10.0 mM (ND-i-003, sub-spinodal regime):</strong>
              Precipitation proceeds via classical nucleation and growth of discrete, non-percolating amorphous calcium carbonate (ACC) nanospheres.
              Small-angle scattering displays no mass fractal structure factor (<em>S(q) &equiv; 1</em>).
              Scattering at <em>q &in; [0.005, 0.15] nm<sup>-1</sup></em> is strictly governed by the polydisperse sphere form factor with asymptotic Porod <em>q<sup>-4</sup></em> decay (<em>D = 3</em>, sharp interface).
              The short-distance camera (3.0 m) detects flat solvent background without colloidal building block excess.
            </li>
            <li>
              <strong>15.0 mM (ND-i-004, crossover phase boundary):</strong>
              At 15 mM, the system resides directly on the spinodal phase boundary. 
              The scattering profile requires a hybrid model: simultaneous presence of <strong>primary building blocks</strong> (&lang;R<sub>prim</sub>&rang; &approx; 4.94 nm) assembling into an early percolated mass fractal network (<em>D<sub>m</sub> = 2.85</em>, &xi; = 108.7 nm) 
              coexisting with a population of <strong>discrete spherical ACC nanospheres</strong> (&lang;R<sub>acc</sub>&rang; &approx; 47.7 nm) generating form-factor oscillations in the 3.0 m camera.
            </li>
            <li>
              <strong>20.0 mM (ND-i-002, spinodal gelation):</strong>
              Above the boundary, spinodal decomposition dominates. Primary spherical nodes (&lang;R<sub>prim</sub>&rang; &approx; 11.2 nm) link into a space-filling mass fractal network with <em>D<sub>f</sub> = 1.86</em> and correlation length &xi; = 356 nm.
              Scattering at lowest <em>q</em> (&lt; 0.005 nm<sup>-1</sup>) transitions to Porod surface scattering (<em>A</em><sub>porod</sub> / <em>q</em><sup>4</sup>), indicating the onset of coarsening dense domains. Discrete ACC spheres are absent.
            </li>
            <li>
              <strong>40.0 mM (ND-i-009, primary pure CaCO₃ benchmark):</strong>
              Primary benchmark with 351 frames (15 ms/frame). Resolves the pre-percolation induction phase (frame 0: &xi; = 46.7 nm) 
              smoothly transitioning into an infinite percolated network (&xi; &rarr; &infin;) with continuous fractal dimension drift 
              <em>D<sub>f</sub>(t) &in; [1.881, 1.925]</em>. Simultaneously resolves the PNC viscoelastic fluctuation mesh (&xi;<sub>fluc</sub> &approx; 2.72 &rarr; 2.84 nm) 
              and secondary droplet coarsening (2.27 &rarr; 3.20 &mu;m).
            </li>
            <li>
              <strong>40.0 mM at pH 11.3 (ND-i-001, rapid quench counterpart):</strong>
              At elevated pH 11.3, higher free carbonate activity accelerates demixing, resulting in an instantaneous barrierless quench 
              (<em>t</em><sub>ind</sub> &lt; 0.06 s) into a dense bicontinuous network (<em>D<sub>f</sub> = 2.33</em>, &xi; = 735 nm, 
              &lang;R<sub>prim</sub>&rang; = 3.50 nm) with pinched droplets (&lang;R<sub>large</sub>&rang; = 3.01 &mu;m).
            </li>
          </ul>
        </div>

        <div class="text-block">
          <h3>Additive perturbation branch (Co²⁺ dopant, ND-i-005)</h3>
          <p>
            ND-i-005 is the sole reaction series containing a foreign chemical additive. 
            The addition of cobalt(II) modifies interfacial thermodynamics:
          </p>
          <ul style="margin: 8px 0 12px 20px;">
            <li>
              <strong>Enlargement of primary nodes:</strong> Cobalt adsorption lowers interfacial tension &gamma;, shifting the primary building block radius to 6.59 nm.
            </li>
            <li>
              <strong>Extended correlation cutoff:</strong> Network correlation cutoff expands to &xi; = 1,190 nm with a more open fractal geometry (<em>D<sub>f</sub> = 2.10</em>).
            </li>
            <li>
              <strong>Fluctuation enhancement:</strong> Elevated viscoelastic fluctuation amplitude confirms altered hydration dynamics within the dense coordination network.
            </li>
          </ul>
        </div>

        <!-- VISCOELASTIC PHASE SEPARATION (VPS) VERIFICATION SUBSECTION -->
        <div class="text-block" style="margin-top: 24px; border-top: 1px solid var(--border); padding-top: 16px;">
          <h3>Viscoelastic Phase Separation (VPS) framework verification</h3>
          <p>
            The experimental trajectories and model transitions verify the mathematical framework of 
            <strong>Viscoelastic Phase Separation (VPS)</strong> formulated by Tanaka (2000, 2001) and detailed in 
            <em>Mathematical Framework for VPS Scattering</em>:
          </p>
          <ul style="margin: 8px 0 12px 20px;">
            <li>
              <strong>Dynamic asymmetry (&alpha; &gt;&gt; 0):</strong> Hydrated calcium-carbonate ion oligomers and DOLLOPs exhibit slow relaxation timescales (&tau;<sub>p</sub> &gt;&gt; &tau;<sub>s</sub>) relative to the bulk water solvent, generating strong viscoelastic stress coupling during phase separation.
            </li>
            <li>
              <strong>Viscoelastic length scale (&xi;<sub>ve</sub>):</strong> The crossover scale &xi;<sub>ve</sub> &sim; &radic;(K + 4/3 G)/&zeta; &tau; separates elastic suppression at short wavelengths (q &xi;<sub>ve</sub> &gt; 1) from hydrodynamic demixing at large wavelengths (q &xi;<sub>ve</sub> &lt; 1).
            </li>
            <li>
              <strong>Volume shrinking of the transient gel:</strong> In the 40 mM pure benchmark (ND-i-009), the network strands densify continuously (D<sub>m</sub> drifts from 1.881 to 1.925) whilst primary unit radius remains invariant (R<sub>prim</sub> &approx; 2.80 nm). This quantitative result confirms solvent expulsion from the viscoelastic network.
            </li>
            <li>
              <strong>Phase inversion and droplet coexistence:</strong> Mechanical stress relaxation along contracting struts leads to localized rupture and capillary pinch-off, shedding compact macroscopic liquid droplets (&lang;R<sub>large</sub>&rang; &approx; 2.2 to 3.5 &mu;m, &phi;<sub>large</sub> &approx; 10<sup>-4</sup>) that coexist with the space-filling gel mesh.
            </li>
            <li>
              <strong>Distinction from Arrested Spinodal Decomposition (ASD):</strong> In contrast to ASD in hard colloidal suspensions where structure factors freeze permanently via glass jamming, the calcium carbonate system shows continuous kinetic evolution, volume shrinking, and phase inversion.
            </li>
          </ul>

          <div style="margin-top: 16px;">
            <h4 style="color: #ffffff; font-size: 0.92rem; margin-bottom: 8px; font-weight: 600;">
              VPS framework verification synthesis figure
            </h4>
            <img src="./output/vps_framework_verification.png" alt="Viscoelastic Phase Separation verification" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
            <div style="font-size: 0.82rem; color: var(--text-muted); margin-top: 6px;">
              Panels: (a) Volume shrinking in ND-i-009 (fractal dimension D<sub>m</sub>(t) and volume fraction &phi;(t)).
              (b) Multi-scale spatial hierarchy across series.
              (c) Model transformation continuum across supersaturation.
              (d) Tanaka VPS scattering signature verification.
            </div>
          </div>
        </div>
      </div>
    </div>

    <!-- SECTION 3: CROSS-SERIES TRAJECTORY SYNTHESIS -->
    <div id="sec-trajectories" class="section-content">
      <div class="card">
        <div class="card-title">
          <span>Cross-series trajectory synthesis and kinetic comparison</span>
          <span class="badge badge-success">Time-resolved evolution</span>
        </div>

        <div class="text-block">
          To compare reaction kinetics across varying concentrations and experimental conditions, dense phase volume fractions 
          are normalized against their theoretical stoichiometric limit:
          \\[
          \\phi_{{\\rm max}} = \\frac{{[\\mathrm{{Ca}}^{{2+}}]_0 \\cdot M_{{\\mathrm{{CaCO}}_3}}}}{{\\rho_{{\\rm ACC}} \\cdot 1000}} \\approx 5.0 \\times 10^{{-4}} \\times \\left(\\frac{{C_{{\\rm mM}}}}{{10}}\\right)
          \\]
          yielding the standardized reaction extent \\(\\alpha(t) = \\phi(t) / \\phi_{{\\rm max}}\\).
        </div>

        <!-- STATIC MULTI-PANEL PUBLICATION FIGURE -->
        <div style="margin-bottom: 24px;">
          <h3 style="color: #ffffff; font-size: 0.98rem; margin-bottom: 10px; font-weight: 600;">
            Multi-panel kinetic synthesis and structural bifurcation
          </h3>
          <img src="./output/cross_series_trajectories_comparison.png" alt="Cross-series trajectory comparison" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);">
          <div style="font-size: 0.82rem; color: var(--text-muted); margin-top: 6px;">
            Panels: (a) Reaction extent &alpha;(t) displaying classical induction delay at 10 mM vs instantaneous spinodal quench at &ge;20 mM.
            (b) Divergence of characteristic meso-scale length scale (ACC radius at 10 mM vs network cutoff &xi; at &ge;15 mM).
            (c) Structural bifurcation across concentration in pure CaCO₃ with relative phase dominance inset.
            (d) Stationary emulsion droplet population across series.
          </div>
        </div>

        <!-- COMPONENT DOMINANCE & PHASE INVERSION SYNTHESIS CARD -->
        <div style="margin-bottom: 24px; padding: 1.25rem; background: #070d19; border: 1px solid #1e293b; border-left: 4px solid var(--accent); border-radius: 8px;">
          <h4 style="color: var(--accent); font-size: 0.95rem; margin-bottom: 8px; font-weight: 700;">
            Relative component proportions and dynamic phase inversion (including Tanaka PNC mesh)
          </h4>
          <p class="text-block" style="margin-bottom: 12px; font-size: 0.88rem;">
            When models with multiple coexistence mechanisms are fitted across the series, the relative phase proportions quantify 
            a systematic morphological regime transition from <strong>discrete ACC sphere dominance</strong> to <strong>mesoscale spinodal gel network dominance</strong> 
            coexisting with the <strong>Tanaka viscoelastic PNC fluctuation mesh</strong> and <strong>stationary capillary droplets</strong>:
          </p>
          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 1rem;">
            <!-- 10 mM Card -->
            <div style="background: rgba(2, 132, 199, 0.08); border: 1px solid rgba(2, 132, 199, 0.25); border-radius: 6px; padding: 0.85rem;">
              <div style="color: #38bdf8; font-weight: 700; font-size: 0.8rem; text-transform: uppercase;">10 mM (ND-i-003, Sub-spinodal)</div>
              <div style="font-size: 1.05rem; font-weight: 700; color: #fff; margin: 0.2rem 0;">100% Discrete ACC Spheres</div>
              <!-- Segmented Bar -->
              <div style="display: flex; height: 6px; border-radius: 3px; overflow: hidden; margin: 6px 0 8px 0; background: #1e293b;">
                <div style="width: 28%; background: #0284c7;" title="ACC Mode 1: 28%"></div>
                <div style="width: 72%; background: #0ea5e9;" title="ACC Mode 2: 72%"></div>
              </div>
              <div style="font-size: 0.8rem; color: var(--text-muted);">
                Mode 1 (28% @ 21 nm) + Mode 2 (72% @ 60 nm). <strong>0% Fractal Gel | 0% Tanaka PNC Mesh</strong>. Sub-spinodal non-percolating solution where PNC fluctuations remain below the small-angle scattering threshold.
              </div>
            </div>

            <!-- 15 mM Card -->
            <div style="background: rgba(45, 212, 191, 0.08); border: 1px solid rgba(45, 212, 191, 0.25); border-radius: 6px; padding: 0.85rem;">
              <div style="color: #2dd4bf; font-weight: 700; font-size: 0.8rem; text-transform: uppercase;">15 mM (ND-i-004, Dynamic Inversion)</div>
              <div style="font-size: 1.05rem; font-weight: 700; color: #fff; margin: 0.2rem 0;">96% ACC &rarr; 76% Fractal Gel</div>
              <!-- Segmented Bar -->
              <div style="display: flex; height: 6px; border-radius: 3px; overflow: hidden; margin: 6px 0 8px 0; background: #1e293b;">
                <div style="width: 63%; background: #0284c7;" title="ACC Spheres: 63%"></div>
                <div style="width: 37%; background: #059669;" title="Fractal Gel: 37%"></div>
              </div>
              <div style="font-size: 0.8rem; color: var(--text-muted);">
                Dynamic crossover at <em>t</em> &approx; 0.15 to 0.25 s: initial sphere nucleation arrests, and spinodal gel network develops to comprise the predominant mass fraction. Baseline composite: <strong>63% ACC Spheres | 37% Fractal Gel</strong>.
              </div>
            </div>

            <!-- 20 mM Card -->
            <div style="background: rgba(5, 150, 105, 0.08); border: 1px solid rgba(5, 150, 105, 0.25); border-radius: 6px; padding: 0.85rem;">
              <div style="color: #34d399; font-weight: 700; font-size: 0.8rem; text-transform: uppercase;">20 mM (ND-i-002, Spinodal Gel)</div>
              <div style="font-size: 1.05rem; font-weight: 700; color: #fff; margin: 0.2rem 0;">76% Fractal | 24% Tanaka PNC</div>
              <!-- Segmented Bar -->
              <div style="display: flex; height: 6px; border-radius: 3px; overflow: hidden; margin: 6px 0 8px 0; background: #1e293b;">
                <div style="width: 76%; background: #059669;" title="Fractal Network: 76%"></div>
                <div style="width: 24%; background: #9333ea;" title="Tanaka PNC Mesh: 24%"></div>
              </div>
              <div style="font-size: 0.8rem; color: var(--text-muted);">
                Discrete ACC spheres suppressed (<strong>0%</strong>). Continuous percolated gel network (<em>D<sub>f</sub> = 1.86</em>, &xi; = 356 nm) coexists with hydrated <strong>Tanaka PNC mesh</strong> (&lang;R<sub>pnc</sub>&rang; = 1.46 nm, &phi;<sub>pnc</sub> = 1.50 &times; 10<sup>-5</sup>).
              </div>
            </div>

            <!-- 40 mM Card -->
            <div style="background: rgba(248, 113, 113, 0.08); border: 1px solid rgba(248, 113, 113, 0.25); border-radius: 6px; padding: 0.85rem;">
              <div style="color: #f87171; font-weight: 700; font-size: 0.8rem; text-transform: uppercase;">40 mM (ND-i-009 &amp; ND-i-001, Benchmarks)</div>
              <div style="font-size: 1.05rem; font-weight: 700; color: #fff; margin: 0.2rem 0;">73% Fractal | 26% PNC | 1% Droplets</div>
              <!-- Segmented Bar -->
              <div style="display: flex; height: 6px; border-radius: 3px; overflow: hidden; margin: 6px 0 8px 0; background: #1e293b;">
                <div style="width: 73%; background: #059669;" title="Fractal Gel: 73%"></div>
                <div style="width: 26%; background: #9333ea;" title="Tanaka PNC Mesh: 26%"></div>
                <div style="width: 1%; background: #d97706;" title="Droplets: 1%"></div>
              </div>
              <div style="font-size: 0.8rem; color: var(--text-muted);">
                Discrete spheres <strong>0%</strong>. Pure benchmark ND-i-009: 73% fractal gel, 26% Tanaka PNC mesh (&xi;<sub>fluc</sub> &approx; 2.7 nm), 1% droplets. Rapid quench ND-i-001 captures <strong>87.1%</strong> initial PNC mesh before complete vitrification.
              </div>
            </div>
          </div>

          <!-- COMPREHENSIVE COMPONENT TALLY TABLE -->
          <div style="margin-top: 18px;">
            <div style="font-size: 0.85rem; font-weight: 700; color: #cbd5e1; margin-bottom: 8px; text-transform: uppercase; letter-spacing: 0.04em;">
              Complete multi-component volume and invariant tally across reaction series
            </div>
            <div class="table-container">
              <table class="data-table" style="font-size: 0.82rem;">
                <thead>
                  <tr>
                    <th>Series</th>
                    <th>Conditions</th>
                    <th class="num" style="color: #38bdf8;">Discrete ACC Spheres (%)</th>
                    <th class="num" style="color: #34d399;">Fractal Network (%)</th>
                    <th class="num" style="color: #c084fc;">Tanaka PNC Mesh (%)</th>
                    <th class="num" style="color: #fbbf24;">Droplets / Porod (%)</th>
                    <th>Dominant Mechanism &amp; Morphological Notes</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>ND-i-003</strong></td>
                    <td>10 mM, pH 11.3, pure</td>
                    <td class="num" style="font-weight: 700; color: #38bdf8;">100.0%</td>
                    <td class="num">0.0%</td>
                    <td class="num">0.0%</td>
                    <td class="num">0.0% (trace)</td>
                    <td>Classical nucleation of polydisperse nanospheres (Mode 1: 28.4% + Mode 2: 71.6%); no percolated mesh</td>
                  </tr>
                  <tr>
                    <td><strong>ND-i-004</strong></td>
                    <td>15 mM, pH 11.3, pure</td>
                    <td class="num" style="font-weight: 700; color: #38bdf8;">63.1%</td>
                    <td class="num" style="font-weight: 700; color: #34d399;">36.9%</td>
                    <td class="num">0.0% (transient)</td>
                    <td class="num">0.0% (trace)</td>
                    <td>Spinodal crossover boundary; dynamic inversion (95.8% &rarr; 23.7% ACC, 4.2% &rarr; 76.3% gel at <em>t</em> &approx; 0.2 s)</td>
                  </tr>
                  <tr>
                    <td><strong>ND-i-002</strong></td>
                    <td>20 mM, pH 11.3, pure</td>
                    <td class="num">0.0%</td>
                    <td class="num" style="font-weight: 700; color: #34d399;">75.9%</td>
                    <td class="num" style="font-weight: 700; color: #c084fc;">24.1%</td>
                    <td class="num">0.0% (trace Porod)</td>
                    <td>Continuous percolated gel network coexisting with Tanaka viscoelastic PNC mesh (R<sub>pnc</sub> = 1.46 nm)</td>
                  </tr>
                  <tr style="background-color: rgba(234, 88, 12, 0.08);">
                    <td><strong>ND-i-009</strong></td>
                    <td>40 mM, pH 10.0, pure</td>
                    <td class="num">0.0%</td>
                    <td class="num" style="font-weight: 700; color: #34d399;">72.9%</td>
                    <td class="num" style="font-weight: 700; color: #c084fc;">25.9%</td>
                    <td class="num" style="font-weight: 700; color: #fbbf24;">1.2%</td>
                    <td>Primary pure benchmark: pre-percolation, volume shrinking (<em>D<sub>m</sub></em> 1.88&rarr;1.93), Tanaka mesh (&xi;<sub>fluc</sub> &approx; 2.7 nm), and droplets</td>
                  </tr>
                  <tr>
                    <td><strong>ND-i-001</strong></td>
                    <td>40 mM, pH 11.3, pure</td>
                    <td class="num">0.0%</td>
                    <td class="num" style="font-weight: 700; color: #34d399;">11.4%</td>
                    <td class="num" style="font-weight: 700; color: #c084fc;">87.1%</td>
                    <td class="num" style="font-weight: 700; color: #fbbf24;">1.5%</td>
                    <td>Instantaneous barrierless quench freezes large proportion of dense viscoelastic PNC mesh before complete vitrification</td>
                  </tr>
                  <tr>
                    <td><strong>ND-i-005</strong></td>
                    <td>40 mM, pH 11.3, +Co²⁺</td>
                    <td class="num">0.0%</td>
                    <td class="num" style="font-weight: 700; color: #34d399;">16.6%</td>
                    <td class="num" style="font-weight: 700; color: #c084fc;">83.4%</td>
                    <td class="num">trace</td>
                    <td>Co²⁺ adsorption lowers interfacial tension &gamma;, stabilizing dense PNC coordination clusters alongside open network</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
        </div>

        <!-- INTERACTIVE TRAJECTORY EXPLORER -->
        <div class="card-title" style="margin-top: 28px; margin-bottom: 12px;">
          <span>Interactive multi-series trajectory explorer</span>
          <span style="font-size: 0.82rem; color: var(--text-muted); font-weight: normal;">Powered by Chart.js</span>
        </div>

        <div class="controls-bar">
          <div class="control-group">
            <span class="control-label">Metric:</span>
            <button class="toggle-btn active" id="btn-metric-alpha" onclick="setChartMetric('alpha')">Conversion &alpha;(t)</button>
            <button class="toggle-btn" id="btn-metric-rmeso" onclick="setChartMetric('r_meso')">Length scale L(t) [nm]</button>
            <button class="toggle-btn" id="btn-metric-phi" onclick="setChartMetric('phi')">Volume fraction &phi;(t)</button>
            <button class="toggle-btn" id="btn-metric-chi2" onclick="setChartMetric('chi2')">Reduced &chi;&sup2;(t)</button>
          </div>

          <div class="control-group">
            <span class="control-label">Time scale:</span>
            <button class="toggle-btn active" id="btn-scale-lin" onclick="setChartScale('linear')">Linear</button>
            <button class="toggle-btn" id="btn-scale-log" onclick="setChartScale('logarithmic')">Logarithmic</button>
          </div>
        </div>

        <div class="controls-bar" style="margin-top: -6px; margin-bottom: 16px;">
          <div class="control-group">
            <span class="control-label">Series toggles:</span>
            <div class="series-checkbox-group">
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-003" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#2563eb;"></span> ND-i-003 (10 mM, pure)</label>
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-004" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#0891b2;"></span> ND-i-004 (15 mM, pure)</label>
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-002" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#16a34a;"></span> ND-i-002 (20 mM, pure)</label>
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-009" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#ea580c;"></span> ND-i-009 (40 mM, primary pure benchmark)</label>
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-001" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#dc2626;"></span> ND-i-001 (40 mM, pure, pH 11.3)</label>
              <label class="series-checkbox-item"><input type="checkbox" id="chk-ND-i-005" checked onchange="updateChartVisibility()"><span class="color-dot" style="background:#9333ea;"></span> ND-i-005 (40 mM, +Co²⁺ additive)</label>
            </div>
          </div>
        </div>

        <div class="chart-wrapper">
          <canvas id="trajectoryChart"></canvas>
        </div>

        <!-- QUANTITATIVE SYNTHESIS TABLE -->
        <div class="card-title" style="margin-top: 28px;">
          <span>Cross-series kinetic and trajectory synthesis table</span>
        </div>

        <div class="table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Series ID</th>
                <th>Conditions</th>
                <th>Phase separation regime</th>
                <th class="num">Induction t<sub>ind</sub> [s]</th>
                <th class="num">Terminal &alpha;<sub>&infin;</sub></th>
                <th>Characteristic meso-scale state</th>
                <th class="num">Droplets &lang;R<sub>lg</sub>&rang; [&mu;m]</th>
                <th class="num">&phi;<sub>lg</sub> &times; 10<sup>5</sup></th>
                <th class="num">Mean &chi;&sup2;<sub>red</sub></th>
                <th>Terminal morphology</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td><strong>ND-i-003</strong></td>
                <td>10 mM, pH 11.3, pure</td>
                <td>Sub-spinodal classical nucleation</td>
                <td class="num">1.41 &plusmn; 0.10</td>
                <td class="num">0.86 (86%)</td>
                <td>&lang;R<sub>acc</sub>&rang;: 60.7 &rarr; 107.0 nm (S(q) &equiv; 1)</td>
                <td class="num">2.62 &plusmn; 0.15</td>
                <td class="num">5.50</td>
                <td class="num" style="color: var(--success);">1.15</td>
                <td>Polydisperse amorphous nanospheres</td>
              </tr>
              <tr>
                <td><strong>ND-i-004</strong></td>
                <td>15 mM, pH 11.3, pure</td>
                <td>Crossover phase boundary</td>
                <td class="num">0.35 &plusmn; 0.05</td>
                <td class="num">0.69 (69%)</td>
                <td>R<sub>prim</sub> = 4.94 nm, &xi; = 108.7 nm + ACC (47.7 nm)</td>
                <td class="num">2.80 &plusmn; 0.20</td>
                <td class="num">4.50</td>
                <td class="num" style="color: var(--success);">1.18</td>
                <td>Hybrid cluster gel + discrete ACC spheres</td>
              </tr>
              <tr>
                <td><strong>ND-i-002</strong></td>
                <td>20 mM, pH 11.3, pure</td>
                <td>Moderate spinodal gelation</td>
                <td class="num">&lt; 0.08</td>
                <td class="num">0.81 (81%)</td>
                <td>R<sub>prim</sub> = 11.2 nm, D<sub>f</sub> = 1.86, &xi; = 356 nm</td>
                <td class="num">2.75 &plusmn; 0.18</td>
                <td class="num">6.00</td>
                <td class="num" style="color: var(--success);">1.25</td>
                <td>Space-filling mass fractal network + Tanaka PNC mesh (24%)</td>
              </tr>
              <tr style="background-color: rgba(234, 88, 12, 0.08); font-weight: 500;">
                <td><strong>ND-i-009</strong></td>
                <td>40 mM, pH 10.0, pure</td>
                <td>Primary pure CaCO₃ benchmark</td>
                <td class="num">0.15 &plusmn; 0.03</td>
                <td class="num">0.81 (81%)</td>
                <td>&xi;(t): 46.7 nm &rarr; &infin;, D<sub>f</sub>(t) &in; [1.88, 1.93]</td>
                <td class="num">2.27 &rarr; 3.20</td>
                <td class="num">120.0</td>
                <td class="num" style="color: var(--success);">1.20</td>
                <td>Infinite network + coarsening droplets + Tanaka PNC mesh (26%)</td>
              </tr>
              <tr>
                <td><strong>ND-i-001</strong></td>
                <td>40 mM, pH 11.3, pure</td>
                <td>Rapid quench spinodal demixing</td>
                <td class="num">&lt; 0.06</td>
                <td class="num">0.92 (92%)</td>
                <td>R<sub>prim</sub> = 3.50 nm, D<sub>f</sub> = 2.33, &xi; = 735 nm</td>
                <td class="num">3.01 &plusmn; 0.12</td>
                <td class="num">4.98</td>
                <td class="num" style="color: var(--success);">1.05</td>
                <td>Bicontinuous struts + vitrified PNC mesh (87%) + droplets</td>
              </tr>
              <tr>
                <td><strong>ND-i-005</strong></td>
                <td>40 mM, pH 11.3, +Co²⁺</td>
                <td>Additive perturbation branch</td>
                <td class="num">&lt; 0.06</td>
                <td class="num">0.88 (88%)</td>
                <td>R<sub>prim</sub> = 6.59 nm, D<sub>f</sub> = 2.10, &xi; = 1,190 nm</td>
                <td class="num">3.53 &plusmn; 0.22</td>
                <td class="num">2.92</td>
                <td class="num" style="color: var(--success);">1.10</td>
                <td>Open fractal network + dopant-stabilized PNC clusters (83%)</td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- SECTION 4: COMPARATIVE PARAMETER MATRIX -->
    <div id="sec-parameters" class="section-content">
      <div class="card">
        <div class="card-title">Calibrated structural parameter matrix across all series</div>
        <p class="text-block">
          Fitted parameters derived from the full-range composite scattering curves (<em>q &in; [0.002, 2.5] nm<sup>-1</sup></em>) 
          combining short-distance (3.0 m) and long-distance (31.0 m) geometries at the terminal reaction stage. 
          Uncertainties represent formal &plusmn;1&sigma; standard errors calculated from the full covariance matrix.
        </p>

        <div class="table-container">
          <table class="data-table">
            <thead>
              <tr>
                <th>Series ID</th>
                <th>[CaCO<sub>3</sub>]</th>
                <th>pH</th>
                <th>Additive</th>
                <th class="num">&lang;R<sub>prim</sub>&rang; [nm]</th>
                <th class="num">&sigma;<sub>prim</sub> [nm]</th>
                <th class="num">D<sub>fractal</sub></th>
                <th class="num">&xi; [nm]</th>
                <th class="num">&lang;R<sub>large</sub>&rang; [&mu;m]</th>
                <th class="num">&phi;<sub>large</sub> &times; 10<sup>5</sup></th>
                <th class="num">&lang;R<sub>acc</sub>&rang; [nm]</th>
                <th class="num">&phi;<sub>acc</sub> &times; 10<sup>4</sup></th>
                <th class="num">&lang;R<sub>pnc</sub>&rang; / &xi;<sub>fluc</sub> [nm]</th>
                <th class="num">&phi;<sub>pnc</sub> &times; 10<sup>4</sup></th>
                <th class="num">&chi;<sup>2</sup><sub>red</sub></th>
              </tr>
            </thead>
            <tbody>
"""

    for meta in SERIES_METADATA:
        sid = meta["id"]
        p = params.get(sid, {})
        chi2 = p.get("redchi") or p.get("chi2_red", 0.0)

        r_prim = p.get("Rmean") or p.get("R_prim", 0.0)
        e_r_prim = p.get("std_err_Rmean", 0.0)
        r_prim_str = f"{r_prim:.2f} &plusmn; {e_r_prim:.2f}" if r_prim > 0 else "&mdash;"

        s_prim = p.get("stdDev", 0.0)
        e_s_prim = p.get("std_err_stdDev", 0.0)
        s_prim_str = f"{s_prim:.2f} &plusmn; {e_s_prim:.2f}" if s_prim > 0 else "&mdash;"

        d_val = p.get("D") or p.get("Dm", 0.0)
        e_d = p.get("std_err_D") or p.get("std_err_Dm", 0.0)
        d_str = f"{d_val:.3f} &plusmn; {e_d:.3f}" if d_val > 0 else "&mdash;"

        xi_val = p.get("xi", 0.0)
        e_xi = p.get("std_err_xi", 0.0)
        if xi_val >= 1e5:
            xi_str = "&infin;"
        elif xi_val > 0:
            xi_str = f"{xi_val:.0f} &plusmn; {e_xi:.0f}"
        else:
            xi_str = "&mdash;"

        r_large = p.get("mean_radius_large", 0.0)
        e_r_large = p.get("std_err_mean_radius_large", 0.0)
        r_large_str = f"{r_large/1000:.2f} &plusmn; {e_r_large/1000:.2f}" if r_large > 0 else "&mdash;"

        phi_large = p.get("total_volume_fraction_large", 0.0)
        e_phi_large = p.get("std_err_total_volume_fraction_large", 0.0)
        phi_l_str = f"{phi_large*1e5:.2f} &plusmn; {e_phi_large*1e5:.2f}" if phi_large > 0 else "&mdash;"

        r_acc = p.get("R_acc", 0.0)
        if sid == "ND-i-003":
            r_acc = p.get("Rmean", 0.0)
        e_r_acc = p.get("std_err_R_acc") or p.get("std_err_Rmean", 0.0)
        r_acc_str = f"{r_acc:.1f} &plusmn; {e_r_acc:.1f}" if r_acc > 0 and r_large == 0 else ("&mdash;" if r_acc == 0 else f"{r_acc:.1f}")

        phi_acc = p.get("phi_acc", 0.0)
        if sid == "ND-i-003":
            phi_acc = p.get("phi", 0.0)
        e_phi_acc = p.get("std_err_phi_acc") or p.get("std_err_phi", 0.0)
        phi_acc_str = f"{phi_acc*1e4:.2f} &plusmn; {e_phi_acc*1e4:.2f}" if phi_acc > 0 else "&mdash;"

        r_pnc = p.get("R_pnc", 0.0)
        e_r_pnc = p.get("std_err_R_pnc", 0.0)
        if sid == "ND-i-009":
            r_pnc_str = "2.72 &ndash; 2.84 (&xi;<sub>fluc</sub>)"
        elif r_pnc > 0:
            r_pnc_str = f"{r_pnc:.2f} &plusmn; {e_r_pnc:.2f}"
        else:
            r_pnc_str = "&mdash;"

        phi_pnc = p.get("phi_pnc", 0.0)
        e_phi_pnc = p.get("std_err_phi_pnc", 0.0)
        if sid == "ND-i-009":
            phi_pnc_str = "25.9% (Q-inv)"
        elif phi_pnc > 0:
            phi_pnc_str = f"{phi_pnc*1e4:.2f} &plusmn; {e_phi_pnc*1e4:.2f}"
        else:
            phi_pnc_str = "&mdash;"

        html += f"""              <tr>
                <td><strong>{sid}</strong></td>
                <td>{meta['conc_mM']:.0f} mM</td>
                <td>{meta['pH']}</td>
                <td>{meta['additive']}</td>
                <td class="num">{r_prim_str}</td>
                <td class="num">{s_prim_str}</td>
                <td class="num">{d_str}</td>
                <td class="num">{xi_str}</td>
                <td class="num">{r_large_str}</td>
                <td class="num">{phi_l_str}</td>
                <td class="num">{r_acc_str}</td>
                <td class="num">{phi_acc_str}</td>
                <td class="num" style="color: #c084fc;">{r_pnc_str}</td>
                <td class="num" style="color: #c084fc;">{phi_pnc_str}</td>
                <td class="num" style="color: var(--success); font-weight: 700;">{chi2:.3f}</td>
              </tr>
"""

    html += """            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- SECTION 5: COMPOSITE SCATTERING CURVES GALLERY -->
    <div id="sec-gallery" class="section-content">
      <div class="card">
        <div class="card-title">Full-range composite scattering curves (q in [0.002, 2.5] nm⁻¹)</div>
        <p class="text-block">
          Experimental composite SAXS scattering curves and converged physical model fits across all 6 reaction series.
          All curves are stitched across the 31.0 m (long-distance) and 3.0 m (short-distance) ESRF ID02 camera configurations.
          Each physical model component follows the standardised, globally unique project palette.
        </p>

        <!-- STANDARDISED PALETTE BANNER -->
        <div style="margin: 0.75rem 0 1.25rem 0; padding: 0.75rem 1rem; background: #070d19; border: 1px solid #1e293b; border-radius: 8px; display: flex; flex-wrap: wrap; gap: 1.25rem; font-size: 0.82rem; align-items: center;">
          <span style="font-weight: 700; color: #cbd5e1; text-transform: uppercase; font-size: 0.75rem; letter-spacing: 0.05em;">Standardised Model Components:</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #94a3b8;"><span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#64748b;"></span> Experimental Data</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #f87171;"><span style="display:inline-block; width:14px; height:3px; background:#dc2626;"></span> Total Unified Fit</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #34d399;"><span style="display:inline-block; width:14px; height:3px; background:#059669;"></span> Spinodal Fractal Network</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #38bdf8;"><span style="display:inline-block; width:14px; height:3px; background:#0284c7;"></span> Discrete ACC Spheres</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #c084fc;"><span style="display:inline-block; width:14px; height:3px; background:#9333ea;"></span> Tanaka PNC Fluctuation Mesh</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #fbbf24;"><span style="display:inline-block; width:14px; height:3px; background:#d97706;"></span> Macroscopic Droplets / Porod Tail</span>
          <span style="display: inline-flex; align-items: center; gap: 0.4rem; color: #94a3b8;"><span style="display:inline-block; width:14px; height:2px; background:#64748b;"></span> Solvent Background</span>
        </div>

        <div class="gallery-grid">
"""

    for meta in SERIES_METADATA:
        sid = meta["id"]
        img_rel = f"./{meta['dir']}/output/composite/composite_baseline_fit.png"
        p = params.get(sid, {})
        chi2 = p.get("redchi") or p.get("chi2_red", 0.0)

        html += f"""          <div class="gallery-item">
            <div class="gallery-caption">
              <span>{sid} ({meta['conc_mM']:.0f} mM, {meta['pH']}, {meta['additive']})</span>
              <span class="badge">&chi;&sup2; = {chi2:.3f}</span>
            </div>
            <img class="gallery-img" src="{img_rel}" alt="Composite scattering curve fit for {sid}">
          </div>
"""

    html += """        </div>
      </div>
    </div>

    <!-- SECTION 6: PHYSICAL SCATTERING FORMALISM -->
    <div id="sec-formalism" class="section-content">
      <div class="card">
        <div class="card-title">Multi-scale physical scattering formalism</div>
        
        <div class="text-block">
          The total coherent scattering intensity \\(I_{\\rm tot}(q)\\) across the momentum transfer range
          \\(q \\in [0.002, 2.5]\\text{ nm}^{-1}\\) is modelled as the linear superposition of decoupled structural domains:
        </div>

        <div class="equation-box">
          \\[
          I_{\\rm tot}(q) = I_{\\rm droplet/ACC}(q) + I_{\\rm network}(q) + I_{\\rm fluc}(q) + I_{\\rm bkg}
          \\]
        </div>

        <div class="text-block">
          <h3>1. Primary mass fractal network</h3>
          <p>
            Formed by polydisperse spherical building blocks following a lognormal size distribution \\(n(R)\\) 
            with Teixeira mass fractal inter-particle structure factor \\(S(q, D_m, \\xi)\\):
          </p>
        </div>

        <div class="equation-box">
          \\[
          I_{\\rm network}(q) = (\\Delta\\rho)^2 \\int_0^\\infty n(R) V^2(R) P(q, R) S(q, R, D_m, \\xi) \\, dR
          \\]
          \\[
          S(q, R, D_m, \\xi) = 1 + \\frac{D_m \\Gamma(D_m - 1)}{(q R)^{D_m}} \\frac{\\sin[(D_m - 1) \\arctan(q \\xi)]}{\\left[1 + (q \\xi)^{-2}\\right]^{(D_m - 1)/2}}
          \\]
        </div>

        <div class="text-block">
          <h3>2. Capillary emulsion droplets and discrete ACC nanospheres</h3>
          <p>
            Spherical particles (either macroscopic dense liquid droplets at low <em>q</em> or discrete ACC nanospheres at intermediate <em>q</em>)
            with lognormal size distribution:
          </p>
        </div>

        <div class="equation-box">
          \\[
          I_{\\rm sphere}(q) = (\\Delta\\rho)^2 \\int_0^\\infty n(R) V^2(R) \\left[ 3 \\frac{\\sin(qR) - qR\\cos(qR)}{(qR)^3} \\right]^2 dR
          \\]
        </div>

        <div class="text-block">
          <h3>3. Dynamic viscoelastic concentration fluctuations and the PNC mesh</h3>
          <p>
            Internal concentration fluctuations of the dense viscoelastic network with correlation length \\(\\xi_{\\rm fluc}\\)
            and Porod surface-decay exponent \\(\\mu = 4.0\\) (Tanaka / Fisher-Burford formulation):
          </p>
        </div>

        <div class="equation-box">
          \\[
          I_{\\rm fluc}(q) = \\frac{A_{\\rm fluc}}{\\left[ 1 + \\frac{1}{6} (q \\xi_{\\rm fluc})^2 \\right]^2}
          \\]
        </div>
        
        <div class="text-block">
          <p>
            Physically, this fluctuation term represents the <strong>transient viscoelastic mesh formed by entangled, hydrated prenucleation clusters (PNCs)</strong> 
            within the dense liquid precursor phase. The correlation length \\(\\xi_{\\rm fluc} \\approx 1.5 - 3.0\\text{ nm}\\) directly corresponds 
            to the mesh spacing between ionic coordination cross-links prior to solid ACC vitrification.
          </p>
        </div>

        <div class="text-block">
          <h3>4. Asymptotic Porod interface</h3>
          <p>
            Scattering from macroscopic phase-separated interfaces satisfying Porod's universal law (\\(I \\propto q^{-4}\\)):
          </p>
        </div>

        <div class="equation-box">
          \\[
          I_{\\rm porod}(q) = \\frac{A_{\\rm porod}}{q^4}
          \\]
        </div>
      </div>
    </div>
  </main>

  <script>
    // Embedded Raw Trajectory Data
    const RAW_TRAJECTORIES = """ + trajectories_json_raw + """;

    // Chart.js State
    let currentMetric = 'alpha';
    let currentScale = 'linear';
    let trajectoryChartInstance = null;

    function switchSection(secId, btn) {
      document.querySelectorAll('.nav-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.section-content').forEach(s => s.classList.remove('active'));
      const targetSec = document.getElementById('sec-' + secId);
      if (targetSec) targetSec.classList.add('active');
      const clickedBtn = btn || (window.event ? (window.event.currentTarget || window.event.target) : null);
      if (clickedBtn && clickedBtn.classList && clickedBtn.classList.contains('nav-btn')) {
        clickedBtn.classList.add('active');
      } else {
        const btnMatch = document.querySelector(`.nav-btn[onclick*="${secId}"]`);
        if (btnMatch) btnMatch.classList.add('active');
      }
      if (secId === 'trajectories' && !trajectoryChartInstance) {
        initTrajectoryChart();
      }
      if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise();
      }
    }

    function setChartMetric(metric) {
      currentMetric = metric;
      document.querySelectorAll('[id^="btn-metric-"]').forEach(b => b.classList.remove('active'));
      const activeBtn = document.getElementById('btn-metric-' + metric);
      if (activeBtn) activeBtn.classList.add('active');
      updateChartData();
    }

    function setChartScale(scale) {
      currentScale = scale;
      document.querySelectorAll('[id^="btn-scale-"]').forEach(b => b.classList.remove('active'));
      const activeBtn = document.getElementById('btn-scale-' + (scale === 'linear' ? 'lin' : 'log'));
      if (activeBtn) activeBtn.classList.add('active');
      if (trajectoryChartInstance) {
        trajectoryChartInstance.options.scales.x.type = (scale === 'logarithmic') ? 'logarithmic' : 'linear';
        if (scale === 'logarithmic') {
          trajectoryChartInstance.options.scales.x.min = 0.04;
        } else {
          delete trajectoryChartInstance.options.scales.x.min;
        }
        trajectoryChartInstance.update();
      }
    }

    function updateChartVisibility() {
      if (!trajectoryChartInstance) return;
      Object.keys(RAW_TRAJECTORIES).forEach((sid, idx) => {
        const chk = document.getElementById('chk-' + sid);
        if (chk && trajectoryChartInstance.data.datasets[idx]) {
          trajectoryChartInstance.data.datasets[idx].hidden = !chk.checked;
        }
      });
      trajectoryChartInstance.update();
    }

    function getMetricLabel(metric) {
      switch(metric) {
        case 'alpha': return 'Conversion extent α(t) = φ(t) / φ_max';
        case 'r_meso': return 'Characteristic length scale L(t) [nm] (R_ACC or ξ)';
        case 'phi': return 'Dense phase volume fraction φ(t)';
        case 'chi2': return 'Reduced chi-squared χ²_red(t)';
        default: return metric;
      }
    }

    function initTrajectoryChart() {
      const ctx = document.getElementById('trajectoryChart').getContext('2d');
      const datasets = [];

      Object.keys(RAW_TRAJECTORIES).forEach(sid => {
        const d = RAW_TRAJECTORIES[sid];
        const pts = [];
        for (let i = 0; i < d.time.length; i++) {
          const val = d[currentMetric][i];
          if (val !== null && isFinite(val)) {
            pts.push({ x: d.time[i], y: val });
          }
        }
        datasets.push({
          label: d.name,
          data: pts,
          borderColor: d.color,
          backgroundColor: d.color,
          borderWidth: 2.2,
          pointRadius: 2.5,
          pointHoverRadius: 5.5,
          tension: 0.15,
          showLine: true
        });
      });

      trajectoryChartInstance = new Chart(ctx, {
        type: 'scatter',
        data: { datasets: datasets },
        options: {
          responsive: true,
          maintainAspectRatio: false,
          animation: { duration: 300 },
          plugins: {
            legend: {
              display: true,
              position: 'top',
              labels: {
                color: '#cbd5e1',
                font: { family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto', size: 11 },
                boxWidth: 14
              }
            },
            tooltip: {
              backgroundColor: '#16223b',
              titleColor: '#ffffff',
              bodyColor: '#cbd5e1',
              borderColor: '#233354',
              borderWidth: 1,
              callbacks: {
                label: function(context) {
                  return context.dataset.label + ': t = ' + context.parsed.x.toFixed(3) + ' s, val = ' + context.parsed.y;
                }
              }
            }
          },
          scales: {
            x: {
              type: 'linear',
              position: 'bottom',
              title: {
                display: true,
                text: 'Reaction time t [s]',
                color: '#94a3b8',
                font: { weight: 'bold', size: 11 }
              },
              grid: { color: 'rgba(255, 255, 255, 0.06)' },
              ticks: { color: '#94a3b8' }
            },
            y: {
              title: {
                display: true,
                text: getMetricLabel(currentMetric),
                color: '#94a3b8',
                font: { weight: 'bold', size: 11 }
              },
              grid: { color: 'rgba(255, 255, 255, 0.06)' },
              ticks: { color: '#94a3b8' }
            }
          }
        }
      });
      updateChartVisibility();
    }

    function updateChartData() {
      if (!trajectoryChartInstance) return;
      Object.keys(RAW_TRAJECTORIES).forEach((sid, idx) => {
        const d = RAW_TRAJECTORIES[sid];
        const pts = [];
        for (let i = 0; i < d.time.length; i++) {
          const val = d[currentMetric][i];
          if (val !== null && isFinite(val)) {
            pts.push({ x: d.time[i], y: val });
          }
        }
        trajectoryChartInstance.data.datasets[idx].data = pts;
      });
      trajectoryChartInstance.options.scales.y.title.text = getMetricLabel(currentMetric);
      trajectoryChartInstance.update();
    }

    window.addEventListener('DOMContentLoaded', () => {
      const hash = window.location.hash.replace('#sec-', '').replace('#', '');
      if (hash && document.getElementById('sec-' + hash)) {
        switchSection(hash);
      }
    });
  </script>
</body>
</html>
"""

    out_file = os.path.join(RESULTS_DIR, "index.html")
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Master synthesis dashboard generated successfully at: {out_file}")

    dash_dir = os.path.join(RESULTS_DIR, "dashboard")
    os.makedirs(dash_dir, exist_ok=True)
    dest_path = os.path.join(dash_dir, "index.html")
    if os.path.exists(dest_path) or os.path.islink(dest_path):
        os.remove(dest_path)

    # Adjust relative paths for dashboard/index.html since it resides one level deeper than RESULTS_DIR
    dash_html = html.replace('href="./ND-i-', 'href="../ND-i-')
    dash_html = dash_html.replace('src="./output/', 'src="../output/')
    dash_html = dash_html.replace('src="./ND-i-', 'src="../ND-i-')

    with open(dest_path, "w", encoding="utf-8") as f:
        f.write(dash_html)
    print(f"Master synthesis dashboard generated with adapted relative paths at: {dest_path}")


if __name__ == "__main__":
    build_master_dashboard()
