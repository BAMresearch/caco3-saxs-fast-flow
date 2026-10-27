#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPS Verification Report Generator
---------------------------------
Data Processing Responsible: Tomasz Stawski, tomasz.stawski@bam.de

Compiles the Viscoelastic Phase Separation (VPS) Framework Verification Report
by dynamically fetching standalone Markdown narratives from `results/docs/narratives/`
and inserting calculated parameter synthesis tables from the analysis engine.

Pure reporting module:
- No hardcoded narrative prose or scientific text in Python code.
- Fetches standalone Markdown sections dynamically.
- Injects dynamic parameter synthesis tables and statistical metrics.

Output:
  results/output/vps_framework_verification_report.md
"""

import os
import sys

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
NARRATIVES_DIR = os.path.join(RESULTS_DIR, 'docs', 'narratives')
OUT_DIR = os.path.join(RESULTS_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(RESULTS_DIR, 'scripts', 'analysis'))
from compute_vps_metrics import load_all_parameters, compute_parameter_summary_table


def load_markdown_section(filename, strip_h1=True):
    """Load content from a standalone markdown file."""
    filepath = os.path.join(NARRATIVES_DIR, filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Narrative file not found: {filepath}")
    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read().strip()
    if strip_h1 and text.startswith('# '):
        # Remove the top-level H1 header and any immediate blank lines
        lines = text.split('\n')
        text = '\n'.join(lines[1:]).lstrip('\n')
    return text


def build_parameter_table_markdown(summary_rows):
    """Format parameter synthesis rows into GitHub-flavored Markdown table."""
    lines = [
        "| Series ID | Concentration | pH | Physical regime | Primary radius $R_{\\text{prim}}$ [nm] | Fractal dimension $D_m$ | Network cutoff $\\xi$ [nm] | Droplet radius $\\langle R_{\\text{large}} \\rangle$ [$\\mu\\text{m}$] | Reduced $\\chi^2$ |",
        "| :--- | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: |"
    ]
    for r in summary_rows:
        row_str = (
            f"| **{r['series_id']}** | {r['concentration_mM']} mM | {r['pH']} | "
            f"{r['regime']} | {r['r_prim']} | {r['dm']} | {r['xi']} | "
            f"{r['r_large_um']} | {r['chi2_red']} |"
        )
        lines.append(row_str)
    return "\n".join(lines)


def generate_report():
    """Dynamically assemble the verification report from standalone MD and data."""
    # 1. Load standalone narratives
    exec_summary_md = load_markdown_section('executive_summary_vps.md')
    vps_theory_md = load_markdown_section('vps_theory.md')
    regimes_md = load_markdown_section('regime_descriptions.md')

    # 2. Compute dynamic parameter table
    data = load_all_parameters()
    summary_rows = compute_parameter_summary_table(data)
    param_table_md = build_parameter_table_markdown(summary_rows)

    # 3. Assemble document sections
    report = []
    report.append("# Viscoelastic phase separation (VPS) framework verification report")
    report.append("")
    report.append("**Data processing responsible:** Tomasz M. Stawski (tomasz.stawski@bam.de)")
    report.append("**Theoretical framework:** Viscoelastic phase separation (VPS)")
    report.append("")
    report.append("---")
    report.append("")
    report.append("## 1. Executive summary")
    report.append("")
    report.append(exec_summary_md)
    report.append("")
    report.append("---")
    report.append("")
    report.append("## 2. Parameter synthesis across regimes")
    report.append("")
    report.append(param_table_md)
    report.append("")
    report.append("---")
    report.append("")
    report.append("## 3. Mathematical verification of Tanaka VPS concepts")
    report.append("")
    report.append(vps_theory_md)
    report.append("")
    report.append("---")
    report.append("")
    report.append("## 4. Logical transition continuum across experimental regimes")
    report.append("")
    report.append(regimes_md)
    report.append("")

    full_report_text = "\n".join(report)
    out_file = os.path.join(OUT_DIR, 'vps_framework_verification_report.md')
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(full_report_text)
    print(f"[Reporting] Saved VPS verification report to: {out_file}")
    return out_file


def main():
    print("=" * 80)
    print("  RUNNING: VPS Verification Report Generator")
    print("=" * 80)
    out = generate_report()
    print(f"[Success] Verification report generated at: {out}")


if __name__ == '__main__':
    main()
