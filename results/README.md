# Multi-Scale SAXS Kinetic Analysis Suite: Architecture & Workflow Guide

Data Processing Responsible: Tomasz M. Stawski (`tomasz.stawski@bam.de`)  
Federal Institute for Materials Research and Testing (BAM)

---

## 1. Architectural Principles

This library provides a modular, decoupled framework for synchrotron time-resolved small-angle X-ray scattering (SAXS) data processing, kinetic fitting, visualization, reporting, and dashboard generation.

### Key Architectural Separations:
1. **Decoupling Narrative Content from Code**:
   - Scientific prose, theoretical definitions, and regime interpretations live as **standalone Markdown files** in [results/docs/narratives/](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/docs/narratives/).
   - Scripts do **not** hardcode lengthy commentary strings. Instead, reporting and dashboard tools fetch the standalone Markdown documents dynamically and inject computed tables/metrics.
   - You can edit theoretical descriptions and physical interpretations in pure Markdown without touching Python code.

2. **Functional Decomposition of Scripts**:
   - **Analysis / Calculation (`scripts/analysis/`)**: Handles data stitching, composite fitting, time-series fitting, parameter covariance estimation, and cross-series metric extraction. Free of Matplotlib rendering and HTML templating.
   - **Plotting / Visualization (`scripts/plotting/`)**: Focuses strictly on figures, diagnostic multi-frame fit plots, publication 4-panel trajectory syntheses, and SVG diagrams.
   - **Reporting (`scripts/reporting/`)**: Merges dynamic calculation outputs with standalone Markdown narratives to generate structured reports (e.g. [vps_framework_verification_report.md](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/output/vps_framework_verification_report.md)).
   - **Dashboards (`scripts/dashboard/`)**: Compiles interactive web dashboards (e.g. [results/index.html](file:///home/tomek/Insync/tomasz.stawski@gmail.com/Google%20Drive/WORK/Projects%20current/project%20-%20CaCO3/results/index.html) and per-dataset dashboards).

---

## 2. Directory Structure

```
results/
├── docs/
│   └── narratives/                          # Standalone scientific narratives (Markdown)
│       ├── executive_summary_vps.md         # Executive findings
│       ├── vps_theory.md                    # Tanaka VPS theory & mathematical models
│       ├── regime_descriptions.md           # 5 physical regimes across 10-40 mM
│       └── scattering_models.md             # Model equations (fractal, form factor, etc.)
│
├── scripts/                                 # Collective scripts (cross-series)
│   ├── run_collective_pipeline.py           # Unified collective pipeline orchestrator
│   ├── analysis/                            # Cross-series computations
│   │   ├── compute_vps_metrics.py
│   │   └── extract_cross_series_trajectories.py
│   ├── plotting/                            # Collective figures & diagrams
│   │   ├── plot_vps_framework.py
│   │   ├── plot_cross_series_trajectories.py
│   │   ├── plot_standardized_composite_baselines.py
│   │   └── generate_model_transformation_diagram.py
│   ├── reporting/                           # Dynamic report assembly
│   │   └── generate_vps_verification_report.py
│   └── dashboard/                           # Master dashboard builder
│       └── build_master_dashboard.py
│
├── ND-i-001_40mM_pH11.3/                    # Dataset folders (ND-i-001 to ND-i-009)
│   ├── data/                                # 1D SAXS curves (raw & merged)
│   ├── models/                              # Mathematical form & structure factors
│   ├── output/                              # Fits, parameters, plots, dashboards
│   └── scripts/
│       ├── run_unified_pipeline.py          # Dataset pipeline orchestrator
│       ├── analysis/                        # Fitting & data prep scripts
│       ├── plotting/                        # Frame & trajectory plot generators
│       ├── dashboard/                       # Dataset interactive dashboard builder
│       └── [00_.. to 05_.. wrappers]       # Backward-compatible shims
```

---

## 3. Workflow Execution

### A. Collective Workflows (Summaries & Cross-Series Comparisons)
Using the collective orchestrator `results/scripts/run_collective_pipeline.py`:

```bash
# 1. Run only analysis (metrics calculation and trajectory extraction):
python results/scripts/run_collective_pipeline.py --fit

# 2. Run only plotting (trajectory comparisons, VPS 4-panel, SVG diagrams):
python results/scripts/run_collective_pipeline.py --plot

# 3. Re-compile the VPS verification report from standalone narratives:
python results/scripts/run_collective_pipeline.py --report

# 4. Re-compile the master interactive HTML dashboard:
python results/scripts/run_collective_pipeline.py --dashboard

# 5. Run everything end-to-end:
python results/scripts/run_collective_pipeline.py --all
```

### B. Dataset-Specific Workflows (Individual Series)
In any dataset directory (e.g., `results/ND-i-004_15mM_pH11.3/scripts/`):

```bash
# 1. Tweak and re-run fitting models without generating plots or dashboards:
python scripts/run_unified_pipeline.py --fit

# 2. Re-render only diagnostic frame plots and parameter trajectories:
python scripts/run_unified_pipeline.py --plot

# 3. Re-build only the dataset interactive HTML dashboard:
python scripts/run_unified_pipeline.py --dashboard

# 4. Run an individual step (backward-compatible):
python scripts/run_unified_pipeline.py --step 2

# 5. Run end-to-end:
python scripts/run_unified_pipeline.py --all
```

---

## 4. Editing Scientific Interpretations

To update theoretical commentary, mathematical derivations, or physical conclusions:
1. Open the relevant file in `results/docs/narratives/` (e.g. `vps_theory.md` or `regime_descriptions.md`).
2. Edit directly in standard GitHub-Flavored Markdown.
3. Run `python results/scripts/run_collective_pipeline.py --report` to re-assemble `vps_framework_verification_report.md` with updated text and freshly calculated parameters.
