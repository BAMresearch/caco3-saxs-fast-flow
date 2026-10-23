# Unified SAXS analysis suite: dataset ND-i-006

> [!WARNING]
> **EXCLUDED DATASET : DO NOT USE**
> This dataset has been excluded from analysis and comparison suites. Any preliminary fits have been removed.

## 1. Experimental conditions
- **System**: Pure CaCO3 (40.0 mM, pH 9.75)
- **Beamline**: ESRF ID02 High-Brilliance SAXS/WAXS Beamline
- **Detectors / Geometries**:
  - Short Distance: 3.0 m camera ($q \in [0.022, 2.54]\text{ nm}^{-1}$), 50 frames ($t = 0.060\text{ s}$ to $10.877\text{ s}$)
  - Long Distance: 31.0 m camera ($q \in [0.002, 0.22]\text{ nm}^{-1}$), 50 frames ($t = 0.060\text{ s}$ to $10.877\text{ s}$)
- **Regime**: Spinodal Mass Fractal Network & Secondary Large Droplets / Spheres.
  - **Morphology Confirmation**: Contains **only mass fractal networks and secondary large spheres/droplets**, but **no discrete ACC particles** (no intermediate 30–50 nm particle form-factor knee).

## 2. Directory structure
```
results/ND-i-006_40mM_pH9.75/
├── data/
│   ├── dataNDi006_fitting_results_short_distance.csv
│   ├── dataNDi006_fitting_results_long_distance.csv
│   ├── dataNDi006_fitting_results_composite.csv
│   ├── composite_final.dat
│   └── frame_metadata.csv
├── models/
│   ├── __init__.py
│   ├── fitting_engine.py
│   └── scattering_models.py
├── output/
│   ├── composite/
│   ├── short_distance/
│   ├── long_distance/
│   ├── unified_trajectories/
│   └── dashboard/
│       └── index.html
└── scripts/
    ├── 00_prepare_data_and_composite.py
    ├── 00_fit_composite_baseline.py
    ├── 01_fit_short_distance_models.py
    ├── 02_fit_long_distance_models.py
    ├── 03_generate_all_frame_plots.py
    ├── 04_compare_parameter_trajectories.py
    ├── 05_build_unified_dashboard.py
    └── run_unified_pipeline.py
```

## 3. Execution
Run the full unified pipeline end-to-end:
```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --all
```
