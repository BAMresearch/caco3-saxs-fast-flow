# Unified multi-scale SAXS fitting suite (dataset ND-i-003)

This repository contains the clean, unified SAXS kinetic fitting suite for dataset **ND-i-003** ($10.0\text{ mM }\text{CaCO}_3$, pH 11.3, pure uninhibited system).

It integrates the validated theoretical concepts across short distance ($3.0\text{ m}$ camera, 50 frames), long distance ($31.0\text{ m}$ camera, 50 frames), and full-range composite scattering curves.

---

## Physical context and scientific models

### 1. Low supersaturation regime and absence of spinodal gelation
In contrast to high-supersaturation regimes ($40\text{ mM}$ in ND-i-001 and ND-i-009) where barrierless phase separation drives the instantaneous formation of a bicontinuous mass fractal gel network ($D_m \approx 1.84$), the $10\text{ mM}$ scenario in **ND-i-003** resides well below the spinodal boundary:
* **Suppression of Mass Fractal Gelation**: The Teixeira mass-fractal structure factor $S(q)$ and pre-nucleation cluster percolation are completely inapplicable.
* **Classical Homogeneous Nucleation**: Phase transformation proceeds via the classical nucleation and growth of discrete, non-percolating amorphous calcium carbonate (ACC) nanospheres.
* **Textbook Porod Scattering**: The scattering at intermediate and high $q$ displays asymptotic $q^{-4}$ power-law decay ($D=3$), confirming compact, dense 3D spherical particles with sharp interfaces rather than diffuse mass fractals ($q^{-D_m}$ with $D_m < 3$).

### 2. Short distance ($q \in [0.022, 2.54]\text{ nm}^{-1}$, 50 frames, $t = 0.060 \to 10.877\text{ s}$)
* **Pure Flat Solvent Background**:
  - The short camera distance ($3.0\text{ m}$) detects no primary building blocks, high-$q$ excess, or prenucleation cluster fluctuations.
  - The signal is rigorously modeled by a flat solvent background:
    $$I(q) = I_{\rm bkg}$$
  - Fitted solvent background hovers at $I_{\rm bkg} \approx (0.7 - 1.0) \times 10^{-4}\text{ cm}^{-1}$.
  - The mean subtracted intensity is centered symmetrically around zero ($\pm 0.0002\text{ cm}^{-1}$), with $\approx 55 - 60\%$ positive data points (pure Gaussian noise).
  - **Quality of Fit**: Reduced $\chi^2_{\rm red} \approx 0.09$, verifying that no sub-$10\text{ nm}$ primary nodes exist.

### 3. Long distance ($q \in [0.002, 0.25]\text{ nm}^{-1}$, 50 frames, $t = 0.060 \to 10.877\text{ s}$)
* **Unified Multi-Component Scattering Model**:
  - Unexplained low-$q$ intensities ($q \in [0.002, 0.01]\text{ nm}^{-1}$) exhibit a clean $q^{-4}$ power-law decay ($I \sim 10^2 - 10^3\text{ cm}^{-1}$ at $q = 0.0025\text{ nm}^{-1}$), arising from asymptotic surface scattering off macroscopic interfaces ($R \gg 500\text{ nm}$, e.g. capillary wall wetting, meniscus, coarse precipitates).
  - The total scattering is rigorously modeled across all frames by:
    $$I_{\rm tot}(q) = I_{\rm sphere}(q) + I_{\rm porod}(q) + I_{\rm bkg}$$
    where:
    $$I_{\rm sphere}(q) = (\Delta\rho)^2 \int_0^\infty n(r) V(r)^2 P(q, r)\, dr$$
    $$I_{\rm porod}(q) = \frac{A_{\rm porod}}{q^4}$$
  - **Induction Period ($t < 2.5\text{ s}$, Frames 0–32)**:
    - Pure macroscopic Porod interface scattering ($I(q) = A_{\rm porod}/q^4 + I_{\rm bkg}$) with $A_{\rm porod} \approx (0.5 - 1.0) \times 10^{-8}\text{ cm}^{-1}\text{nm}^{-4}$.
    - Reduced $\chi^2_{\rm red} \in [0.07, 0.28]$, capturing the low-$q$ upturn while $\phi_{\rm ACC} = 0$.
  - **Nucleation & Growth ($t = 2.45 \to 10.877\text{ s}$, Frames 33–49)**:
    - Fitted via backward evolutionary propagation ($49 \to 33$), guaranteeing smooth, monotonic convergence:
    - Mean radius grows smoothly from $\langle R \rangle \approx 64.4\text{ nm}$ to $101.8\text{ nm}$ (diameter $129 \to 204\text{ nm}$).
    - Size dispersion evolves from $\sigma_R \approx 16.1\text{ nm}$ to $23.5\text{ nm}$ ($\approx 23\%$ relative polydispersity).
    - ACC volume fraction $\phi(t)$ increases monotonically from $0.6 \times 10^{-4}$ to $(4.41 \pm 0.03) \times 10^{-4}$.
    - Porod amplitude scales smoothly to $A_{\rm porod} \approx (2.76 \pm 0.06) \times 10^{-8}\text{ cm}^{-1}\text{nm}^{-4}$.
  - **Quantitative Stoichiometric Agreement**:
    For $10.0\text{ mM }\text{CaCO}_3$ ($0.010\text{ mol/L}$), assuming $\rho_{\rm ACC} \approx 2.0\text{ g/cm}^3$ and $M = 100.09\text{ g/mol}$, the maximum theoretical volume fraction is:
    $$\phi_{\rm max} = \frac{0.010 \times 100.09}{2.0 \times 1000} = 5.0 \times 10^{-4}$$
    The fitted plateau $\phi = 4.41 \times 10^{-4}$ represents $\approx 88.2\%$ of the total stoichiometric precursor budget (and composite baseline yields $\phi = 4.42 \times 10^{-4}$), confirming quantitative conversion of precursor ions into discrete ACC nanospheres.
  - **Quality of Fit**: Reduced $\chi^2_{\rm red} \sim 0.08 - 0.19$ across the entire growth trajectory (substantially reduced from $\approx 0.75$ in earlier sphere-only fits).

### 4. Rigorous parameter uncertainties and visualization
* All fitted parameters include formal standard errors ($\pm 1\sigma$ confidence limits) propagated from the covariance matrix of `lmfit`.
* All trajectory figures feature shaded $\pm 1\sigma$ uncertainty bands.
* Diagnostic frame plots show full component decompositions ($I_{\rm total}$, $I_{\rm sphere}$, $I_{\rm porod}$, $I_{\rm bkg}$) with legends placed at `loc='upper right'` to prevent obstructing scattering curves.
* The interactive dashboard displays live $\pm 1\sigma$ error margins on all parameters across both camera distances.

### 5. Companion non-parametric validation: McSAS3 suite
* In addition to the parametric lognormal models here, an **independent non-parametric Monte Carlo size distribution analysis** has been executed in the companion directory [`mcsas_fitting_ND-i-003/`](../mcsas_fitting_ND-i-003/).
* Both analyses coexist in parallel: the parametric suite provides rigorous analytical volume fractions ($\phi$), baseline solvent metrics, and composite scattering curves, while the McSAS3 suite provides shape-free Monte Carlo validation of $\langle R \rangle(t)$ and $D_V(R)$.

---

## Directory organisation

```
unified_fitting_ND-i-003/
├── README.md                                         # Project documentation & physical framework
├── data/
│   ├── dataNDi003_fitting_results_short_distance.csv # 3.0 m camera dataset (50 frames, 1638 q-pts)
│   ├── dataNDi003_fitting_results_long_distance.csv  # 31.0 m camera dataset (50 frames, 1645 q-pts)
│   ├── dataNDi003_fitting_results_composite.csv      # Full-range composite curve (1940 q-pts)
│   ├── composite_final.dat                           # Stitched 3-column SAXS curve (q, I, dI)
│   └── frame_metadata.csv                            # Real timestamps (0.060 to 10.877 s) & source files
├── models/
│   ├── __init__.py
│   ├── scattering_models.py                          # Vectorized sphere form factor & lognormal distribution
│   └── fitting_engine.py                             # Sequential parameter regularizers & lmfit builders
├── scripts/
│   ├── 00_prepare_data_and_composite.py              # Extraction of raw ID02 scans & composite stitching
│   ├── 00_fit_composite_baseline.py                  # Full-range composite baseline fit
│   ├── 01_fit_short_distance_models.py               # Short-distance flat solvent background fits
│   ├── 02_fit_long_distance_models.py                # Long-distance discrete ACC nanosphere fits
│   ├── 03_generate_all_frame_plots.py                # Publication-grade diagnostic plots (100 PNGs)
│   ├── 04_compare_parameter_trajectories.py          # Unified multi-scale trajectory synthesis
│   ├── 05_build_unified_dashboard.py                 # Master interactive HTML reviewer
│   └── run_unified_pipeline.py                       # End-to-end pipeline runner
└── output/
    ├── composite/                                    # Full-range composite baseline parameters & plots
    ├── merged/                                       # Mirror alias of composite baseline
    ├── short_distance/
    │   ├── parameters.csv                            # Fitted background, chi2, and noise statistics
    │   ├── fitted_curves.csv                         # Model evaluation across all q
    │   ├── trajectories.png                          # Background & noise metric trajectories (+/- 1 sigma)
    │   └── fits/                                     # Diagnostic plots (fit_frame_000.png to 049.png)
    ├── long_distance/
    │   ├── parameters.csv                            # Fitted <R>, sigma_R, phi, chi2, and induction tags
    │   ├── fitted_curves.csv                         # Model evaluation across all q
    │   ├── trajectories.png                          # Kinetic growth trajectories (+/- 1 sigma)
    │   └── fits/                                     # Diagnostic plots (fit_frame_000.png to 049.png)
    ├── unified_trajectories/
    │   └── unified_cross_scale_trajectories.png      # Length scales, stoichiometry, and cross-scale divergence
    └── dashboard/
        └── index.html                                # Standalone interactive review application
```

---

## Pipeline execution

To run the complete workflow from the project directory:

```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --all
```

To run individual steps:
```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 0   # Composite baseline fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 1   # Short-distance background fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 2   # Long-distance ACC nanosphere fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 3   # Generate diagnostic frame plots
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 4   # Multi-scale parameter trajectories
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 5   # Build interactive HTML dashboard
```

To open the master interactive reviewer:
```bash
xdg-open output/dashboard/index.html
```
