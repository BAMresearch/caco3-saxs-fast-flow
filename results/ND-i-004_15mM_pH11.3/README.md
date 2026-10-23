# Unified multi-scale SAXS fitting suite (dataset ND-i-004)

This repository contains the complete, unified SAXS kinetic fitting suite for dataset **ND-i-004** ($15.0\text{ mM }\text{CaCO}_3$, pH 11.3, pure system without additive, ESRF ID02 synchrotron time-resolved SAXS).

It implements the validated **Hybrid Multi-Scale Physical Model** combining:
1. **Mesoscale Spinodal Mass-Fractal Network with Primary Building Blocks** ($R_{\rm prim} = 3.5\text{ nm}$, $\sigma_{\rm prim} = 1.5\text{ nm}$, $D_m \approx 2.33 \to 2.48$, $\xi \approx 20 \to 197\text{ nm}$).
2. **Primary Building Block Form Factor Across Both Scales**: At long distance ($q \in [0.002, 0.22]\text{ nm}^{-1}$), inter-particle correlations produce the mass-fractal power law $q^{-D_m}$. At short distance ($q \in [0.022, 2.54]\text{ nm}^{-1}$), $S(q) \to 1$, naturally reducing the network scattering to the physical form factor $P(q, R_{\rm prim})$ of the primary building blocks.
3. **Discrete Amorphous Calcium Carbonate (ACC) Nanosphere Nucleation & Growth** ($\langle R_{\rm ACC} \rangle \approx 15.9 \to 43.8\text{ nm}$, narrow polydispersity $\sigma \approx 7.4 \to 8.1\text{ nm}$, $\sim 17\%$), giving clear form-factor Bessel oscillations.
4. **Giant Phase-Separated Coarsening Domains** represented via analytical asymptotic Porod surface scattering ($A_{\rm porod} / q^4$).

---

## Scientific rationale and the hybrid model at 15 mM

At $15.0\text{ mM}$ supersaturation and $\text{pH }11.3$, precipitation lies directly at the physical crossover between:
* **Spinodal Liquid-Liquid Demixing / Gelation** (characteristic of the high-concentration $40\text{ mM}$ series ND-i-001 and ND-i-009, where low-$q$ mass fractal scattering $S(q, D_m, \xi)$ dominates).
* **Classical Nucleation & Growth of Polydisperse ACC Droplets** (characteristic of the $19.2\text{ mM}$ series ND-i-003, where high-$q$ $P(q, R, \sigma)$ form-factor oscillations and steep $-4$ Porod slopes dominate).

### Quantitative proof and backward evolutionary strategy:
1. **Unified Building Block Architecture**: The Spinodal Mass Fractal Network does not exist without a physical building block form factor. By anchoring $R_{\rm prim} = 3.5\text{ nm}$ ($\sigma = 1.5\text{ nm}$) across both cameras, the network scattering at long distance and the high-$q$ plateau at short distance are unified under a single physical framework.
2. **Composite / Merged Baseline Fit**: Refitted across all four decades of $q \in [0.002, 2.5]\text{ nm}^{-1}$ using the unified hybrid model (Primary Building Blocks + Teixeira Structure Factor + ACC Nanospheres + Porod Tail + Background), resolving the full multi-scale hierarchy.
3. **Backward Evolutionary Propagation ($49 \to 0$)**: Forward fitting from $t=0$ suffers from low SNR and nascent droplet ambiguity, causing optimizers to inflate polydispersity ($\sigma \to 12.5\text{ nm}$, smearing out oscillations). Backward propagation seeds from the mature state at Frame 49 ($t = 10.88\text{ s}$), where SNR is highest and Bessel oscillations ($q \in [0.08, 0.35]\text{ nm}^{-1}$) are sharp. Propagating backward smoothly tracks droplet contraction to early nuclei ($R \approx 15.9\text{ nm}$) and network contraction ($\xi: 197 \to 20\text{ nm}$) without boundary artifacts.
4. **Short Distance ($3.0\text{ m}$, $q \in [0.022, 2.54]\text{ nm}^{-1}$)**: Combining discrete polydisperse ACC nanospheres ($R = 15.9 \to 43.8\text{ nm}$, $\sigma = 7.4 \to 8.1\text{ nm}$) with asymptotic Porod scattering ($A_{\rm porod} / q^4$) and primary building blocks yields exceptionally low residuals ($\chi^2_{\rm red} = 0.095 \to 0.360$) and captures all short-distance oscillations. All diagnostic plot legends are positioned at `loc='upper right'` to prevent obscuring data points.
5. **Long Distance ($31.0\text{ m}$, $q \in [0.002, 0.22]\text{ nm}^{-1}$)**: The low-$q$ power law matches the Teixeira mass-fractal structure factor $S(q, D_m, \xi)$, while the high-$q$ end seamlessly matches the transferred ACC nanosphere population. Fitting this combined model backward yields $\chi^2_{\rm red} = 0.064 \to 0.212$ across all 50 frames.

---

## Model mathematical formulation

The total observed scattering intensity is given by:

$$I_{\rm tot}(q, t) = I_{\rm net}(q, t) + I_{\rm porod}(q, t) + I_{\rm ACC}(q, t) + I_{\rm bkg}(t)$$

### 1. Meso-scale percolated network (Teixeira fractal with building blocks)
$$I_{\rm net}(q) = \phi_{\rm net} (\Delta \rho)^2 \int_0^\infty n(R) V(R)^2 P(q, R) S(q, R, D_m, \xi) \, dR$$
$$S(q, R, D_m, \xi) = 1 + \frac{D_m \Gamma(D_m - 1)}{(q R)^{D_m}} \left[ 1 + \frac{1}{(q \xi)^2} \right]^{-\frac{D_m - 1}{2}} \sin\left[ (D_m - 1) \arctan(q \xi) \right]$$
*For short distance ($q \gg \xi^{-1}$ and $q R_{\rm prim} \approx 1$), $S(q) \to 1$, reducing $I_{\rm net}(q)$ to the building block form factor.*

### 2. Discrete ACC nanospheres (lognormal distribution)
$$I_{\rm ACC}(q) = (\Delta \rho)^2 \int_0^\infty \left(\frac{4}{3}\pi R^3\right)^2 \left[ 3 \frac{\sin(qR) - qR\cos(qR)}{(qR)^3} \right]^2 N(R)\, dR$$
$$N(R) = \frac{\phi_{\rm ACC}}{\frac{4}{3}\pi R_0^3 \exp(9\sigma_{\rm log}^2 / 2)} \frac{1}{\sqrt{2\pi}\sigma_{\rm log} R} \exp\left( - \frac{(\ln R - \ln R_0)^2}{2\sigma_{\rm log}^2} \right)$$

### 3. Asymptotic Porod surface scattering (giant coarsening domains, $q R \gg 1$)
$$I_{\rm porod}(q) = \frac{A_{\rm porod}}{q^4}$$
Eliminates Bessel function oscillation aliasing and accelerates frame convergence by $> 20\times$.

---

## Key kinetic results summary (ND-i-004)

| Reaction time $t$ | $\langle R_{\rm ACC} \rangle$ [nm] | $\sigma_{\rm ACC}$ [nm] | $\phi_{\rm ACC} \times 10^5$ | $\phi_{\rm net} \times 10^5$ | $D_m$ (network) | $\xi$ (network) [nm] | $A_{\rm porod} \times 10^6$ (SD) | $\chi^2_{\rm red}$ (SD) | $\chi^2_{\rm red}$ (LD) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.060 s** | $15.89 \pm 0.35$ | $7.45 \pm 0.32$ | $0.37 \pm 0.04$ | $0.01 \pm 0.00$ | $2.48 \pm 0.04$ | $20.2 \pm 3.1$ | $0.00 \pm 0.01$ | **0.095** | **0.064** |
| **0.269 s** | $36.61 \pm 0.05$ | $8.07 \pm 0.08$ | $2.47 \pm 0.03$ | $0.07 \pm 0.01$ | $2.47 \pm 0.03$ | $55.4 \pm 3.9$ | $1.65 \pm 0.05$ | **0.149** | **0.125** |
| **0.733 s** | $41.75 \pm 0.05$ | $7.90 \pm 0.07$ | $3.18 \pm 0.04$ | $0.14 \pm 0.01$ | $2.47 \pm 0.02$ | $70.1 \pm 4.2$ | $3.51 \pm 0.06$ | **0.268** | **0.152** |
| **1.855 s** | $42.62 \pm 0.06$ | $7.98 \pm 0.08$ | $2.95 \pm 0.04$ | $0.18 \pm 0.01$ | $2.42 \pm 0.02$ | $112.1 \pm 5.8$ | $5.76 \pm 0.07$ | **0.329** | **0.170** |
| **4.686 s** | $43.33 \pm 0.06$ | $7.62 \pm 0.07$ | $2.83 \pm 0.03$ | $0.24 \pm 0.02$ | $2.37 \pm 0.02$ | $158.2 \pm 7.9$ | $6.93 \pm 0.08$ | **0.360** | **0.174** |
| **10.877 s** | $43.77 \pm 0.06$ | $7.48 \pm 0.07$ | $2.78 \pm 0.03$ | $0.29 \pm 0.02$ | $2.33 \pm 0.02$ | $196.9 \pm 9.8$ | $7.40 \pm 0.08$ | **0.343** | **0.212** |

---

## Directory organisation

```
unified_fitting_ND-i-004/
├── README.md                                          # Comprehensive documentation (this file)
├── data/
│   ├── dataNDi004_fitting_results_composite.csv      # Full composite scattering curve (q in [0.002, 2.5] nm^-1)
│   ├── dataNDi004_fitting_results_short_distance.csv # 3.0 m camera wide CSV (50 frames)
│   ├── dataNDi004_fitting_results_long_distance.csv  # 31.0 m camera wide CSV (50 frames)
│   ├── frame_metadata.csv                            # Frame timings & exposure info
│   ├── short_distance/                               # 50 regular .dat files
│   ├── long_distance/                                # 50 regular .dat files
│   └── final_stage/                                  # 4 regular .dat files
├── models/
│   ├── __init__.py
│   ├── scattering_models.py                          # Fully vectorized analytical form & structure factors
│   └── fitting_engine.py                             # Non-linear least-squares engine with uncertainties
├── scripts/
│   ├── 00_prepare_data_and_composite.py              # Wide-format CSV generation & composite curve stitching
│   ├── 00_fit_composite_baseline.py                  # Baseline anchor fit across full q-range
│   ├── 01_fit_short_distance_models.py               # Time-resolved 3.0 m fit (ACC + Fluctuations)
│   ├── 02_fit_long_distance_models.py                # Time-resolved 31.0 m fit (Spinodal + Porod + ACC)
│   ├── 03_generate_all_frame_plots.py                # Parallel rendering of 100 diagnostic fit plots
│   ├── 04_compare_parameter_trajectories.py          # 4-panel cross-scale kinetic synthesis figure
│   ├── 05_build_unified_dashboard.py                 # Self-contained interactive HTML review dashboard
│   └── run_unified_pipeline.py                       # Master orchestrator script
└── output/
    ├── composite/                                    # Baseline anchor parameters & fit plot
    ├── short_distance/                               # parameters.csv, fitted_curves.csv, trajectories.png, fits/ (50 PNGs)
    ├── long_distance/                                # parameters.csv, fitted_curves.csv, trajectories.png, fits/ (50 PNGs)
    ├── unified_trajectories/                         # unified_cross_scale_trajectories.png (300 DPI)
    └── dashboard/                                    # index.html (Master interactive dashboard)
```

---

## Execution guide

To reproduce or re-run the entire pipeline from scratch:

```bash
cd "unified_fitting_ND-i-004"
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --all
```

To run individual steps:
```bash
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 0  # Prepare data & composite
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 1  # Composite baseline fit
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 2  # Short-distance fit (3.0 m)
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 3  # Long-distance fit (31.0 m)
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 4  # Render 100 diagnostic frame plots
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 5  # Cross-scale trajectory synthesis
/home/tomek/anaconda3/bin/python scripts/run_unified_pipeline.py --step 6  # Build master interactive dashboard
```
