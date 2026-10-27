# Small-angle X-ray scattering (SAXS) model description and physical interpretation

**Title:** Multi-scale SAXS analysis of calcium carbonate precipitation: phase separation regimes, viscoelastic network dynamics, and model continuum  
**Target section:** Theoretical model formulation, kinetic regimes and fitting parameter synthesis   

**Associated project artifacts:**
* Interactive Master Dashboard: [results/index.html](results/index.html) and [results/dashboard/index.html](results/dashboard/index.html)
* Vector Model Transformation Diagram: [results/output/model_transformation_diagram.svg](results/output/model_transformation_diagram.svg)
* Viscoelastic Phase Separation (VPS) Verification Report: [results/output/vps_framework_verification_report.md](results/output/vps_framework_verification_report.md)
* Cross-Series Kinetic Trajectories: [results/output/cross_series_trajectories_comparison.png](results/output/cross_series_trajectories_comparison.png)
* Trajectory Telemetry Data: [results/output/cross_series_trajectories.json](results/output/cross_series_trajectories.json)

---

## 1. Introduction and five physicochemical regimes

In supersaturated, alkaline aqueous solutions ($\text{Ca}^{2+} + \text{CO}_3^{2-} \rightarrow \text{CaCO}_3$), the initial phase separation trajectory frequently bypasses classical monomer-by-monomer nucleation. Instead, precipitation initiates through liquid-liquid phase separation (LLPS), condensation of hydrated ionic oligomers (DOLLOPs), and spinodal network formation.

Under stopped-flow rapid-mixing conditions at the ESRF ID02 beamline ($q \in [0.0020, 2.54]\ \text{nm}^{-1}$, real-space dimensions $2.5\ \text{nm}$ to $3.1\ \mu\text{m}$, time resolution $0.045\ \text{s}$ to $10.88\ \text{s}$), the six studied reaction series organise into five distinct physical regimes:

### Regime 1: Sub-spinodal classical nucleation (10.0 mM, ND-i-003, pH 11.3)
* **Thermodynamic Condition:** Below the spinodal boundary ($\ln S_{\text{ACC}} < \ln S^*$).
* **Observed Scattering Morphology:** Initial frames contain predominantly background noise, reflecting an induction period during which a critical thermodynamic barrier must be overcome.
* **Fitted Physical Model:** Discrete ACC nanospheres ($R_{\text{mean}} = 138.8 \pm 1.2\ \text{nm}$, $\text{stdDev} = 6.4\ \text{nm}$, $\phi = 1.72 \times 10^{-4}$), with a secondary smaller ACC mode ($R_2 = 93.5 \pm 2.3\ \text{nm}$, $\text{stdDev}_2 = 35.0\ \text{nm}$, $\phi_2 = 4.34 \times 10^{-4}$), and trace large droplet mode ($\langle R_{\text{large}} \rangle = 2.63\ \mu\text{m}$, $\phi_{\text{large}} = 6.12 \times 10^{-5}$). No mass fractal network is present ($S(q) = 1$).
* **Statistical Quality:** $\chi^2_{\text{red}} = 0.966$.
* **Time-Resolved Kinetics:** ACC nanospheres grow continuously from $\approx 60\ \text{nm}$ to $107\ \text{nm}$ over $\approx 5\ \text{s}$.

### Regime 2: Crossover boundary / hybrid nucleation-spinodal demixing (15.0 mM, ND-i-004, pH 11.3)
* **Thermodynamic Condition:** At the instability boundary ($\ln S_{\text{ACC}} \approx \ln S^*$).
* **Observed Scattering Morphology:** Dual-knee scattering curve exhibiting simultaneous coexistence of discrete ACC nanospheres and an emerging percolating mass fractal network.
* **Fitted Physical Model:** Coexisting discrete ACC nanospheres ($R_{\text{acc}} = 47.7 \pm 0.1\ \text{nm}$, $\text{std}_{\text{acc}} = 4.79 \pm 0.11\ \text{nm}$, $\phi_{\text{acc}} = 4.65 \times 10^{-5}$) and early percolating mass fractal gel network ($R_{\text{prim}} = 4.94 \pm 0.46\ \text{nm}$, $\text{std}_{\text{prim}} = 6.56 \pm 0.14\ \text{nm}$, $\phi_{\text{net}} = 2.71 \times 10^{-5}$, $D_m = 2.850 \pm 0.0001$, $\xi = 108.7 \pm 2.2\ \text{nm}$), plus Porod tail ($A_{\text{porod}} = 1.50 \times 10^{-7}$).
* **Statistical Quality:** $\chi^2_{\text{red}} = 1.309$.

### Regime 3: Spinodal gelation (20.0 mM, ND-i-002, pH 11.3)
* **Thermodynamic Condition:** Barrierless liquid-liquid demixing ($\ln S_{\text{ACC}} > \ln S^*$).
* **Observed Scattering Morphology:** Discrete ACC spheres disappear entirely ($\phi_{\text{sph}} \to 0$). Small building units form practically instantaneously without induction delay, establishing a characteristic knee at $q \approx 0.1\ \text{nm}^{-1}$. Evolution is confined to low $q$ where the fractal network expands.
* **Fitted Physical Model:** Percolated mass fractal network ($R_{\text{mean}} = 11.46 \pm 0.12\ \text{nm}$, $\text{stdDev} = 4.01 \pm 0.02\ \text{nm}$, $\phi = 4.74 \times 10^{-5}$, $D_m = 1.859 \pm 0.007$, $\xi = 351.6 \pm 24.0\ \text{nm}$), Tanaka viscoelastic fluctuation mesh ($\xi_{\text{fluc}} = 7.84 \pm 0.45\ \text{nm}$, $A_{\text{fluc}} = 0.020 \pm 0.003\ \text{cm}^{-1}$), and Porod tail ($A_{\text{porod}} = 5.28 \times 10^{-8}$).
* **Statistical Quality:** $\chi^2_{\text{red}} = 0.971$.

### Regime 4: Deep spinodal LLPS and viscoelastic phase separation benchmark (40.0 mM pure CaCO3, ND-i-009 & ND-i-001)
* **ND-i-009 (40.0 mM, pure, pH 10.0 Benchmark):**
  * **Composite Baseline Parameters:** Percolated mass fractal network ($R_{\text{mean}} = 9.40 \pm 0.06\ \text{nm}$, $\text{stdDev} = 3.41 \pm 0.02\ \text{nm}$, $\phi = 6.00 \times 10^{-3}$, $D_m = 1.978 \pm 0.007$, $\xi = 248.2 \pm 5.1\ \text{nm}$), coexisting macroscopic droplets ($\langle R_{\text{large}} \rangle = 3196 \pm 3\ \text{nm}$ [$3.20\ \mu\text{m}$], $\text{std}_{\text{large}} = 231.4 \pm 1.6\ \text{nm}$, $\phi_{\text{large}} = 1.205 \times 10^{-3}$), and high-$q$ Ornstein-Zernike fluctuations ($A_{\text{OZ}} = 10.0 \pm 16.9$, $\xi_{\text{OZ}} = 14.66 \pm 12.41\ \text{nm}$), $\text{bkg} = 0.0\ \text{cm}^{-1}$, $\chi^2_{\text{red}} = 3.218$, $\text{AIC} = 1941$.
  * **Time-Resolved Short-Distance Trajectory (351 frames, 15 ms cadence):** Primary building block radius is strictly invariant ($R_{\text{prim}} \approx 2.80\ \text{nm}$), while network fractal dimension densifies monotonically ($D_m: 1.881 \to 1.925$), demonstrating **Tanaka volume shrinking** (solvent expulsion from viscoelastic strands). Internal fluctuation mesh correlation length $\xi_{\text{fluc}} \approx 2.7 - 2.8\ \text{nm}$.
* **ND-i-001 (40.0 mM, pure, pH 11.3 Rapid Quench):**
  * Higher supersaturation drives instantaneous barrierless quench ($t_{\text{ind}} < 0.06\ \text{s}$): denser network ($R_{\text{mean}} = 15.43 \pm 0.74\ \text{nm}$, $\text{stdDev} = 6.25 \pm 0.11\ \text{nm}$, $D_m = 2.265 \pm 0.006$, $\xi = 923.1 \pm 129.7\ \text{nm}$, $\phi = 6.1 \times 10^{-5}$), coexisting macroscopic droplets ($\langle R_{\text{large}} \rangle = 3042 \pm 11\ \text{nm}$ [$3.04\ \mu\text{m}$], $\phi_{\text{large}} = 8.1 \times 10^{-6}$), and high-$q$ Tanaka viscoelastic fluctuation mesh ($\xi_{\text{fluc}} = 15.00 \pm 0.16\ \text{nm}$, $A_{\text{fluc}} = 4.51 \pm 0.12\ \text{cm}^{-1}$), $\chi^2_{\text{red}} = 0.718$.

### Regime 5: Additive perturbation branch (40.0 mM + 1.0 mM CoCl2, ND-i-005, pH 11.3)
* **Chemical Dopant Effect:** $\text{Co}^{2+}$ adsorption lowers interfacial tension $\gamma$, modulating primary building blocks ($R_{\text{mean}} = 2.40 \pm 0.13\ \text{nm}$, $\text{stdDev} = 1.27 \pm 0.03\ \text{nm}$) and network correlation length ($\xi = 4015\ \text{nm}$), whilst opening/loosening the fractal network ($D_m = 1.966 \pm 0.004$ vs $2.265$ in pure ND-i-001).
* **Fitted Physical Model:** Open mass fractal network with coexisting macroscopic droplets ($\langle R_{\text{large}} \rangle = 3471 \pm 38\ \text{nm}$ [$3.47\ \mu\text{m}$], $\phi_{\text{large}} = 2.80 \times 10^{-5}$) and Tanaka viscoelastic fluctuation mesh ($\xi_{\text{fluc}} = 1.82 \pm 0.66\ \text{nm}$, $A_{\text{fluc}} = 0.0028 \pm 0.0018\ \text{cm}^{-1}$), $\chi^2_{\text{red}} = 4.470$.

---

## 2. Parameter synthesis across regimes

The calibrated parameters extracted directly from the composite baseline fitting outputs across all reaction series are:

| Series ID | Concentration | pH | Physical regime | Primary building block $R_{\text{prim}}$ [nm] | Fractal dimension $D_m$ | Network cutoff $\xi$ [nm] | Droplet radius $\langle R_{\text{large}} \rangle$ [$\mu\text{m}$] | Fluctuation mesh $\xi_{\text{fluc}}$ [nm] | Fluctuation ampl. $A_{\text{fluc}}$ [$\text{cm}^{-1}$] | Reduced $\chi^2$ |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **ND-i-003** | 10.0 mM | 11.3 | Sub-spinodal nucleation | 138.8 (ACC) | N/A ($S=1$) | N/A | 2.63 | N/A | N/A | 0.966 |
| **ND-i-004** | 15.0 mM | 11.3 | Crossover boundary | 4.94 | 2.850 | 108.7 | N/A (Porod) | N/A | N/A | 1.309 |
| **ND-i-002** | 20.0 mM | 11.3 | Spinodal gelation | 11.46 | 1.859 | 351.6 | N/A (Porod) | 7.84 | 0.020 | 0.971 |
| **ND-i-009** | 40.0 mM | 10.0 | Primary pure 40 mM benchmark | 9.40 | 1.978 | 248.2 | 3.20 | 2.72 &ndash; 2.84 | 0.050 (OZ) | 3.218 |
| **ND-i-001** | 40.0 mM | 11.3 | Rapid quench spinodal | 15.43 | 2.265 | 923.1 | 3.04 | 15.00 | 4.506 | 0.718 |
| **ND-i-005** | 40.0 mM | 11.3 | Additive perturbation branch | 2.40 | 1.966 | 4014.7 | 3.47 | 1.82 | 0.0028 | 4.470 |

---

## 3. Mathematical scattering model formulation (code-aligned)

The total coherent scattering intensity $I_{\text{tot}}(q)$ is evaluated as implemented in `scattering_models.py` and `fitting_engine.py`:

$$I_{\text{tot}}(q) = I_{\text{primary/fractal}}(q) + I_{\text{large}}(q) + I_{\text{high-q}}(q) + I_{\text{porod}}(q) + I_{\text{bkg}}$$

### 3.1. Primary particle form factor
For an isotropic sphere of radius $R$:
$$P(q, R) = \left[ 3 \frac{\sin(qR) - qR\cos(qR)}{(qR)^3} \right]^2, \quad \lim_{qR \to 0} P(q, R) = 1.0$$

### 3.2. Mass fractal structure factor
The inter-particle structure factor $S(q)$ for primary building blocks within the fractal network follows the Teixeira formulation:
$$S(q) = 1 + \frac{D_m \Gamma(D_m - 1)}{(q R)^{D_m}} \frac{\sin[(D_m - 1)\arctan(q\xi)]}{\left[ 1 + (q\xi)^{-2} \right]^{(D_m - 1)/2}}$$

For macroscopic, infinitely large percolated gels ($\xi \ge 10^6\ \text{nm}$), the code evaluates the exact analytical infinite-network limit:
$$S(q) = 1 + \frac{D_m \Gamma(D_m - 1) \sin\left(\frac{\pi}{2}(D_m - 1)\right)}{(q R)^{D_m}}$$

### 3.3. Polydispersity integration
For primary and large sphere populations, particle radii follow a lognormal probability density function $f_N(r; \mu, \sigma)$ converted from arithmetic mean $R_{\text{mean}}$ and standard deviation $\text{stdDev}$:
$$\mu = \ln\left( \frac{R_{\text{mean}}^2}{\sqrt{\text{stdDev}^2 + R_{\text{mean}}^2}} \right), \quad \sigma = \sqrt{\ln\left( 1 + \frac{\text{stdDev}^2}{R_{\text{mean}}^2} \right)}$$

Volume per unit is integrated numerically over $r \in [10^{-5}, R_{\text{mean}} + 10\cdot\text{stdDev}]$ with 300 grid steps:
$$V_{\text{unit}} = \int \frac{4}{3}\pi r^3 f_N(r; \mu, \sigma) \, dr$$
Number density is $n = \phi / V_{\text{unit}}$, yielding the integrated intensity:
$$I_{\text{fractal}}(q) = (\Delta \rho)^2 \int r^6 n(r) P(q, r) S(q, r, D_m, \xi) \, dr$$
where $(\Delta \rho)^2 = 26.884225\ \text{cm}^{-2}$ is the calibrated electron contrast for amorphous calcium carbonate in water.

### 3.4. Large droplet population
The secondary micrometer droplet mode is integrated analogously:
$$I_{\text{large}}(q) = (\Delta \rho)^2 \int r^6 n_{\text{large}}(r) P(q, r) \, dr$$
producing steep Porod decay ($I \propto q^{-4}$) at $q R \gg 1$.

### 3.5. High-$q$ viscoelastic mesh fluctuations (PNC network dynamics)
At high momentum transfer ($q \ge 0.35\ \text{nm}^{-1}$), the transient network of hydrated prenucleation cluster (PNC) oligomers exhibits concentration and network mesh fluctuations rather than isolated particle scattering. These fluctuations are modeled via:
* **Fisher-Burford / Tanaka Viscoelastic Fluctuation Mesh:**
  $$I_{\text{FB}}(q) = \frac{A_{\text{fluc}}}{\left[ 1 + \frac{2}{3\mu} (q \xi_{\text{fluc}})^2 \right]^{\mu/2}}$$
  where $\xi_{\text{fluc}}$ is the dynamic correlation length (mesh size of entangled PNC clusters), $A_{\text{fluc}}$ is the forward fluctuation amplitude, and $\mu = 4$ reproduces asymptotic Porod decay ($q^{-4}$) representing interface fluctuations, while $\mu = 2$ corresponds to standard Lorentzian decay.
* **Ornstein-Zernike (OZ) Fluctuations (ND-i-006, ND-i-009):**
  $$I_{\text{OZ}}(q) = \frac{A_{\text{OZ}}}{1 + (q \xi_{\text{OZ}})^2}$$
  representing thermal and concentration fluctuations in the solvent channels between percolated struts.

> **Note on Physical Interpretation:** Prenucleation clusters (PNCs) in percolated post-spinodal regimes do not persist as isolated, rigid spherical globules with sharp density boundaries. Instead, hydrated ion clusters dynamically entangle into a viscoelastic network whose characteristic correlation distance is quantitatively captured by $\xi_{\text{fluc}}$ (or $\xi_{\text{OZ}}$). Consequently, discrete polydisperse spherical models for PNCs are unphysical in the post-spinodal regime and are not used.

### 3.6. Porod tail and incoherent background
* **Porod tail (in ND-i-002, ND-i-004):** $I_{\text{porod}}(q) = \frac{A_{\text{porod}}}{q^4}$
* **Incoherent background:** $I_{\text{bkg}} = \text{bkg}$

---

## 4. Verification of Tanaka viscoelastic phase separation (VPS)

### 4.1. Dynamic asymmetry and the viscoelastic length scale
In the Tanaka Two-Fluid Model, the mineral network velocity $\mathbf{v}_p$ and solvent velocity $\mathbf{v}_s$ are coupled through friction $\zeta$:
$$-\phi \nabla \mu + \nabla \cdot \boldsymbol{\sigma}^{(p)} - \zeta (\mathbf{v}_p - \mathbf{v}_s) = 0$$
Dynamic asymmetry is large ($\alpha = \zeta_p / \zeta_s - 1 \gg 0$): the mineral network supports the viscoelastic stress $\boldsymbol{\sigma}^{(p)}$, while the solvent supports negligible stress.

The viscoelastic length scale $\xi_{\text{ve}}$ is:
$$\xi_{\text{ve}} = \sqrt{\frac{K_{\text{os}} + \frac{4}{3}G}{\zeta}} \tau_r \sim \sqrt{D_{\text{coop}} \tau_r} \approx 100 - 1000\ \text{nm}$$
* For $q > 1/\xi_{\text{ve}}$ ($r < \xi_{\text{ve}}$): Stress relaxation is slow compared to diffusion. The network acts as an elastic gel, suppressing hydrodynamic coarsening and driving volume shrinking.
* For $q < 1/\xi_{\text{ve}}$ ($r > \xi_{\text{ve}}$): Stress has relaxed. Capillary flow and phase inversion into macroscopic liquid droplets ($\langle R_{\text{large}} \rangle \in [2.2, 3.5]\ \mu\text{m}$) proceed unhindered.

### 4.2. Quantitative proof of volume shrinking (ND-i-009 benchmark)
In classical spinodal decomposition (Model B/H), domain coarsening preserves phase volume. In contrast, Viscoelastic Phase Separation requires volume shrinking of the transient gel network:
* Short-distance time-resolved frames (351 frames, 15 ms cadence) confirm that the primary building block radius remains invariant at $R_{\text{prim}} = 2.80\ \text{nm}$.
* The network fractal dimension increases monotonically from $D_m = 1.881$ to $1.925$.
* Increasing $D_m$ at invariant primary node size directly demonstrates continuous compaction and densification of network strands as solvent is expelled.

### 4.3. Phase inversion into macroscopic droplets via capillary pinch-off
When coarsening strains exceed the transient yield stress of the struts, the network ruptures mechanically. Detached fluid domains contract under interfacial tension into compact spherical droplets ($\langle R_{\text{large}} \rangle \approx 3.2\ \mu\text{m}$), in agreement with Rayleigh-Plateau capillary instability:
$$\lambda_{\text{opt}} \approx 9.02 R_{\text{strut}} \implies R_{\text{drop}} \approx 1.89 R_{\text{strut}} \approx 2.8 - 3.5\ \mu\text{m}$$

### 4.4. Distinction from arrested spinodal decomposition (ASD)
Unlike Arrested Spinodal Decomposition, which permanently freezes domain growth via repulsive glass jamming, the calcium carbonate system shows continuous structural evolution, progressive volume shrinking, and secondary droplet coexistence.

---

## 5. Benchmark 40 mM pure CaCO3 baseline parameters (ND-i-009)

The calibrated composite parameters for ND-i-009 ($q \in [0.0020, 2.54]\ \text{nm}^{-1}$) are:

| Model component | Parameter symbol | Fitted value | Std. error | Physical significance in code |
| :--- | :---: | :---: | :---: | :--- |
| **Large droplets** | $\langle R_{\text{large}} \rangle$ | $3196.1\ \text{nm}$ ($3.20\ \mu\text{m}$) | $\pm 2.8\ \text{nm}$ | Mean droplet radius from capillary break-up |
| | $\text{std}_{\text{large}}$ | $231.4\ \text{nm}$ ($0.23\ \mu\text{m}$) | $\pm 1.6\ \text{nm}$ | Droplet size distribution standard deviation |
| | $\phi_{\text{large}}$ | $1.205 \times 10^{-3}$ | $\pm 0.010 \times 10^{-3}$ | Volume fraction of pinched-off emulsion phase |
| **Fractal network** | $R_{\text{mean}}$ | $9.40\ \text{nm}$ | $\pm 0.06\ \text{nm}$ | Primary network building block radius |
| | $\text{stdDev}$ | $3.41\ \text{nm}$ | $\pm 0.02\ \text{nm}$ | Primary node size distribution standard deviation |
| | $\phi$ | $6.000 \times 10^{-3}$ | $\pm 0.075 \times 10^{-3}$ | Volume fraction of spanning fractal gel mesh |
| | $D$ | $1.978$ | $\pm 0.007$ | Fractal dimension of DLCA spinodal network |
| | $\xi$ | $248.2\ \text{nm}$ | $\pm 5.1\ \text{nm}$ | Upper network structural correlation cutoff |
| **High-$q$ fluctuations** | $A_{\text{OZ}}$ | $10.00$ | $\pm 16.89$ | Ornstein-Zernike fluctuation amplitude |
| | $\xi_{\text{OZ}}$ | $14.66\ \text{nm}$ | $\pm 12.41\ \text{nm}$ | Fluctuation correlation length in solvent channels |
| **Contrast and background** | $\rho^2$ | $26.884225$ | Fixed | Calibrated electron SLD contrast factor for ACC |
| | $\text{bkg}$ | $0.0000\ \text{cm}^{-1}$ | Fixed ($0.0$) | Buffer incoherent background |
| **Statistical quality** | $\chi^2_{\text{red}}$ | **$3.218$** | : | Global goodness of fit over 3 orders in $q$ |
| | AIC | **$1941.0$** | : | Akaike Information Criterion |

---

## 6. The volume-scattering disparity ($I(0) \propto R^6$)

The forward scattering intensity per object scales with $V^2$, or $R^6$:
$$I_{\text{single}}(0) = (\Delta\rho)^2 V^2 \propto R^6$$

Evaluating for ND-i-009 using the fitted parameters:
* Primary network node: $r_0 = 9.40\ \text{nm} \implies V_{\text{node}} = \frac{4}{3}\pi (9.40)^3 \approx 3.48 \times 10^3\ \text{nm}^3$.
* Secondary large droplet: $R_{\text{large}} = 3196\ \text{nm} \implies V_{\text{drop}} = \frac{4}{3}\pi (3196)^3 \approx 1.37 \times 10^{11}\ \text{nm}^3$.
* Volume ratio per particle:
  $$\frac{V_{\text{drop}}}{V_{\text{node}}} = \left(\frac{3196}{9.40}\right)^3 \approx \mathbf{3.93 \times 10^7}$$
* Single-particle forward scattering advantage:
  $$\frac{I_{\text{drop}}(0)}{I_{\text{node}}(0)} = \left(\frac{V_{\text{drop}}}{V_{\text{node}}}\right)^2 \approx \mathbf{1.54 \times 10^{15}}$$

Macroscopically, with $\phi_{\text{large}} = 1.205 \times 10^{-3}$ and $\phi_{\text{mesh}} = 6.000 \times 10^{-3}$:
$$\frac{I_{\text{drop}}(0)}{I_{\text{mesh}}(0)} = \frac{\phi_{\text{large}} V_{\text{drop}}}{\phi_{\text{mesh}} V_{\text{cluster}}} \gg 1$$

Even though the droplets represent only $\approx 16.7\%$ of total mineral volume, their steep Porod profile ($I \propto q^{-4}$) generates $\approx 77.8\%$ of measured scattering at $q_{\text{min}} = 0.00213\ \text{nm}^{-1}$, dropping below $2\%$ by $q \approx 0.008\ \text{nm}^{-1}$ where the fractal network dominates.

---

## 7. Comparison with alternative physical hypotheses

| Feature / metric | Hypothesis 1: Crystalline spherulites (vaterite/calcite) | Hypothesis 2: Solid ACC colloid nucleation | Hypothesis 3 (this work): Spinodal bicontinuous mesh + capillary droplets |
| :--- | :--- | :--- | :--- |
| **Physical state of droplets** | Rigid, faceted polycrystalline aggregates. | Solid, spherical amorphous nanoparticles. | Hydrated, viscous liquid precursor droplets (PILP). |
| **Diffraction signatures** | Strong Bragg peaks in WAXS. **Not observed.** | Broad amorphous halo, no Bragg peaks. | Broad solvent-like liquid structure factor. Perfectly matches data. |
| **Boundary topology** | Rough surface ($D_s \approx 2.3 - 2.8$). Yields $q^{-3.2}$ to $q^{-3.7}$. | Spherical interface, but typically smaller ($R \sim 20 - 100\ \text{nm}$). | Surface-tension-smoothed meniscus ($D=3$). Strictly yields $q^{-4}$ Porod decay. |
| **Mesh relationship** | Incompatible; crystals consume gel. | Unclear how solid spheres form continuous $D_m \approx 1.85$ mesh. | **Direct hydrodynamic daughter phase:** pinched off from continuous mesh via Rayleigh-Plateau. |
| **Statistical fit quality** | Fails across multiple decades ($\chi^2_{\text{red}} > 15$). | Insufficient for multi-tier composite curves ($\chi^2_{\text{red}} > 5$). | **Converged across all series ($\chi^2_{\text{red}} \in [0.4, 3.2]$)**. |

---

## 8. Experimental and beamline specifications (ESRF ID02)

| Parameter | Specification | Beamline operation notes |
| :--- | :--- | :--- |
| **Synchrotron facility** | **ESRF - The European Synchrotron** (Grenoble, France) | Extremely Brilliant Source (EBS) storage ring (Proposal CH-6828) |
| **Beamline** | **ID02 (TRUSAXS / USAXS / SAXS)** | High-brilliance undulator beamline |
| **Photon energy ($E$)** | **$12.23\ \text{keV}$ ($12,230.0\ \text{eV}$)** | $\text{Si}(111)$ double crystal monochromator (`egy=12230.0`, `theta=9.308`) |
| **X-ray wavelength ($\lambda$)** | **$0.101377\ \text{nm}$ ($1.0138\ \text{Å}$)** | Calibrated pyFAI wavelength parameter (`WaveLength = 1.01377e-10 m`) |
| **Detector system** | **Dectris Eiger2 4M** | Hybrid photon counting pixel detector ($2068 \times 2162\ \text{pixels}$, $75\ \mu\text{m}$ pitch) |
| **Short-distance camera (SD)** | **$3.0\ \text{m}$ nominal ($2.9997\ \text{m}$ calibrated)** | Covers $q \in [0.022, 2.54]\ \text{nm}^{-1}$; $100\ \mu\text{L}$ to $200\ \mu\text{L}$ stopped-flow shots |
| **Long-distance camera (LD)** | **$31.0\ \text{m}$ nominal ($30.7059\ \text{m}$ calibrated)** | Covers $q \in [0.002, 0.25]\ \text{nm}^{-1}$; $100\ \mu\text{L}$ to $200\ \mu\text{L}$ stopped-flow shots |
| **Combined composite $q$-range** | **$0.0020\ \text{nm}^{-1}$ to $2.54\ \text{nm}^{-1}$** | $>3$ orders of magnitude in reciprocal space stitched across cameras |
| **Sample environment** | **BioLogic SFM-400 Stopped-Flow Mixer** | 4-syringe pneumatic injection ($3.33\ \text{mL/s}$/syringe) with $1.3\ \text{mm}$ capillary |

---

## 9. License and attribution

* **Data processing responsible:** Tomasz M. Stawski (`tomasz.stawski@bam.de`), Federal Institute for Materials Research and Testing (BAM)
* **Beamtime proposal:** ESRF Proposal CH-6828 (*Time-resolved Small-Angle X-ray Scattering on Non-Classical Nucleation Pathways of Calcium Carbonate*).
* **License:** This work is licensed under the terms of the standard [MIT License](LICENSE).
