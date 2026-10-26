# SAXS Scattering Model Formulations

The total coherent scattering intensity $I_{\text{tot}}(q)$ is evaluated as implemented across the analysis suite:

$$I_{\text{tot}}(q) = I_{\text{primary/fractal}}(q) + I_{\text{large}}(q) + I_{\text{high-q}}(q) + I_{\text{porod}}(q) + I_{\text{bkg}}$$

### 1. Primary Building Blocks & Mass Fractal Structure Factor
For spherical building blocks within a mass fractal network, the inter-particle structure factor $S(q)$ follows the Teixeira formulation:

$$S(q) = 1 + \frac{D_m \Gamma(D_m - 1)}{(q R)^{D_m}} \frac{\sin[(D_m - 1)\arctan(q\xi)]}{\left[ 1 + (q\xi)^{-2} \right]^{(D_m - 1)/2}}$$

For macroscopic, space-filling percolated gels ($\xi \ge 10^6\text{ nm}$), the analytical infinite-network limit is evaluated:

$$S(q) = 1 + \frac{D_m \Gamma(D_m - 1) \sin\left(\frac{\pi}{2}(D_m - 1)\right)}{(q R)^{D_m}}$$

### 2. Polydispersity Integration
For primary particles and large droplets, particle radii follow a lognormal probability density function $f_N(r; \mu, \sigma)$ converted from arithmetic mean $R_{\text{mean}}$ and standard deviation $\text{stdDev}$:

$$\mu = \ln\left( \frac{R_{\text{mean}}^2}{\sqrt{\text{stdDev}^2 + R_{\text{mean}}^2}} \right), \quad \sigma = \sqrt{\ln\left( 1 + \frac{\text{stdDev}^2}{R_{\text{mean}}^2} \right)}$$

The integrated intensity incorporates the square of particle volume:

$$I_{\text{fractal}}(q) = (\Delta \rho)^2 \int r^6 n(r) P(q, r) S(q, r, D_m, \xi) \, dr$$

where $(\Delta \rho)^2 = 26.884225\text{ cm}^{-2}$ is the calibrated electron contrast for amorphous calcium carbonate in water.

### 3. Macroscopic Droplets & Porod Tail
- **Large Droplets:** Monodisperse or lognormal spherical droplets ($\langle R_{\text{large}} \rangle \in [2.2, 3.5]\ \mu\text{m}$) capturing capillary-sheared emulsion domains.
- **Asymptotic Porod Tail:** $I_{\text{porod}}(q) = A_{\text{porod}} / q^4$ representing scattering from smooth interfaces of giant domains ($q R \gg 1$).

### 4. High-$q$ Viscoelastic Mesh Fluctuations (PNCs / Tanaka Framework)
The excess scattering at high $q$ arises from dynamic concentration and network density fluctuations of the transient viscoelastic mesh formed by entangled, hydrated ionic oligomers (prenucleation clusters / DOLLOPs), modeled via the generalized Fisher-Burford / Tanaka formulation:

$$I_{\text{high-q}}(q) \equiv I_{\text{fluc}}(q) = \frac{A_{\text{fluc}}}{\left[ 1 + \frac{2}{3\mu_{\text{fluc}}} (q \xi_{\text{fluc}})^2 \right]^{\mu_{\text{fluc}}/2}}$$

where:
- $\xi_{\text{fluc}} \in [1.5, 4.5]\text{ nm}$ is the characteristic correlation mesh size of the entangled PNC oligomer network.
- $A_{\text{fluc}}$ is the fluctuation amplitude scaling with the local density contrast of the transient strands.
- $\mu_{\text{fluc}} = 4.0$ corresponds to Porod-type boundary decay across viscoelastic strands (with $\mu_{\text{fluc}} = 2.0$ representing classical Ornstein-Zernike critical concentration fluctuations).


