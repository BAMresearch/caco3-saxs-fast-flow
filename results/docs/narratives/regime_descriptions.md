# Logical Transition Continuum Across Experimental Regimes

Across reactant supersaturation and chemical conditions, the six experimental series form a continuum spanning five distinct thermodynamic regimes. Throughout these descriptions, three distinct length scales are differentiated:
1. **Primary Building Block Radius ($R_{\text{prim}}$):** The radius of the structural colloidal nodes forming the struts of the mass fractal network, probed by the primary form-factor knee ($P(q, R_{\text{prim}})$).
2. **Entangled PNC Fluctuation Mesh ($\xi_{\text{fluc}}$):** The internal correlation length of dynamic concentration fluctuations within the entangled, hydrated ionic oligomer (PNC / DOLLOP) mesh inside the strands, modeled by the generalized Fisher-Burford term ($I_{\text{fluc}}(q)$).
3. **Mass Fractal Network Cutoff Length ($\xi$):** The overall spatial correlation extent of the percolating cluster network, modeled via the Teixeira structure factor ($S(q, D_m, \xi)$).

---

### 1. Regime 1: Sub-spinodal Nucleation (10.0 mM, ND-i-003, pH 11.3)
- **Thermodynamic condition:** $\ln S_{\text{ACC}} < \ln S^*$, below the spinodal boundary.
- **Physical model:** Classical nucleation of discrete amorphous calcium carbonate (ACC) nanospheres ($R_{\text{acc}} \approx 64 \to 102\text{ nm}$ over 10 s, composite $138.8\text{ nm}$) with $S(q) = 1$ (absence of a percolated mass fractal network).
- **High-$q$ behavior:** Pure flat solvent background without sub-10 nm primary nodes or detectable PNC fluctuation mesh. Intermediate $q$ displays $q^{-4}$ Porod decay from smooth spherical interfaces coexisting with macroscopic surface scattering ($A_{\text{porod}}/q^4$).
- **Transition trigger:** Increasing supersaturation towards $C^* \approx 15\text{ mM}$ lowers the classical nucleation barrier until barrierless spinodal demixing becomes thermodynamically accessible.

---

### 2. Regime 2: Crossover Boundary (15.0 mM, ND-i-004, pH 11.3)
- **Thermodynamic condition:** $\ln S_{\text{ACC}} \approx \ln S^*$, directly at the threshold of thermodynamic and mechanical instability.
- **Physical model:** Scattering curve displaying simultaneous coexistence of discrete ACC nanospheres ($\langle R_{\text{acc}} \rangle \approx 15.9 \to 43.8\text{ nm}$, $\sigma \approx 7.5\text{ nm}$) and an emerging percolating mass fractal gel network ($R_{\text{prim}} \approx 3.5 - 4.9\text{ nm}$, $D_m \approx 2.48 \to 2.33$, $\xi \approx 20 \to 197\text{ nm}$).
- **High-$q$ behavior:** Scattering in the high-$q$ regime is dominated by the form factor of the discrete ACC nanospheres and the $R_{\text{prim}} \approx 3.5\text{ nm}$ building block nodes, with PNC mesh fluctuations being less pronounced.
- **Transition trigger:** Suppression of discrete ACC droplet nucleation as barrierless spinodal demixing takes over the entire solution volume.

---

### 3. Regime 3: Spinodal Gelation (20.0 mM, ND-i-002, pH 11.3)
- **Thermodynamic condition:** $\ln S_{\text{ACC}} > \ln S^*$, barrierless liquid-liquid phase separation (LLPS).
- **Physical model:** Discrete ACC nanospheres disappear entirely ($I_{\text{acc}} \to 0$). Continuous mass fractal gel network forms practically instantaneously ($D_m \approx 1.30 \to 1.97$, composite $D_m = 1.859$, $\xi \approx 350 - 400\text{ nm}$) coexisting with macroscopic phase-separated droplets ($\langle R_{\text{large}} \rangle \approx 5.2 \to 5.9\ \mu\text{m}$).
- **Structural node vs. fluctuation mesh:**
  - Primary building block nodes: $R_{\text{prim}} \approx 8.5 \to 11.5\text{ nm}$ (mean $10.65\text{ nm}$).
  - Entangled PNC fluctuation mesh: $\xi_{\text{fluc}} \approx 13.5 \to 15.0\text{ nm}$ ($A_{\text{fluc}} \approx 0.03 \to 0.16\text{ cm}^{-1}$), capturing the high-$q$ excess.
- **Transition trigger:** Further increase in supersaturation accelerates phase separation kinetics, requiring higher time resolution to resolve the early pre-percolation.

---

### 4. Regime 4: Deep Spinodal LLPS and VPS Benchmark (40.0 mM Pure CaCO3)
- **Thermodynamic condition:** Deep spinodal quench with extreme dynamic asymmetry between ionic oligomers and solvent.
- **ND-i-009 Benchmark (pH 10.0, 351 frames, 15 ms cadence):**
  - Resolves time-resolved pre-percolation ($\xi: 46.7\text{ nm} \to \infty$).
  - **Proof of Tanaka volume shrinking:** Monotonic densification of the mass fractal network ($D_m: 1.881 \to 1.925$) at a strictly invariant primary building block radius ($R_{\text{prim}} \approx 8.35 \to 8.69\text{ nm}$, mean $8.34\text{ nm}$, composite $9.40\text{ nm}$).
  - **Entangled PNC fluctuation mesh:** High-$q$ excess captures the internal viscoelastic correlation mesh $\xi_{\text{fluc}}(t)$, which drifts smoothly from $2.72\text{ nm}$ to $2.84\text{ nm}$ (mean $2.76\text{ nm}$), confirming continuous compaction of the entangled DOLLOP oligomers inside the network strands.
  - Secondary macroscopic liquid droplets coarsen continuously ($\langle R_{\text{large}} \rangle = 2.27 \to 3.20\ \mu\text{m}$) via capillary stress relaxation and strut rupture.
- **ND-i-001 Quench (pH 11.3, 50 frames):**
  - Higher supersaturation produces instantaneous barrierless quench ($t_{\text{ind}} < 0.06\text{ s}$).
  - Dense network ($D_m \approx 1.86 \to 1.98$, composite $2.265$, $\xi \approx 735 - 923\text{ nm}$) with primary building blocks $\langle R_{\text{prim}} \rangle \approx 6.04 \to 6.85\text{ nm}$ (mean $6.39\text{ nm}$).
  - Internal PNC fluctuation mesh: $\xi_{\text{fluc}} \approx 3.03 \to 3.06\text{ nm}$.
  - Coexisting macroscopic droplets: $\langle R_{\text{large}} \rangle \approx 3.28 \to 3.91\ \mu\text{m}$.

---

### 5. Regime 5: Additive Perturbation Branch (40.0 mM + 1.0 mM Co2+, ND-i-005, pH 11.3)
- **Thermodynamic condition:** Same nominal supersaturation as ND-i-001, perturbed by a transition metal ion ($\text{Co}^{2+}$).
- **Physical model:** Specific coordination of $\text{Co}^{2+}$ into nascent mineral clusters alters local interfacial energetics and network relaxation.
  - **Primary structural nodes:** $\langle R_{\text{prim}} \rangle \approx 6.91 \to 10.18\text{ nm}$ (mean $8.25\text{ nm}$, vs $6.39\text{ nm}$ in pure ND-i-001).
  - **Disruption of PNC fluctuation mesh:** Dramatically expands the internal oligomer fluctuation correlation length to $\xi_{\text{fluc}} \approx 8.04 \to 9.06\text{ nm}$ (mean $8.65\text{ nm}$, vs $3.04\text{ nm}$ in pure ND-i-001), reflecting significant steric/coordination swelling of the hydrated ionic mesh.
  - **Extended mesoscale correlation:** Extends the network correlation length to $\xi \approx 1000 - 1100\text{ nm}$ with fractal dimension $D_m \approx 2.15 \to 2.39$.
  - **Macroscopic droplets:** Preserves phase-inverted droplet coexistence ($\langle R_{\text{large}} \rangle \approx 3.1 \to 3.8\ \mu\text{m}$).
