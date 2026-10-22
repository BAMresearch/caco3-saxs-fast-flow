#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
Fitting script for a single population of polydisperse spherical ACC nanoparticles (log-normal distribution)
Tailored for the "10 mM pH 11.3 ND-i-003 selected" dataset (low supersaturation regime).

In this low supersaturation scenario:
- Spinodal decomposition and mass-fractal percolation are suppressed.
- No primary nodes or continuous gel mesh form (short distance data is at water background level).
- Precipitation proceeds via classical nucleation and growth of discrete, compact amorphous calcium carbonate (ACC) nanospheres.
- The scattering profiles in the long camera distance exhibit Porod decay (I ~ q^-4) and a Guinier knee corresponding to growing spheres (<R> ~ 75-105 nm).
"""

import os
import glob
from datetime import datetime
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from lmfit import Parameters, minimize
import matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_cache'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ID02 standard 50-frame time series timestamps (seconds)
ID02_TIMESTAMPS = np.array([
    0.060, 0.075, 0.091, 0.108, 0.126, 0.146, 0.167, 0.190, 0.214, 0.241,
    0.269, 0.300, 0.334, 0.370, 0.410, 0.453, 0.499, 0.550, 0.606, 0.667,
    0.733, 0.806, 0.885, 0.971, 1.066, 1.169, 1.282, 1.405, 1.539, 1.685,
    1.844, 2.017, 2.205, 2.410, 2.632, 2.873, 3.136, 3.421, 3.731, 4.068,
    4.433, 4.830, 5.261, 5.728, 6.235, 6.785, 7.381, 8.028, 8.729, 10.877
])


class SinglePolydisperseSphereFitter:
    """
    Encapsulates SAXS fitting for discrete spherical ACC nanoparticles
    with a log-normal size distribution.
    """

    def __init__(self, project_name, data_directories, q_range=(0.002, 0.10)):
        self.project_name = project_name
        self.data_directories = data_directories
        self.q_range = q_range

        self.time_stamp = datetime.now().strftime("%d%m%Y_%H%M%S")
        self.output_path = os.path.join(self.project_name, 'output_fits', self.time_stamp)
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(os.path.join(self.output_path, 'frame_fits'), exist_ok=True)
        print(f"Output directory: {self.output_path}")

        self.data = None
        self.file_names = []
        self.time_seconds = []
        self.params = None
        self.fit_results_df = None
        self.fitted_curves_df = None

    # --- 1. Data Loading ---

    def load_data(self):
        """
        Loads data from specified directories. Filters out short-distance files if mixed,
        and interpolates long-distance / final stage data onto a unified q-grid.
        """
        print("Loading and preparing SAXS profiles...")
        all_data = []
        q_min, q_max = -np.inf, np.inf

        for directory in self.data_directories:
            if not os.path.exists(directory):
                print(f"Warning: Directory '{directory}' does not exist. Skipping.")
                continue

            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".dat"):
                    # Exclude short-distance files (e.g. 00176 or 00177) that start at q > 0.02
                    # because they have no particle signal and restrict the q-range
                    filepath = os.path.join(directory, filename)
                    try:
                        df = pd.read_csv(filepath, sep=r'\s+', header=None, names=["q", "I", "dI"], skiprows=1)
                        df["q"] = pd.to_numeric(df["q"], errors="coerce")
                        df["I"] = pd.to_numeric(df["I"], errors="coerce")
                        df["dI"] = pd.to_numeric(df["dI"], errors="coerce")
                        df = df.dropna()

                        if len(df) > 100:
                            # Verify this is a long-distance file (q_min < 0.005 nm^-1)
                            if df["q"].min() < 0.005:
                                all_data.append((filename, df))
                                q_min = max(q_min, df["q"].min())
                                q_max = min(q_max, df["q"].max())
                            else:
                                print(f"Skipping short-distance file: {filename} (q_min = {df['q'].min():.3e})")
                    except Exception as e:
                        print(f"Error reading {filename}: {e}")

        if not all_data:
            print("Error: No valid long-distance *.dat files found.")
            return

        print(f"Loaded {len(all_data)} long-distance profiles. Common q range: [{q_min:.4e}, {q_max:.4e}] nm^-1")
        common_q = np.geomspace(q_min, min(q_max, 0.20), 800)

        data_dict = {'q': common_q}
        self.file_names = []
        self.time_seconds = []

        for idx, (filename, df) in enumerate(all_data):
            f_I = interpolate.interp1d(df["q"], df["I"], kind='linear', bounds_error=False, fill_value=np.nan)
            f_dI = interpolate.interp1d(df["q"], df["dI"], kind='linear', bounds_error=False, fill_value=np.nan)

            interp_I = f_I(common_q)
            interp_dI = np.maximum(f_dI(common_q), 1e-4)

            data_dict[str(idx)] = interp_I
            data_dict['error_' + str(idx)] = interp_dI
            self.file_names.append(filename)

            # Assign physical time
            if idx < len(ID02_TIMESTAMPS):
                self.time_seconds.append(ID02_TIMESTAMPS[idx])
            else:
                # Final stages: later time points (e.g. 60s, 120s)
                extra_idx = idx - len(ID02_TIMESTAMPS) + 1
                self.time_seconds.append(ID02_TIMESTAMPS[-1] + extra_idx * 30.0)

        self.data = pd.DataFrame(data_dict).dropna().reset_index(drop=True)
        self.data.to_csv(os.path.join(self.output_path, f"{self.project_name}_combined_data.csv"), index=False)
        print(f"Data prepared successfully ({len(self.file_names)} profiles, {len(self.data)} q-points).")

    # --- 2. Model Formulation ---

    @staticmethod
    def sphere_form_factor(qr):
        """ Normalized hard sphere form factor P(q, R). """
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (3.0 * (np.sin(qr) - qr * np.cos(qr)) / (qr**3))**2
        return np.where(qr < 1e-5, 1.0, P)

    @staticmethod
    def log_normal_pdf(r, mu, sigma):
        """ Lognormal probability density function. """
        with np.errstate(divide='ignore', invalid='ignore'):
            pdf = (1.0 / (r * sigma * np.sqrt(2.0 * np.pi))) * np.exp(-((np.log(r) - mu)**2) / (2.0 * sigma**2))
        return np.where(r <= 0, 0.0, np.nan_to_num(pdf, nan=0.0))

    @staticmethod
    def calc_sphere_intensity(q, Rmean, stdDev, phi, bkg, rho2=1.0, steps=250):
        """
        Calculates intensity for a lognormal distribution of spheres:
        I(q) = rho2 * int [ r^6 * n(r) * P(q, r) ] dr + bkg
        """
        if phi <= 0 or Rmean <= 0:
            return np.full_like(q, max(0.0, bkg))

        # Convert arithmetic mean and std to lognormal parameters
        mu = np.log(Rmean**2 / np.sqrt(stdDev**2 + Rmean**2))
        sigma = np.sqrt(np.log(1.0 + (stdDev**2 / Rmean**2)))

        r_max = Rmean + 8.0 * stdDev
        r_grid = np.linspace(1e-4, r_max, steps)
        dr = r_grid[1] - r_grid[0]

        # Unit volume
        vol_unit = np.sum((4.0 / 3.0) * np.pi * (r_grid**3) * SinglePolydisperseSphereFitter.log_normal_pdf(r_grid, mu, sigma) * dr)
        vol_unit = max(vol_unit, 1e-12)
        number_density = phi / vol_unit

        qr = np.outer(r_grid, q)
        P = SinglePolydisperseSphereFitter.sphere_form_factor(qr)
        dist = SinglePolydisperseSphereFitter.log_normal_pdf(r_grid, mu, sigma) * number_density

        integrand = (r_grid[:, None]**6) * dist[:, None] * P * dr
        np.nan_to_num(integrand, copy=False, nan=0.0)

        return rho2 * np.sum(integrand, axis=0) + bkg

    def _model_func(self, params, q):
        Rmean = params['Rmean'].value
        stdDev = params['stdDev'].value
        phi = params['phi'].value
        bkg = params['bkg'].value
        rho2 = params['rho2'].value if 'rho2' in params else 1.0
        return self.calc_sphere_intensity(q, Rmean, stdDev, phi, bkg, rho2=rho2)

    def _residuals(self, params, q, data, uncertainty):
        model = self._model_func(params, q)
        return (data - model) / uncertainty

    # --- 3. Parameter Setup ---

    def setup_parameters(self):
        """ Initializes lmfit parameters with realistic physical bounds for ACC nanoparticles. """
        params = Parameters()
        params.add('rho2', value=1.0, vary=False)
        params.add('phi', value=5.0e-4, min=1.0e-12, max=0.01, vary=True)
        # Rmean: realistic radius bounds for ACC nanospheres (10 to 250 nm)
        params.add('Rmean', value=85.0, min=10.0, max=250.0, vary=True)
        # stdDev: lognormal standard deviation (1 to 100 nm)
        params.add('stdDev', value=28.0, min=1.0, max=100.0, vary=True)
        params.add('bkg', value=0.01, min=0.0, max=1.0, vary=True)

        self.params = params
        print("Parameters initialized: Rmean in [10, 250] nm, stdDev in [1, 100] nm.")

    # --- 4. Fitting Execution ---

    def run_fit(self, method='least_squares'):
        """
        Executes sequential fitting backwards from the latest frame (highest SNR)
        to the earliest frame.
        """
        if self.data is None or self.params is None:
            print("Error: Load data and setup parameters before running fit.")
            return None

        print(f"\nStarting sequential backward fit across {len(self.file_names)} profiles...")
        x = self.data['q'].values
        n_profiles = len(self.file_names)

        # Apply fitting q-range
        mask = (x >= self.q_range[0]) & (x <= self.q_range[1])
        x_fit = x[mask]

        fitted_records = []
        fits_df = pd.DataFrame({'q': x})

        # Initialize with the robust parameters from final/late stages
        current_params = self.params.copy()

        # Iterate backwards from last frame to first
        for pos in reversed(range(n_profiles)):
            frame_name = self.file_names[pos]
            t_sec = self.time_seconds[pos]
            y_all = self.data[str(pos)].values
            dy_all = self.data['error_' + str(pos)].values

            y_fit = y_all[mask]
            dy_fit = np.maximum(dy_all[mask], 0.01)

            # Check if this frame has detectable scattering above baseline
            # Measure signal around the characteristic Guinier knee (q ~ 0.012 - 0.025 nm^-1)
            knee_mask = (x_fit >= 0.012) & (x_fit <= 0.025)
            knee_signal = np.median(y_fit[knee_mask]) if np.any(knee_mask) else 0.0
            is_induction = knee_signal < 6.0  # signal threshold for true particle formation

            if is_induction:
                # During induction period, fix Rmean and stdDev to the earliest detected particle size
                # to prevent unconstrained divergence in pure noise, and fit only phi and bkg
                fit_params = current_params.copy()
                fit_params['Rmean'].vary = False
                fit_params['stdDev'].vary = False
                out = minimize(self._residuals, fit_params, args=(x_fit, y_fit, dy_fit), method=method)
            else:
                fit_params = current_params.copy()
                fit_params['Rmean'].vary = True
                fit_params['stdDev'].vary = True
                out = minimize(self._residuals, fit_params, args=(x_fit, y_fit, dy_fit), method=method)
                # Inherit parameters for earlier frames
                current_params = out.params.copy()

            # Record fit metrics
            record = {
                'frame_index': pos,
                'filename': frame_name,
                'time_seconds': t_sec,
                'Rmean': out.params['Rmean'].value,
                'std_err_Rmean': out.params['Rmean'].stderr if out.params['Rmean'].stderr is not None else 0.0,
                'stdDev': out.params['stdDev'].value,
                'std_err_stdDev': out.params['stdDev'].stderr if out.params['stdDev'].stderr is not None else 0.0,
                'phi': out.params['phi'].value,
                'std_err_phi': out.params['phi'].stderr if out.params['phi'].stderr is not None else 0.0,
                'bkg': out.params['bkg'].value,
                'std_err_bkg': out.params['bkg'].stderr if out.params['bkg'].stderr is not None else 0.0,
                'chisqr': out.chisqr,
                'redchi': out.redchi,
                'is_induction': is_induction
            }
            fitted_records.append(record)

            # Store fitted full curve
            model_full = self._model_func(out.params, x)
            fits_df[f'fit_{frame_name}'] = model_full

            # Plot individual frame fit
            self._plot_frame_fit(x, y_all, dy_all, model_full, out.params, frame_name, pos, t_sec, out.redchi)

            if pos % 10 == 0 or pos == n_profiles - 1:
                print(f"Frame {pos:2d} ({frame_name[:28]}): <R> = {record['Rmean']:5.1f} nm, "
                      f"sigma = {record['stdDev']:4.1f} nm, phi = {record['phi']:8.2e}, redchi = {record['redchi']:5.2f}")

        # Reverse records to chronological order
        fitted_records.reverse()
        self.fit_results_df = pd.DataFrame(fitted_records)
        self.fit_results_df.to_csv(os.path.join(self.output_path, 'Polydisperse_Sphere_fit_parameters.csv'), index=False)

        self.fitted_curves_df = fits_df
        self.fitted_curves_df.to_csv(os.path.join(self.output_path, 'Fitted_curves_full_model.csv'), index=False)

        # Generate summary figures
        self._plot_parameter_trajectories()
        self._plot_representative_fits()

        print(f"\nFitting complete! All results saved in: {self.output_path}")
        return self.fit_results_df

    # --- 5. Visualization Helpers ---

    def _plot_frame_fit(self, x, y, dy, y_model, params, frame_name, pos, t_sec, redchi):
        fig, ax = plt.subplots(figsize=(6, 5))
        mask_pos = y > 0
        ax.errorbar(x[mask_pos], y[mask_pos], yerr=dy[mask_pos], fmt='o', markersize=2,
                    color='#3b528b', alpha=0.6, label='Data (ID02, 31m)')
        ax.plot(x, y_model, color='#d95f02', linewidth=2,
                label=f'Fit (<R>={params["Rmean"].value:.1f} nm, $\\sigma$={params["stdDev"].value:.1f} nm)')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(self.q_range[0] * 0.9, self.q_range[1] * 1.5)
        ax.set_ylim(bottom=1e-3, top=max(1e3, np.max(y[mask_pos]) * 2.0))
        ax.set_xlabel(r'$q$ (nm$^{-1}$)', fontsize=12)
        ax.set_ylabel(r'Intensity (cm$^{-1}$)', fontsize=12)
        ax.set_title(f"ND-i-003 Frame {pos} ({frame_name})\n$t = {t_sec:.3f}$ s, $\\chi^2_{{\\rm red}} = {redchi:.2f}$", fontsize=10)
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, which='both', linestyle=':', color='gray', alpha=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_path, 'frame_fits', f"fit_frame_{pos:03d}_{frame_name}.png"), dpi=120)
        plt.close(fig)

    def _plot_parameter_trajectories(self):
        """ 4-panel kinetic trajectory figure showing Rmean, stdDev, phi, and chi2 vs time. """
        df = self.fit_results_df
        t = df['time_seconds'].values
        mask_grow = ~df['is_induction']
        t_onset = df.loc[mask_grow, 'time_seconds'].min() if np.any(mask_grow) else 0.0

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        fig.suptitle('ND-i-003: Discrete ACC Nanosphere Growth Kinetics (10 mM, pH 11.3)', fontsize=14, fontweight="bold")

        for ax in axes.flat:
            if t_onset > 0:
                ax.axvspan(0, t_onset, color='lightgray', alpha=0.35, label='Induction Period' if ax == axes[0, 0] else None)

        # 1. Mean Radius <R>
        ax = axes[0, 0]
        ax.errorbar(t[mask_grow], df.loc[mask_grow, 'Rmean'], yerr=df.loc[mask_grow, 'std_err_Rmean'],
                    fmt='o-', color='#1f77b4', ecolor='#aec7e8', capsize=3, markersize=5, label=r'Particle Radius $\langle R \rangle$ (nm)')
        ax.set_ylabel(r'Mean Radius $\langle R \rangle$ (nm)', fontsize=12)
        ax.set_ylim(20, max(120, df.loc[mask_grow, 'Rmean'].max() * 1.25))
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(fontsize=10)

        # 2. Size Dispersion stdDev
        ax = axes[0, 1]
        ax.errorbar(t[mask_grow], df.loc[mask_grow, 'stdDev'], yerr=df.loc[mask_grow, 'std_err_stdDev'],
                    fmt='s-', color='#ff7f0e', ecolor='#ffbb78', capsize=3, markersize=5, label=r'Dispersion $\sigma_R$ (nm)')
        poly_ratio = (df.loc[mask_grow, 'stdDev'] / df.loc[mask_grow, 'Rmean']) * 100.0
        ax_twin = ax.twinx()
        ax_twin.plot(t[mask_grow], poly_ratio, '--', color='#2ca02c', alpha=0.7, label=r'Polydispersity $\sigma/\langle R \rangle$ (%)')
        ax_twin.set_ylabel(r'Polydispersity (%)', fontsize=11, color='#2ca02c')
        ax_twin.set_ylim(0, 80)
        ax.set_ylabel(r'Standard Deviation $\sigma_R$ (nm)', fontsize=12)
        ax.set_ylim(10, 60)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(fontsize=10, loc='upper left')

        # 3. Volume Fraction phi
        ax = axes[1, 0]
        ax.errorbar(t, df['phi'], yerr=df['std_err_phi'], fmt='^-', color='#2ca02c',
                    ecolor='#98df8a', capsize=3, markersize=5, label=r'Scale / Volume Fraction $\phi$')
        ax.set_xlabel(r'Reaction Time $ (s)', fontsize=12)
        ax.set_ylabel(r'Volume Fraction $\phi$', fontsize=12)
        ax.set_yscale('log')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(fontsize=10)

        # 4. Reduced Chi-Square redchi
        ax = axes[1, 1]
        ax.plot(t, df['redchi'], 'd-', color='#d62728', markersize=5, label=r'$\chi^2_{\rm red}$')
        ax.axhline(1.0, color='black', linestyle='--', alpha=0.5, label=r'Ideal $\chi^2_{\rm red} = 1$')
        ax.set_xlabel(r'Reaction Time $ (s)', fontsize=12)
        ax.set_ylabel(r'Reduced $\chi^2$', fontsize=12)
        ax.set_ylim(0, max(5.0, min(20.0, df.loc[mask_grow, 'redchi'].max() * 1.2)))
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(fontsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_path, 'kinetic_trajectories.png'), dpi=300)
        plt.savefig(os.path.join(self.output_path, 'kinetic_trajectories.pdf'))
        plt.close(fig)

    def _plot_representative_fits(self):
        """ Publication-ready multi-curve comparison plot. """
        df_res = self.fit_results_df
        x = self.data['q'].values

        # Select 5 key milestone frames: early induction, onset, rapid growth, late frame, final stage
        indices = [0, 24, 34, 49, len(df_res) - 1]
        colors = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']

        fig, ax = plt.subplots(figsize=(8, 7))

        for idx, col in zip(indices, colors):
            row = df_res.iloc[idx]
            fname = row['filename']
            t_sec = row['time_seconds']
            R = row['Rmean']
            redchi = row['redchi']

            y = self.data[str(idx)].values
            dy = self.data['error_' + str(idx)].values
            y_model = self.fitted_curves_df[f'fit_{fname}'].values

            mask = y > 0
            ax.errorbar(x[mask], y[mask], yerr=dy[mask], fmt='o', markersize=2,
                        color=col, alpha=0.4, linestyle='none')
            ax.plot(x, y_model, color=col, linewidth=2,
                    label=f"$t = {t_sec:.2f}$ s (Frame {idx+1}): $\\langle R \\rangle = {R:.1f}$ nm ($\\chi^2_{{\\rm red}}={redchi:.1f}$)")

        # Add Porod trend line (q^-4)
        q_ref = np.geomspace(0.025, 0.07, 100)
        I_ref = 3e-6 * (q_ref**-4)
        ax.plot(q_ref, I_ref, 'k--', linewidth=1.5, label=r'Porod $q^{-4}$')
        ax.text(0.045, 1.2e-1, r'$q^{-4}$', fontsize=12, fontweight='bold')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.0019, 0.12)
        ax.set_ylim(1e-2, 1e4)
        ax.set_xlabel(r'Scattering Vector $q$ (nm$^{-1}$)', fontsize=14)
        ax.set_ylabel(r'Intensity $I(q)$ (cm$^{-1}$)', fontsize=14)
        ax.set_title('ND-i-003: Discrete ACC Nanosphere Form Factor Evolution\n(10 mM $\\mathrm{CaCO}_3$, pH 11.3, ESRF ID02)', fontsize=12)
        ax.grid(True, which='both', linestyle=':', color='gray', alpha=0.4)
        ax.legend(fontsize=9, loc='lower left')

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_path, 'representative_fits_comparison.png'), dpi=300)
        plt.savefig(os.path.join(self.output_path, 'representative_fits_comparison.pdf'))
        plt.close(fig)


# --- Main Entry Point ---

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    PROJECT_NAME = '10mM_fitting'
    DATA_DIRS = [
        "longDistance/",
        "finalStages/"
    ]
    Q_RANGE = (0.002, 0.09)

    print("=" * 70)
    print("  ND-i-003: Discrete Polydisperse ACC Nanosphere Kinetic Fitting")
    print("=" * 70)

    fitter = SinglePolydisperseSphereFitter(
        project_name=PROJECT_NAME,
        data_directories=DATA_DIRS,
        q_range=Q_RANGE
    )

    fitter.load_data()
    fitter.setup_parameters()

    if fitter.data is not None:
        results_df = fitter.run_fit()
        if results_df is not None:
            print("\n--- Fitting Completed Successfully ---")
            print(f"Total profiles fitted: {len(results_df)}")
            print("\nLatest 5 profiles:")
            print(results_df[['filename', 'time_seconds', 'Rmean', 'stdDev', 'phi', 'redchi']].tail())
    else:
        print("Data failed to load.")
