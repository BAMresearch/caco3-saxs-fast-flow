# -*- coding: utf-8 -*-
# Author: Tomasz Stawski, tomasz.stawski@bam.de
"""
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import numpy as np

# Define directories (update these paths)
dir1 = "longDistance/"
dir2 = "shortDistance/"
dir3 = "finalStages/"  # Add a third directory

# Ensure columns are numeric
def read_data_from_directory_fixed(directory):
    data = []
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".dat"):
            filepath = os.path.join(directory, filename)
            # Read the .dat file assuming three columns without headers
            df = pd.read_csv(filepath, delim_whitespace=True, header=None, names=["q_nm^-1", "intensity_cm^-1", "dI"])
            # Convert columns to numeric
            df["q_nm^-1"] = pd.to_numeric(df["q_nm^-1"], errors="coerce")
            df["intensity_cm^-1"] = pd.to_numeric(df["intensity_cm^-1"], errors="coerce")
            df["dI"] = pd.to_numeric(df["dI"], errors="coerce")
            df = df.dropna()  # Clean up any invalid rows
            data.append((filename, df))
    return data

# Read data for all directories with fixed parsing
data_dir1 = read_data_from_directory_fixed(dir1)
data_dir2 = read_data_from_directory_fixed(dir2)
data_dir3 = read_data_from_directory_fixed(dir3)

# Define color palettes
palette1 = cm.get_cmap("cool", len(data_dir1))
palette2 = cm.get_cmap("autumn", len(data_dir2))
palette3 = cm.get_cmap("Greys", len(data_dir3))

# Set Seaborn style
sns.set_style("whitegrid")  # You can use "darkgrid", "white", etc.
sns.set(rc={"axes.facecolor": "lightgray"})  # Set background of graph area

# Create subplots for side-by-side graphs
fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [2, 1]})

# Main plot (dir1 and dir2)
ax = axes[0]
for i, (filename, df) in enumerate(data_dir1):
    ax.errorbar(df["q_nm^-1"], df["intensity_cm^-1"], yerr=df["dI"],
                fmt='-', color=palette1(i), alpha=0.8, capsize=2, label=dir1 if i == 0 else None)

for i, (filename, df) in enumerate(data_dir2):
    ax.errorbar(df["q_nm^-1"], df["intensity_cm^-1"], yerr=df["dI"],
                fmt='-', color=palette2(i), alpha=0.8, capsize=2, label=dir2 if i == 0 else None)

# Add a q^-4 trend line
q_min = min([df["q_nm^-1"].min() for _, df in data_dir1 + data_dir2 + data_dir3])
q_max = max([df["q_nm^-1"].max() for _, df in data_dir1 + data_dir2 + data_dir3])
q_values = np.logspace(np.log10(q_min), np.log10(q_max), 500)
trend_line = 1e-2 * q_values**-4  # Adjust scaling for visibility

ax.plot(q_values, trend_line, linestyle=':', color='darkmagenta', linewidth=2, label='$q^{-4}$')
ax.text(q_values[-1], trend_line[-1] * 1.5, "$q^{-4}$", fontsize=14, color='darkmagenta')

# Set main plot properties
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-6, 1e7)  # Fix y-axis range
ax.set_xlabel('$q$ (nm$^{-1}$)', fontsize=18)
ax.set_ylabel('Intensity (cm$^{-1}$)', fontsize=18)
ax.set_title('Main Data (dir1 and dir2)', fontsize=14)
ax.tick_params(axis='both', labelsize=14)
ax.grid(True, which="both", linestyle='--', linewidth=1, color='darkgray')
ax.legend(fontsize=12, loc='lower left', title='Series')

# Final Stages plot (dir3)
ax2 = axes[1]
for i, (filename, df) in enumerate(data_dir3):
    ax2.errorbar(df["q_nm^-1"], df["intensity_cm^-1"], yerr=df["dI"],
                 fmt='-', color=palette3(i), linewidth=2.5, alpha=0.8, capsize=2)

# Add q^-4 trend line to final stages plot
ax2.plot(q_values, trend_line, linestyle=':', color='darkmagenta', linewidth=2, label='$q^{-4}$')

# Set final stages plot properties
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(1e-6, 1e7)  # Fix y-axis range for consistency
ax2.set_xlabel('$q$ (nm$^{-1}$)', fontsize=18)
ax2.set_title('Final Stages (dir3)', fontsize=14)
ax2.tick_params(axis='both', labelsize=14)
ax2.grid(True, which="both", linestyle='--', linewidth=1, color='darkgray')

# Adjust layout and show the plots
plt.tight_layout()
# Save the figure as PNG at 600 dpi
plt.savefig("mygraph.png", dpi=600)
plt.show()
