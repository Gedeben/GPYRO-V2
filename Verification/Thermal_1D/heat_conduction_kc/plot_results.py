#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys

# Determine the current script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# Input files
file_analytical=os.path.join(SCRIPT_DIR, "theoretical_results.csv")
file_gpyro=os.path.join(SCRIPT_DIR, "heat_conduction_kc_summary_01_0001.csv")
results_file_path = os.path.join(SCRIPT_DIR, 'heat_conduction_kc_results.png')


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 26  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

# Load data
try:
    data_th = pd.read_csv(file_analytical, delimiter=',')
    data_gpyro = pd.read_csv(file_gpyro, delimiter=',')
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

# Style parameters for theoretical (analytical) curves
param_th_sets = [
    {'label': 'HEATING (Front)', 'color': 'blue', 'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 1, 'markeredgewidth': 4},
    {'label': 'HEATING (Back)', 'color': 'red', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4},
]

# Style parameters for Gpyro simulation curves (black & white)
param_gpyro_sets = [
    {'label': 'Gpyro (Front)', 'color':'blue', 'linewidth': 3},
    {'label': 'Gpyro (Back)' , 'color':'red', 'linewidth': 3},
]

# Plot theoretical (analytical) data
plt.plot(data_th["Time"].values,
         data_th["cart_surf"].values,
         **param_th_sets[0])

plt.plot(data_gpyro["t"].values,
         data_gpyro["002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"].values,
         **param_gpyro_sets[0])


plt.plot(data_th["Time"].values,
         data_th["cart_back"].values,
         **param_th_sets[1])

plt.plot(data_gpyro["t"].values,
         data_gpyro["004_TEMPERATURE( 0.0100_ 0.0000_ 0.0000)"].values,
         **param_gpyro_sets[1])

# Finalize plot
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.legend(bbox_to_anchor=(0.43, 0.3), loc='lower left', frameon=False)
plt.tick_params()
plt.tight_layout()
plt.savefig(results_file_path)
# plt.show()
#%%

# -------------------------------------------------------------------
# Compute error between analytical and simulation for each curve
# -------------------------------------------------------------------

def compute_mean_absolute_error(sim_times, sim_values, ref_times, ref_values, num_points=200):
    """
    Interpolates both series onto a common uniform time grid and computes
    the maximum absolute error.
    """
    # Ensure overlapping time domain
    t_min = max(min(sim_times), min(ref_times))
    t_max = min(max(sim_times), max(ref_times))
    common_time = np.linspace(t_min, t_max, num_points)

    sim_interp = np.interp(common_time, sim_times, sim_values)
    ref_interp = np.interp(common_time, ref_times, ref_values)

    return np.mean(np.abs(sim_interp - ref_interp))

# Extract relevant series
t_ana = data_th["Time"].values

errors = []
# Front
err_front = compute_mean_absolute_error(
    data_gpyro["t"].values,
    data_gpyro["002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"].values,
    t_ana,
    data_th["cart_surf"].values
)
errors.append(err_front)

# Back
err_back = compute_mean_absolute_error(
    data_gpyro["t"].values,
    data_gpyro["004_TEMPERATURE( 0.0100_ 0.0000_ 0.0000)"].values,
    t_ana,
    data_th["cart_back"].values
)
errors.append(err_back)

# Final validation decision
max_error = max(errors)
threshold = 2.0  # Adjust this value as needed

print(f"Mean absolute error over curves: {max_error:.3f} °C")

if max_error <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
