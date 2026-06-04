import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
# analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')

file_path = os.path.join(SCRIPT_DIR, 'endothermic_reaction_summary_01_0001.csv')

results_file_path = os.path.join(SCRIPT_DIR, 'Reaction_Enthalpy_Endo.png')

#%% color and style parameters for the plots in black and white

style_traces = {
    "Gpyro": {"color": "blue", "linewidth": 3},
    "Gpyro V0.8": {"color": "green", "linewidth": 3, 'linestyle' :'-', 'marker': 's', 'markersize': 10, 'markevery': 7,  'markeredgewidth': 2},
    "Analytical": {'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 4,  'markeredgewidth': 4},
}
# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 30
mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['legend.fontsize'] = 30
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['figure.figsize'] = (12, 7)


#%%

# Load data
try:
    data = pd.read_csv(file_path)
    # analytical_data = pd.read_csv(analytical_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)
    data = pd.read_csv(file_path)

RK = data["004_REACTION_RATE_K( 0.0010_ 0.0000_ 0.0000)"].values
S = data["005_QSC( 0.0010_ 0.0000_ 0.0000)"].values
DH=-1000000
Analytical_S = RK * DH


Time = data["t"].values

plt.figure()
plt.plot(Time, S, label='Gpyro', **style_traces["Gpyro"])
plt.plot(Time, Analytical_S, label='Analytical', **style_traces["Analytical"])
plt.xlim([0,2000])
plt.xlabel('Time [s]')
plt.ylabel(r'Source Term [W.m$^{-3}$]')
#plt.title('Source Term vs Time')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(results_file_path)

#%%



# Function to compute absolute error between simulation and analytical solution
def compute_mean_absolute_error(sim_times, sim_values, ref_times, ref_values, num_points=1000):
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



# Compute error
error_value = compute_mean_absolute_error(
    Time, S,
    Time, Analytical_S
)
# Define a threshold for validation
threshold = 0.05 
# Print error and exit accordingly
print(f"Validation error (mean of absolute differences ): {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




