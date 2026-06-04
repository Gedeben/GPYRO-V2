import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys

# Define script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_analytical = os.path.join(SCRIPT_DIR, "analytical_results.csv")


file_gpyro = os.path.join(SCRIPT_DIR, "convective_cooling_summary_01_0001.csv")
results_file_path = os.path.join(SCRIPT_DIR, 'convective_cooling.png')

# Read data
data_exp = pd.read_csv(file_analytical, delimiter=',')
data_gpyro = pd.read_csv(file_gpyro, delimiter=',')

# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 26  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

# color and style parameters for the plots 
param_sets = [
    { 'color':'blue', 'linewidth': 3, 'linestyle': '-'},
    { 'color': 'black', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4}
    ]


# Style parameters for theoretical (analytical) curves
param_th_sets = [
    {'label': 'Analytical (Front)', 'color': 'blue', 'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 1, 'markeredgewidth': 4},
    {'label': 'Analytical (Back)', 'color': 'black', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4},
    {'label': 'Analytical (4 cm)', 'color': 'red', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4}
]

# Style parameters for Gpyro simulation curves (black & white)
param_gpyro_sets = [
    {'label': 'Gpyro (Front)', 'color':'blue', 'linewidth': 3},
    {'label': 'Gpyro (Back)' , 'color':'black', 'linewidth': 3},
    {'label': 'Gpyro (50 cm)' , 'color':'red', 'linewidth': 3}
]

# Load data
try:
    data_exp = pd.read_csv(file_analytical, delimiter=',')
    data_gpyro = pd.read_csv(file_gpyro, delimiter=',')

except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)
#%%
# Create figure
plt.figure()

# Plot analytical data
plt.plot(np.array(data_exp["Time"]), np.array(data_exp["Back"]), ** param_th_sets[1])
plt.plot(np.array(data_gpyro["t"]),  np.array(data_gpyro["004_TEMPERATURE( 1.0000_ 0.0000_ 0.0000)"]), **param_gpyro_sets[1])


plt.plot(np.array(data_exp["Time"]), np.array(data_exp["50 cm"]),**param_th_sets[2])
plt.plot(np.array(data_gpyro["t"]), np.array(data_gpyro['003_TEMPERATURE( 0.5000_ 0.0000_ 0.0000)']), **param_gpyro_sets[2])


plt.plot(np.array(data_exp["Time"]), np.array(data_exp["Front"]), ** param_th_sets[0])

plt.plot(np.array(data_gpyro["t"]), np.array(data_gpyro['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)']), **param_gpyro_sets[0])


# Labels and formatting
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.legend(frameon=False)
plt.tick_params()
plt.tight_layout()

# Save figure
plt.savefig(results_file_path)

#%%



# Function to compute absolute error between simulation and analytical solution
def compute_absolute_error_uniform_grid(sim_depths, sim_temps, ana_depths, ana_temps, num_points=100):
    # Define a uniform grid over the common depth range
    min_depth = max(min(sim_depths), min(ana_depths))
    max_depth = min(max(sim_depths), max(ana_depths))
    uniform_depths = np.linspace(min_depth, max_depth, num_points)

    # Interpolate both datasets on the uniform grid
    sim_interp = np.interp(uniform_depths, sim_depths, sim_temps)
    ana_interp = np.interp(uniform_depths, ana_depths, ana_temps)

    # Compute sum of absolute differences
    error = np.mean(np.abs(sim_interp - ana_interp))
    return error

# Compute error
error_value = compute_absolute_error_uniform_grid(
    np.array(data_gpyro["t"].values),
    np.array(data_gpyro["004_TEMPERATURE( 1.0000_ 0.0000_ 0.0000)"].values),
    np.array(data_exp["Time"].values),
    np.array(data_exp["Back"].values)
)
# Define a threshold for validation
threshold = 1.0  # You can adjust this value based on your tolerance

# Print error and exit accordingly
#print(f"Validation error (sum of absolute differences on uniform grid): {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)














