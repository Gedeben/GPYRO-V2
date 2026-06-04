import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_path = os.path.join(SCRIPT_DIR, 'insulated_steel_plate_profile_01_01_TEMPERATURE.csv')
analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')
results_file_path = os.path.join(SCRIPT_DIR, 'insulated_steel_plate.png')



#%% color and style parameters for the plots 
param_sets = [
    { 'color':'blue', 'linewidth': 3, 'linestyle': '-'},
    { 'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4}
    ]

#%%


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 30  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)




#%%

# Load data
try:
    data = pd.read_csv(file_path, header=None)
    analytical_data = pd.read_csv(analytical_path, delimiter=',', decimal='.')
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

# Extract temperature and depth
Temperature_equilibrium = data.iloc[-1, 1:].values
depths = data.iloc[0, 1:].astype(float).values




x_analytical = np.array(analytical_data['depth'])
y_analytical = np.array(analytical_data['temp'])

# Plot the temperature profile
plt.plot(x_analytical, y_analytical, label='Analytical', **param_sets[1])
plt.plot(depths, Temperature_equilibrium, label='Gpyro', **param_sets[0])

plt.xlabel('Depth (m)')
plt.ylabel("Temperature [°C]")
plt.legend(frameon=False)
plt.tick_params()
plt.tight_layout()
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
    depths,
    Temperature_equilibrium,
    analytical_data['depth'].values,
    analytical_data['temp'].values
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




