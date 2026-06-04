import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_path = os.path.join(SCRIPT_DIR, 'simple_TGA_summary_01_0001.csv')
analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')
results_file_path = os.path.join(SCRIPT_DIR, 'TGA_simple.png')



#%% color and style parameters for the plots in black and white
param_sets = [
    { 'color':'blue', 'linewidth': 3, 'linestyle': '-'},
    { 'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 5,  'markeredgewidth': 4} ]


#%%



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
    data = pd.read_csv(file_path)
    analytical_data = pd.read_csv(analytical_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)
    data = pd.read_csv(file_path)


MLR_an= analytical_data["MLR = -dm/dt (g/s)"].values
t_an= analytical_data['Time (s)'].values

MLR= data["002_MLR( 0.0000_ 0.0000_ 0.0000)"].values
t= data['t'].values

Temperature= data["001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"]


# Plot the temperature profile
plt.plot(t, MLR, **param_sets[0], label = 'Gpyro')
plt.plot(t_an, MLR_an, **param_sets[1], label ='Analytical')

plt.xlabel('Time [s]')
plt.ylabel("MLR [g/m²/s]")
plt.xlim([0,5000])
plt.legend(frameon=False)
plt.tick_params()
plt.tight_layout()
plt.savefig(results_file_path)

#%%



# Function to compute absolute error between simulation and analytical solution
def compute_absolute_error_uniform_grid(sim_depths, sim_temps, ana_depths, ana_temps, num_points=1000):
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
    t, MLR,t_an, MLR_an)  

# Define a threshold for validation
threshold = 0.001  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error (mean of absolute differences ): {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




