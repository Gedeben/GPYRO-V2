import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
# analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')

file_path = os.path.join(SCRIPT_DIR, 'sensible_enthalpy_summary_01_0001.csv')

results_file_path = os.path.join(SCRIPT_DIR, 'Reaction_Enthalpy_sensible.png')


#%% color and style parameters for the plots in black and white
style_traces = {
    "Gpyro": {"color": "blue", "linewidth": 3},
    "Gpyro V0.8": {"color": "green", "linewidth": 3, 'linestyle' :'-', 'marker': 's', 'markersize': 10, 'markevery': 7,  'markeredgewidth': 2},
    "Analytical": {'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 1,  'markeredgewidth': 4},
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

# --- Parameters from Gpyro input file ---
# Reaction
delta_H = 0  # J/kg, from DHS(1)

# Species Properties
rho_A = 1000.0  # kg/m³, Density of reactant 'PVC' (R0(1))
rho_B = 400.0   # kg/m³, Density of product 'CHAR' (R0(2))
cp_A = 1200.0   # J/kg-K, Specific heat of 'PVC' (C0(1))
cp_B = 1700.0   # J/kg-K, Specific heat of 'CHAR' (C0(2))

# reference enthalpy
T_ref = 200.0   # K, adum temperature 

# --- Calculations ---
# Solid Fraction (SF) for CHI=1
theta =(rho_B / rho_A)

# Reaction rate and temperature from simulation data
RK = data["004_REACTION_RATE_K( 0.0010_ 0.0000_ 0.0000)"].values
T = data["001_TEMPERATURE( 0.0010_ 0.0000_ 0.0000)"].values + 273.15
Yi = data["003_YI( 0.0010_ 0.0000_ 0.0000)"].values
QSC = data["005_QSC( 0.0010_ 0.0000_ 0.0000)"].values




sensible_A=1*RK*(cp_A*(T-T_ref))
sensible_B= theta*RK*(cp_B*(T-T_ref))
# Total analytical source term
dH=-sensible_A+sensible_B



Time = data["t"].values

plt.figure()
plt.plot(Time, QSC, label='Gpyro', **style_traces["Gpyro"])
plt.plot(Time, dH, label='Analytical', **style_traces["Analytical"])
plt.xlim([0,1000])
plt.xlabel('Time (s)')
plt.ylabel('Source Term (W/m³)')
#plt.title('Source Term vs Time')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(results_file_path)

#%%



# Function to compute absolute error between simulation and analytical solution
def compute_mean_absolute_error(sim_times, sim_values, ref_times, ref_values, num_points=5000):
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
m=np.max(abs(dH))
error_value = compute_mean_absolute_error(
    Time, dH,
    Time, QSC
)
error_value=error_value/m
# Define a threshold for validation
threshold = 0.001  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error (relative mean of absolute differences ): {error_value:.6f} %")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




