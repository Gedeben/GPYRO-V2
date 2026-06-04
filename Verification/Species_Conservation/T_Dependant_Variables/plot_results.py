import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import numpy as np


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_path = os.path.join(SCRIPT_DIR, 'T_Dependant_Variables_summary_01_0001.csv')


compare_to_old_gpyro = False

# Execute only if the specific argument is provided
if len(sys.argv) > 1 and sys.argv[1] == "compare_old_gpyro":
    file_path_old = os.path.join(SCRIPT_DIR, 'Gpyro_V08_T_Dependant_Variables_summary_01_0001.csv')
    data_old = pd.read_csv(file_path_old)
    compare_to_old_gpyro= True

# analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')
# results_file_path = os.path.join(SCRIPT_DIR, 'Controlled_K.png')



#%% color and style parameters for the plots in black and white
param_sets = [
    {'color': 'blue' , 'linewidth': 3, 'linestyle': '-'},
   {'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 10,  'markeredgewidth': 4},
    {'color': 'red' , 'linewidth': 3, 'linestyle': '-'},

    ]

#%%


# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 30
mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['legend.fontsize'] = 30
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['figure.figsize'] = (20, 12)


# Load data
try:
    data = pd.read_csv(file_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

# Time vector
t = data['t'].values


T = 300 + t

C_N = 1550 * (T/300)**0.3

K_N = 0.2 * (T/300)**0.6 + (5.6e-8*0.002*T**3)

RHO_BULK = 1430 * (T/300)**0.25

# Dictionary mapping: {column_name : (plot_title, y_label)}
col_labels = {
    "001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"      : ("Temperature", "T [°C]"),
    "002_THERMAL_CONDUCTIVITY_Z( 0.0000_ 0.0000_ 0.0000)" : ("Thermal conductivity", "k [W/m·K]"),
    "003_SPECIFIC_HEAT_CAPACITY( 0.0000_ 0.0000_ 0.0000)" : ("Specific heat capacity", "Cp [J/kg·K]"),
    "005_BULK_DENSITY( 0.0000_ 0.0000_ 0.0000)"     : ("Bulk density", "ρ [kg/m³]"),
}

# Dictionnaire de solutions analytiques
analytical_solutions = {
    "001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)": T - 273.15,
    "002_THERMAL_CONDUCTIVITY_Z( 0.0000_ 0.0000_ 0.0000)": K_N,
    "003_SPECIFIC_HEAT_CAPACITY( 0.0000_ 0.0000_ 0.0000)": C_N,
    "005_BULK_DENSITY( 0.0000_ 0.0000_ 0.0000)": RHO_BULK,
}



cols_to_plot = list(col_labels.keys())


# --------------------
# Combined subplot figure (3x3 grid)
# --------------------
fig, axes = plt.subplots(2, 2, sharex=True)
axes = axes.flatten()

for i, col in enumerate(cols_to_plot):
    title, ylabel = col_labels[col]
    axes[i].plot(t, data[col], **param_sets[0], label="Gpyro")
    if compare_to_old_gpyro:
        axes[i].plot(t, data_old[col], **param_sets[2], label="Gpyro V0.8")
    if col in analytical_solutions:
        axes[i].plot(t, analytical_solutions[col], **param_sets[1], label='Analytical')
    axes[i].set_title(title)
    axes[i].set_xlabel("Time [s]")
    axes[i].set_ylabel(ylabel)
    axes[i].legend(frameon=False)
    axes[i].ticklabel_format(style="plain", axis="both") 

# Hide unused subplots if any
for j in range(len(cols_to_plot), len(axes)):
    fig.delaxes(axes[j])

fig.tight_layout()
combined_file = Path(SCRIPT_DIR) / "All_Variables_Subplots_TDep.png"
plt.savefig(combined_file, dpi=300)
# plt.close()


#%%


# Function to compute absolute error between simulation and analytical solution
def compute_case_error(t, data, analytical_solutions):
    """
    Compute the maximum mean absolute error among all variables 
    between Gpyro results and analytical solutions.
    """
    errors = {}
    for col, ana_values in analytical_solutions.items():
        sim_values = data[col].values
        error = np.mean(np.abs(sim_values - ana_values))
        errors[col] = error
    
    # Take the worst-case error
    max_error = max(errors.values())
    return max_error, errors


error_value, all_errors = compute_case_error(t, data, analytical_solutions)

# Print details
for col, err in all_errors.items():
    print(f"{col_labels[col][0]}: mean absolute error = {err:.6e}")

print(f"\nOverall validation error (worst among variables): {error_value:.6e}")

# Validation threshold
threshold = 0.75
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)


