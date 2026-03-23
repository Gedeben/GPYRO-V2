import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import sys
import os
from matplotlib.lines import Line2D


# Determine the current script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path = os.path.join(SCRIPT_DIR, 'mlr_2d_3d.png')



# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 30  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

# Load data
try:
    Gpyro1D=pd.read_csv("1D/reference_cc_summary_01_0001.csv")
    Gpyro2D=pd.read_csv("2D/reference_cc_summary_01_0001.csv")
    Gpyro3D=pd.read_csv("3D/reference_cc_summary_01_0001.csv")

except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

style_traces = {
    "Gpyro1D":  {"color": "black", "linestyle": "-", "marker": '+', "markevery":9, "markersize":20, "linewidth": 4},
    "Gpyro2D": {"color": "blue", "linestyle": "--", "linewidth": 4},
    "Gpyro3D": {"color": "red", "linestyle": "--", "linewidth": 3.5},
    "FDS": {"color": "red", "linestyle": "-.", "linewidth": 4}
}




plt.plot(Gpyro1D['t'].values, Gpyro1D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values, label= "Gpyro 1D", **style_traces["Gpyro1D"])
plt.plot(Gpyro2D['t'].values, 10*Gpyro2D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values, label= "Gpyro 2D", **style_traces["Gpyro2D"])
plt.plot(Gpyro3D['t'].values, 100*Gpyro3D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values, label= "Gpyro 3D", **style_traces["Gpyro3D"])




# Mise en forme
plt.xlabel("Time (s)")

plt.ylabel(r"MLR (g.m$^{-2}$.s$^{-1}$)")

plt.xlim([0, 500])
#plt.ylim([0, 25])
plt.grid(False) #, linestyle='--', linewidth=0.5, alpha=0.7)  # Grille en pointillé
plt.legend(frameon=False)  # Légende sans cadre
plt.tight_layout()  # Ajuste automatiquement la disposition


plt.savefig(results_file_path)






#%%

# -------------------------------------------------------------------
# Compute error between analytical and simulation for each curve
# -------------------------------------------------------------------

def compute_max_absolute_error(sim_times, sim_values, ref_times, ref_values, num_points=500):
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

    return np.max(np.abs(sim_interp - ref_interp))



# Compute error
error_value2D = compute_max_absolute_error(
    Gpyro1D["t"].values,
    Gpyro1D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values,
    Gpyro2D["t"].values,
    10*Gpyro2D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values
    )

error_value3D = compute_max_absolute_error(
    Gpyro1D["t"].values,
    Gpyro1D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values,
    Gpyro3D["t"].values,
    100*Gpyro3D['001_MLR( 0.0000_ 0.0500_ 0.0500)'].values
    )

error_value= max(error_value2D, error_value3D)

# Define a threshold for validation
threshold = 0.005  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error : {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




