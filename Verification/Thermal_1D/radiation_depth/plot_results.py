import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import sys
import os

# Path setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

Gpyro_path1 = os.path.join(SCRIPT_DIR, 'rad_profile_01_01_S.csv')
Gpyro_path2 = os.path.join(SCRIPT_DIR, 'rad_profile_02_01_S.csv')
Gpyro_path3 = os.path.join(SCRIPT_DIR, 'rad_profile_03_01_S.csv')

results_file_path1= os.path.join(SCRIPT_DIR, 'top_flux_results.png')
results_file_path2= os.path.join(SCRIPT_DIR, 'bottom_flux_results.png')
results_file_path3= os.path.join(SCRIPT_DIR, 'double_flux_results.png')

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
    Gpyro_flux_top=pd.read_csv(Gpyro_path1, header=None)
    Gpyro_flux_bottom=pd.read_csv(Gpyro_path2, header=None)
    Gpyro_flux_both=pd.read_csv(Gpyro_path3, header=None)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

# Define styles
style_traces = {
    "Gpyro": {"color": "blue", "linewidth": 3},
    "Gpyro V0.8": {"color": "green", "linewidth": 3, 'linestyle' :'-', 'marker': 's', 'markersize': 10, 'markevery': 7,  'markeredgewidth': 2},
    "Analytical": {'color': 'blue', 'marker': '+', 'linestyle': 'none', 'markersize': 16, 'markevery': 4,  'markeredgewidth': 4},
}

compare_to_old_gpyro = False
Gpyro_label="Gpyro"

# Execute only if the specific argument is provided
if len(sys.argv) > 1 and sys.argv[1] == "compare_old_gpyro":
    Gpyro_old_path1 = os.path.join(SCRIPT_DIR, 'Gpyro_V08_rad_profile_01_01_S.csv')
    Gpyro_old_path2 = os.path.join(SCRIPT_DIR, 'Gpyro_V08_rad_profile_02_01_S.csv')
    Gpyro_old_path3 = os.path.join(SCRIPT_DIR, 'Gpyro_V08_rad_profile_03_01_S.csv')
    Gpyro_flux_old_top = pd.read_csv(Gpyro_old_path1, header=None)
    Gpyro_flux_old_bottom = pd.read_csv(Gpyro_old_path2, header=None)
    Gpyro_flux_old_both = pd.read_csv(Gpyro_old_path3, header=None)
    compare_to_old_gpyro= True
    Gpyro_label="Gpyro latest"
    

#%% Top flux

K = 24      # Absorption coefficient [1/m]
Qe = 5000   # Incident radiative heat flux [W/m^2]
L = 0.10
depths = Gpyro_flux_top.iloc[0, 1:].astype(float).values

# Calculate the thickness of each layer (dz)
dz = depths[1:] - depths[:-1]
x_start_of_layer = depths[:-1]
analytical_flux_top = Qe * K * np.exp(-K * x_start_of_layer)


plt.figure()

flux_gpyro_top =Gpyro_flux_top.iloc[-1, 1:-1].values

if compare_to_old_gpyro:    
    depths_old = Gpyro_flux_old_top.iloc[0, 1:].astype(float).values
    flux_gpyro_old =Gpyro_flux_old_top.iloc[-1, 1:-1].values
    plt.plot(depths_old[1:], flux_gpyro_old, label="Gpyro V0.8", **style_traces["Gpyro V0.8"])
    
plt.plot(depths[1:], flux_gpyro_top, label=Gpyro_label, **style_traces["Gpyro"])
plt.plot(depths[1:],analytical_flux_top , label="Analytical", **style_traces["Analytical"])


plt.legend(frameon=False)
plt.xlabel("Depth (m)")
plt.ylabel("Flux (W/m³)")
plt.tight_layout()
plt.savefig(results_file_path1)

#%%Bottom flux

x_start_of_layer = depths[:-1]
analytical_flux_bottom = Qe * K * np.exp(-K * (L - x_start_of_layer))

plt.figure()
flux_gpyro_bottom =Gpyro_flux_bottom.iloc[-1, 1:-1].values

if compare_to_old_gpyro:
    depths_old = Gpyro_flux_old_bottom.iloc[0, 1:].astype(float).values
    flux_gpyro_old =Gpyro_flux_old_bottom.iloc[-1, 1:-1].values
    plt.plot(depths_old[1:], flux_gpyro_old, label="Gpyro V0.8", **style_traces["Gpyro V0.8"])
    
plt.plot(depths[1:], flux_gpyro_bottom, label=Gpyro_label, **style_traces["Gpyro"])
plt.plot(depths[1:],analytical_flux_bottom , label="Analytical", **style_traces["Analytical"])


plt.legend(frameon=False)
plt.xlabel("Depth (m)")
plt.ylabel("Flux (W/m³)")
plt.tight_layout()
plt.savefig(results_file_path2)

#%% both face flux
Qe_top = 5000   # Incident radiative heat flux from the top [W/m^2]
Qe_bottom = 5000 # Incident radiative heat flux from the bottom [W/m^2]
dz = depths[1:] - depths[:-1]
x_start_of_layer = depths[:-1]

flux_from_top = Qe_top * K * np.exp(-K * x_start_of_layer)
flux_from_bottom = Qe_bottom * K * np.exp(-K * (L - x_start_of_layer))

analytical_flux_both = (flux_from_top + flux_from_bottom)

plt.figure()
flux_gpyro_both =Gpyro_flux_both.iloc[-1, 1:-1].values

if compare_to_old_gpyro:
    depths_old = Gpyro_flux_old_both.iloc[0, 1:].astype(float).values
    flux_gpyro_old =Gpyro_flux_old_both.iloc[-1, 1:-1].values
    plt.plot(depths_old[1:], flux_gpyro_old, label="Gpyro V0.8", **style_traces["Gpyro V0.8"])
    
plt.plot(depths[1:], flux_gpyro_both, label=Gpyro_label, **style_traces["Gpyro"])
plt.plot(depths[1:],analytical_flux_both , label="Analytical", **style_traces["Analytical"])


plt.legend(frameon=False)
plt.xlabel("Depth (m)")
plt.ylabel("Flux (W/m³)")
plt.tight_layout()
plt.savefig(results_file_path3)

#%%
# -------------------------------------------------------------------
# Compute error between analytical and simulation for each curve
# -------------------------------------------------------------------

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

error_value1 = compute_mean_absolute_error(
    depths[1:], analytical_flux_top,
    depths[1:], flux_gpyro_top)

error_value2 = compute_mean_absolute_error(
    depths[1:], analytical_flux_bottom,
    depths[1:], flux_gpyro_bottom)

error_value3 = compute_mean_absolute_error(
    depths[1:], analytical_flux_both,
    depths[1:], flux_gpyro_both)


error_value=np.max([error_value1,error_value2,error_value2])
# Define a threshold for validation
threshold = 10 # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error : {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
