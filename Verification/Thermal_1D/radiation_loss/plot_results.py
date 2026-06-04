import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import sys
import os
from matplotlib.lines import Line2D


# Determine the current script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path_1 = os.path.join(SCRIPT_DIR, 'radiation_loss_tg0k.png')
results_file_path_2 = os.path.join(SCRIPT_DIR, 'radiation_loss_tg0C.png')
# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 25  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30 
 
mpl.rcParams['figure.figsize'] = (12, 8)

# Load data
try:
    Gpyro1=pd.read_csv("Tgas0K/radiation_loss_summary_01_0001.csv")
    Gpyro2=pd.read_csv("Tgas273K/radiation_loss_summary_01_0001.csv")
    FDS1=pd.read_csv("Tgas0K/FDS/radiation_loss_devc.csv", skiprows=1)
    FDS2=pd.read_csv("Tgas273K/FDS/radiation_loss_devc.csv", skiprows=1)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)



style_traces = {
    "Gpyro": {"linewidth": 3},
    "FDS": {'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 30, 'markeredgewidth': 4},
}

depth_styles = {
    "z=0": {'color': 'blue'},
    "z=40": {'color': 'red'},
}



plt.figure()


 # ----------- GPYRO -----------
plt.plot(Gpyro1['t'].values, Gpyro1['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)'].values,
         label="Gpyro", **style_traces["Gpyro"], **depth_styles["z=0"])
plt.plot(Gpyro1['t'].values, Gpyro1['005_TEMPERATURE( 0.0400_ 0.0000_ 0.0000)'].values,
         **style_traces["Gpyro"], **depth_styles["z=40"])

# ----------- FDS -----------
plt.plot(FDS1['Time'].values, FDS1['Front'].values ,
         label="FDS", **style_traces["FDS"], **depth_styles["z=0"])
plt.plot(FDS1['Time'].values, FDS1['Back'].values ,
         **style_traces["FDS"], **depth_styles["z=40"])


# ----------- LEGENDES -----------

# Légende 1 : codes
legend_codes = [
    Line2D([0], [0], **style_traces[k], label=k, color='gray')
    for k in style_traces]


# Légende 2 : profondeur avec valeurs
depth_labels = {
    "z=0": "Depth: 0 mm",
    "z=40": "Depth: 40 mm"
}
legend_depths = [
    Line2D([0], [0],  lw=3, label=depth_labels[k], **depth_styles[k], marker='s', markersize=16, linestyle='none')
    for k in depth_styles
]


# Affichage côte à côte en bas
first_legend = plt.legend(handles=legend_codes, loc='lower center', bbox_to_anchor=(0.45, 0.7), ncol=1, frameon=False)
second_legend = plt.legend(handles=legend_depths, loc='lower center', bbox_to_anchor=(0.8, 0.7), ncol=1, frameon=False)
plt.gca().add_artist(first_legend)


# ----------- AXES -----------

plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.tight_layout()

plt.savefig(results_file_path_1)

#%%

plt.figure()

#---------- GPYRO -----------
plt.plot(Gpyro2['t'].values, Gpyro2['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)'].values,
         label="Gpyro", **style_traces["Gpyro"], **depth_styles["z=0"])
plt.plot(Gpyro2['t'].values, Gpyro2['005_TEMPERATURE( 0.0400_ 0.0000_ 0.0000)'].values,
         **style_traces["Gpyro"], **depth_styles["z=40"])

# ----------- FDS -----------
plt.plot(FDS2['Time'].values, FDS2['Front'].values ,
         label="FDS", **style_traces["FDS"], **depth_styles["z=0"])
plt.plot(FDS2['Time'].values, FDS2['Back'].values ,
         **style_traces["FDS"], **depth_styles["z=40"])



# ----------- LEGENDES -----------

# Légende 1 : codes
legend_codes = [
    Line2D([0], [0], **style_traces[k], label=k, color='gray')
    for k in style_traces]


# Légende 2 : profondeur avec valeurs
depth_labels = {
    "z=0": "Depth: 0 mm",
    "z=40": "Depth: 40 mm"
}
legend_depths = [
    Line2D([0], [0],  lw=3, label=depth_labels[k], **depth_styles[k], marker='s', markersize=16, linestyle='none')
    for k in depth_styles
]


# Affichage côte à côte en bas
first_legend = plt.legend(handles=legend_codes, loc='lower center', bbox_to_anchor=(0.45, 0.7), ncol=1, frameon=False)
second_legend = plt.legend(handles=legend_depths, loc='lower center', bbox_to_anchor=(0.8, 0.7), ncol=1, frameon=False)
plt.gca().add_artist(first_legend)

# ----------- AXES -----------

plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.tight_layout()

plt.savefig(results_file_path_2)

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
error_value_1_1 = compute_mean_absolute_error(
    FDS1['Time'].values,
    FDS1['Front'].values,
    Gpyro1["t"].values,
    Gpyro1['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)'].values
    )

error_value_1_2 = compute_mean_absolute_error(
    FDS1['Time'].values,
    FDS1['Back'].values,
    Gpyro1["t"].values,
    Gpyro1['005_TEMPERATURE( 0.0400_ 0.0000_ 0.0000)'].values
    )

error_value_2_1 = compute_mean_absolute_error(
    FDS2['Time'].values,
    FDS2['Front'].values,
    Gpyro2["t"].values,
    Gpyro2['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)'].values
    )

error_value_2_2 = compute_mean_absolute_error(
    FDS2['Time'].values,
    FDS2['Back'].values,
    Gpyro2["t"].values,
    Gpyro2['005_TEMPERATURE( 0.0400_ 0.0000_ 0.0000)'].values
    )

error_value = np.max([error_value_1_1, error_value_1_2, error_value_2_1, error_value_2_2])

# Define a threshold for validation
threshold = 0.5 # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error Tg = 0 K : {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
