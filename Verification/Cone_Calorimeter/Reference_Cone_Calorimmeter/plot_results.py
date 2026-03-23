import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import sys
import os
from matplotlib.lines import Line2D


# Determine the current script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path = os.path.join(SCRIPT_DIR, 'reference_cc.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'reference_cc_T.png')
results_file_path3 = os.path.join(SCRIPT_DIR, 'reference_cc_Concentration.png')


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
    Thermakin=pd.read_csv("ThermaKin/output_ThermaKin.csv")
    Gpyro=pd.read_csv("reference_cc_summary_01_0001.csv")
    #Gpyro2=pd.read_csv("GpyroV0.8/PVC_summary_01_0001.csv")
    FDS=pd.read_csv("FDS/reference_cc_devc.csv", skiprows=1)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

style_traces = {
    "Thermakin": {"color": "black", "linestyle": "-", "linewidth": 4},
    "GpyroV2": {"color": "blue", "linestyle": "--", "linewidth": 4},
    "GpyroV1": {"color": "gray", "linestyle": "-", "linewidth": 4},
    "FDS": {"color": "red", "linestyle": "-.", "linewidth": 4}
}




plt.plot(Thermakin["time (s)"].values,1000* Thermakin["MLR (kg/s)"].values,label="ThermaKin", **style_traces["Thermakin"])
plt.plot(Gpyro['t'].values, Gpyro['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, label= "Gpyro", **style_traces["GpyroV2"])
#plt.plot(Gpyro2['t'], Gpyro2['009_MLR( 0.0000_ 0.0000_ 0.0000)'], label= "Gpyro V0.8",**style_traces["GpyroV1"])
plt.plot(FDS['Time'].values, 1000*FDS['MLR'].values, label= 'FDS', **style_traces["FDS"])


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



plt.figure()

# Style par code (couleur uniquement)
style_traces = {
    "Thermakin": {"color": "black", "linewidth": 4},
    "Gpyro": {"color": "blue", "linewidth": 3},
    "FDS": {"color": "red", "linewidth": 2},
}

# Enlève 'label' dans depth_styles
depth_styles = {
    "z=0": {"linestyle": "-", "marker": " "},
    "z=3": {"linestyle": "--", "marker": " "},
    "z=6": {"linestyle": ":", "marker": " "},
}



# ----------- THERMAKIN -----------
plt.plot(Thermakin["time (s)"].values, Thermakin["T(z=0)"].values - 273.15,
         label="ThermaKin", **style_traces["Thermakin"], **depth_styles["z=0"])
plt.plot(Thermakin["time (s)"].values, Thermakin["T(z=3)"].values - 273.15,
         **style_traces["Thermakin"], **depth_styles["z=3"])
plt.plot(Thermakin["time (s)"].values, Thermakin["T(z=6)"].values - 273.15,
         **style_traces["Thermakin"], **depth_styles["z=6"])

# ----------- GPYRO -----------
plt.plot(Gpyro['t'].values, Gpyro['002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)'].values,
         label="Gpyro", **style_traces["Gpyro"], **depth_styles["z=0"])
plt.plot(Gpyro['t'].values, Gpyro['003_TEMPERATURE( 0.0030_ 0.0000_ 0.0000)'].values ,
         **style_traces["Gpyro"], **depth_styles["z=3"])
plt.plot(Gpyro['t'].values, Gpyro['004_TEMPERATURE( 0.0060_ 0.0000_ 0.0000)'].values,
         **style_traces["Gpyro"], **depth_styles["z=6"])

# ----------- FDS -----------
plt.plot(FDS['Time'].values, FDS['Tz0'].values ,
         label="FDS", **style_traces["FDS"], **depth_styles["z=0"])
plt.plot(FDS['Time'].values, FDS['Tz3'].values,
         **style_traces["FDS"], **depth_styles["z=3"])
plt.plot(FDS['Time'].values, FDS['Tz6'].values ,
         **style_traces["FDS"], **depth_styles["z=6"])


mpl.rcParams['legend.fontsize'] = 25  


# ----------- LEGENDES -----------

# Légende 1 : codes
legend_codes = [
    Line2D([0], [0], color=style_traces[k]["color"], lw=3, label=k)
    for k in style_traces
]

# Légende 2 : profondeur avec valeurs
depth_labels = {
    "z=0": "Depth: 0 mm",
    "z=3": "Depth: 3 mm",
    "z=6": "Depth: 6 mm"
}
legend_depths = [
    Line2D([0], [0], color="black", lw=2, label=depth_labels[k], **depth_styles[k])
    for k in depth_styles
]

# Affichage côte à côte en bas
first_legend = plt.legend(handles=legend_codes, loc='lower center', bbox_to_anchor=(0.45, 0.0), ncol=1, frameon=False)
second_legend = plt.legend(handles=legend_depths, loc='lower center', bbox_to_anchor=(0.8, 0.0), ncol=1, frameon=False)
plt.gca().add_artist(first_legend)

# ----------- AXES -----------

plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.xlim([0, 500])
plt.tight_layout()

plt.savefig(results_file_path2)


#%%
style_traces = {
    "Thermakin": {"linestyle": "-", "linewidth": 3, "marker":" ", "markersize":10},
    "Gpyro": {"linestyle": "--", "linewidth": 3}
}

# Enlève 'label' dans depth_styles
species_styles = {
    "CA": {"color": "blue"},
    "CB": {"color": "red"},
    "CC": {"color": "black"},
}
plt.figure()

# ----------- THERMAKIN -----------
plt.plot(Thermakin["time (s)"].values, Thermakin["C_A (kg/m^3)"].values,
         label="ThermaKin", **style_traces["Thermakin"], **species_styles["CA"])
plt.plot(Thermakin["time (s)"].values, Thermakin["C_B (kg/m^3)"].values,
         **style_traces["Thermakin"], **species_styles["CB"])
plt.plot(Thermakin["time (s)"].values, Thermakin["C_C (kg/m^3)"].values,
         **style_traces["Thermakin"], **species_styles["CC"])

# ----------- GPYRO -----------
plt.plot(Gpyro['t'].values, Gpyro['010_CI( 0.0030_ 0.0000_ 0.0000)'].values,
         label="Gpyro", **style_traces["Gpyro"], **species_styles["CA"])
plt.plot(Gpyro['t'].values, Gpyro['011_CI( 0.0030_ 0.0000_ 0.0000)'].values ,
         **style_traces["Gpyro"], **species_styles["CB"])
plt.plot(Gpyro['t'].values, Gpyro['012_CI( 0.0030_ 0.0000_ 0.0000)'].values,
         **style_traces["Gpyro"], **species_styles["CC"])



# ----------- LEGENDES -----------

# Légende 1 : codes
legend_codes = [
    Line2D([0], [0], color= 'black', **style_traces[k], label=k)
    for k in style_traces
]
#
# Légende 2 : profondeur avec valeurs
species_labels = {
    "CA": "MAT A",
    "CB": "MAT B",
    "CC": "MAT C"
}

legend_species = [
    Line2D([0], [0], lw=4, label=species_labels[k], **species_styles[k])
    for k in species_styles
]

# Affichage côte à côte en bas
first_legend = plt.legend(handles=legend_codes, loc='lower center', bbox_to_anchor=(0.8, 0.75), ncol=1, frameon=False)
second_legend = plt.legend(handles=legend_species, loc='lower center', bbox_to_anchor=(0.45, 0.7), ncol=1, frameon=False)
plt.gca().add_artist(first_legend)

# ----------- AXES -----------

plt.xlabel("Time (s)")
plt.ylabel("Concentration (Kg/m³)")
plt.xlim([0, 500])
plt.tight_layout()

plt.savefig(results_file_path3)


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
error_value = compute_max_absolute_error(
    FDS['Time'].values,
    1000*FDS['MLR'].values,
    Gpyro["t"].values,
    Gpyro['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values
    )

# Define a threshold for validation
threshold = 2.5  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error : {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




