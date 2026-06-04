import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.lines import Line2D


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
analytical_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')

file_path_1 = os.path.join(SCRIPT_DIR, 'Oxydation_summary_01_0001.csv')
file_path_2 = os.path.join(SCRIPT_DIR, 'Oxydation_summary_02_0002.csv')
file_path_3 = os.path.join(SCRIPT_DIR, 'Oxydation_summary_03_0003.csv') 
file_path_4 = os.path.join(SCRIPT_DIR, 'Oxydation_summary_04_0004.csv')


results_file_path = os.path.join(SCRIPT_DIR, 'Oxydation_MLR.png')

#%%

style_traces = {
    "Gpyro": {"linewidth": 3},
    "Analytical": { 'marker': '+', 'linestyle': 'none', 'markersize': 12, 'markevery': 10,  'markeredgewidth': 3},
}

case_styles = [
 {'color': 'blue'},
{'color': 'black'},
{'color': 'red'},
{'color': 'green'},
]

# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 25  
mpl.rcParams['axes.labelsize'] = 25  
mpl.rcParams['legend.fontsize'] = 25  
mpl.rcParams['xtick.labelsize'] = 25 
mpl.rcParams['ytick.labelsize'] = 25  
mpl.rcParams['figure.figsize'] = (12, 8)


#%%

# Load data
try:
    data_1 = pd.read_csv(file_path_1)
    data_2 = pd.read_csv(file_path_2)
    data_3 = pd.read_csv(file_path_3)
    data_4 = pd.read_csv(file_path_4)
    analytical_data = pd.read_csv(analytical_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)
    data_1 = pd.read_csv(file_path_1)
    data_2 = pd.read_csv(file_path_2)
    data_3 = pd.read_csv(file_path_3)
    data_4 = pd.read_csv(file_path_4)

MLR_0 = data_1["002_MLR( 0.0000_ 0.0000_ 0.0000)"].values
MLR_43 = data_2["002_MLR( 0.0000_ 0.0000_ 0.0000)"].values
MLR_82 = data_3["002_MLR( 0.0000_ 0.0000_ 0.0000)"].values
MLR_205 = data_4["002_MLR( 0.0000_ 0.0000_ 0.0000)"].values

Temperature_an= analytical_data["Temperature (K)"].values - 273.15
Temperature_0= data_1["001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"]
Temperature_43= data_2["001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"]
Temperature_82= data_3["001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"]
Temperature_205= data_4["001_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"]

plt.figure()
plt.plot(Temperature_0, MLR_0, **style_traces["Gpyro"],**case_styles[0])
plt.plot(Temperature_43, MLR_43, **style_traces["Gpyro"],**case_styles[1])
plt.plot(Temperature_82, MLR_82,  **style_traces["Gpyro"],**case_styles[2])
plt.plot(Temperature_205, MLR_205,  **style_traces["Gpyro"],**case_styles[3])

plt.plot(Temperature_an, analytical_data["MLR 0.0% O2 (g/s)"],**style_traces["Analytical"],**case_styles[0])
plt.plot(Temperature_an, analytical_data["MLR 4.3% O2 (g/s)"], **style_traces["Analytical"],**case_styles[1])
plt.plot(Temperature_an, analytical_data["MLR 8.2% O2 (g/s)"],**style_traces["Analytical"],**case_styles[2])
plt.plot(Temperature_an, analytical_data["MLR 20.5% O2 (g/s)"], **style_traces["Analytical"],**case_styles[3])

plt.xlabel('Temperature (°C)')
plt.ylabel('Mass Loss Rate (g/s)')
#plt.title('Mass Loss Rate vs Temperature')


# ================== LEGENDES ==================
legend_codes = [
    Line2D([0], [0], **style_traces["Gpyro"], label="Gpyro", color="gray"),
    Line2D([0], [0], **style_traces["Analytical"], label="Analytical ", color="gray"),
]
label='Simulation MLR 0.0%',
label='Simulation MLR 4.3%',
label='Simulation MLR 8.2%',
label='Simulation MLR 20.5%',
case_labels = [
    r"$Y_{O2}$=0.0%",
    r"$Y_{O2}$=4.3%",
    r"$Y_{O2}$=8.2%",
    r"$Y_{O2}$=20.5%",
]

legend_cases = [
    Line2D(
        [0], [0],
        lw=3,
        label=case_labels[k],
        **case_styles[k],
        marker="s",
        markersize=12,
        linestyle="none"
    )
    for k in range(len(case_styles))
]

# Affichage côte à côte en bas
first_legend = plt.legend(handles=legend_codes, loc='lower center', bbox_to_anchor=(0.2, 0.7),  ncol=1, frameon=False)
second_legend = plt.legend(handles=legend_cases, loc='lower center', bbox_to_anchor=(0.2, 0.3),  ncol=1, frameon=False)
plt.gca().add_artist(first_legend)
plt.gca().add_artist(second_legend)

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
error_value_1 = compute_mean_absolute_error(
    Temperature_0,
    MLR_0,
    Temperature_an,
    analytical_data["MLR 0.0% O2 (g/s)"].values
)

error_value_2 = compute_mean_absolute_error(
    Temperature_43,
    MLR_43,
    Temperature_an,
    analytical_data["MLR 4.3% O2 (g/s)"].values
)

error_value_3 = compute_mean_absolute_error(
    Temperature_82,
    MLR_82,
    Temperature_an,
    analytical_data["MLR 8.2% O2 (g/s)"].values
)

error_value_4 = compute_mean_absolute_error(
    Temperature_205,
    MLR_205,
    Temperature_an,
    analytical_data["MLR 20.5% O2 (g/s)"].values
)

error_value = np.max([error_value_1,
                      error_value_2,
                      error_value_3,
                      error_value_4])

# Define a threshold for validation
threshold = 0.01  # You can adjust this value based on your tolerance

# Print error and exit accordingly
print(f"Validation error (mean of absolute differences ): {error_value:.6f}")
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)




