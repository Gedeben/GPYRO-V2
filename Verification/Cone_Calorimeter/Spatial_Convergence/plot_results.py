#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
import matplotlib as mpl
import os
import sys
from scipy.stats import linregress
import csv


# Determine the current script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path = os.path.join(SCRIPT_DIR, 'grid_convergence.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'grid_error.png')


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)






# Load data
try:
    Gpyro1mm=pd.read_csv("1mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro05mm=pd.read_csv("0.5mm_dt0.05s/reference_cc_summary_01_0001.csv") 
    Gpyro02mm=pd.read_csv("0.2mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro01mm=pd.read_csv("0.1mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro0025mm=pd.read_csv("0.025mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro001mm=pd.read_csv("0.01mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro0001mm=pd.read_csv("0.001mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro00001mm=pd.read_csv("0.0001mm_dt0.05s/reference_cc_summary_01_0001.csv")

except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)


style_traces = {
    "dz1":     {'label':"dz= 1 mm"    , "color": "black" , "linestyle": "-" , "linewidth": 3},
    "dz05":    {'label':"dz= 0.5 mm"  , "color": "blue"  , "linestyle": "-", "linewidth": 4},
    "dz02":    {'label':"dz=0.2 mm"   , "color": "red"   , "linestyle": "-.", "linewidth": 4},
    "dz01":    {'label':"dz=0.1mm"    , "color": "green" , "linestyle": ":" , "linewidth": 5},
    "dz0025":  {'label':"dz=0.025 mm" , "color": "purple", "linestyle": "-" , "linewidth": 3},
    "dz001":   {'label':"dz=0.01 mm"  , "color": "orange", "linestyle": "--", "linewidth": 3},
    "dz0001":  {'label':"dz=0.001 mm" , "color": "gray"  , "linestyle": "-.", "linewidth": 3},
    "dz00001": {'label':"dz=0.0001 mm", "color": "brown" , "linestyle": ":" , "linewidth": 3}
}



# Premier graphique
plt.plot(Gpyro1mm['t'].values, Gpyro1mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz1"])
plt.plot(Gpyro05mm['t'].values, Gpyro05mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz05"])
plt.plot(Gpyro02mm['t'].values, Gpyro02mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values,  **style_traces["dz02"])
plt.plot(Gpyro01mm['t'].values, Gpyro01mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz01"])
plt.plot(Gpyro0025mm['t'].values, Gpyro0025mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz0025"])
plt.plot(Gpyro001mm['t'].values, Gpyro001mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz001"])
plt.plot(Gpyro0001mm['t'].values, Gpyro0001mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz0001"])
plt.plot(Gpyro00001mm['t'].values, Gpyro00001mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dz00001"])


# Mise en forme
plt.xlabel("Time (s)")
plt.ylabel(r"MLR (g.m$^{-2}$.s$^{-1}$)")
plt.xlim([0, 600])
#plt.ylim([0, 25])
plt.grid(False) #, linestyle='--', linewidth=0.5, alpha=0.7)  # Grille en pointillé
plt.legend(frameon=False)  # Légende sans cadre
plt.tight_layout()  # Ajuste automatiquement la disposition


# Save figure
plt.savefig(results_file_path)



#%%

# -------------------------------------------------------------------
# Compute error between analytical and simulation for each curve
# -------------------------------------------------------------------

def compute_mean_absolute_error(sim_times, sim_values, ref_times, ref_values, num_points=500):
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


cases = {
    1: Gpyro1mm,
    0.5: Gpyro05mm,
    0.2: Gpyro02mm,
    0.1: Gpyro01mm,
    0.025: Gpyro0025mm,
    0.01: Gpyro001mm,
    0.001: Gpyro0001mm,
}

errors = []
dz_values = []

t_ref = Gpyro00001mm['t'].values
mlr_ref = Gpyro00001mm['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values
for dz, df in sorted(cases.items(), reverse=True):
    t = df['t'].values
    mlr = df['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values
    error = compute_mean_absolute_error(t, mlr, t_ref, mlr_ref)
    dz_values.append(dz)
    errors.append(error)
# Conversion en log
log_dz = np.log10(dz_values)
log_err = np.log10(errors)


slope, intercept, r_value, p_value, std_err = linregress(log_dz, log_err)

# Plot fit line
fit_line = 10**(slope * log_dz + intercept)

plt.figure()
plt.loglog(dz_values, errors, marker='o', linestyle='-', linewidth=4, markersize=15, color='k', label='Simulation')
plt.loglog(dz_values, fit_line, '--', color='red', linewidth=4, label=f'Fit: slope={slope:.2f}, $R^2$={r_value**2:.3f}')
plt.xlabel(r"$\Delta Z$ (mm)")
plt.ylabel("Mean error")
plt.legend(frameon=False, fontsize=30)

#plt.title("Temporal convergence of MLR")
plt.grid(True, which="both", ls="--", lw=0.5)

plt.tight_layout()

# Save figure
plt.savefig(results_file_path2)
#%%

timing_files = {
    1:     "1mm_dt0.05s/reference_cc_timing.csv",
    0.5:   "0.5mm_dt0.05s/reference_cc_timing.csv",
    0.2:   "0.2mm_dt0.05s/reference_cc_timing.csv",
    0.1:   "0.1mm_dt0.05s/reference_cc_timing.csv",
    0.025: "0.025mm_dt0.05s/reference_cc_timing.csv",
    0.01:  "0.01mm_dt0.05s/reference_cc_timing.csv",
    0.001: "0.001mm_dt0.05s/reference_cc_timing.csv",
    0.0001:"0.0001mm_dt0.05s/reference_cc_timing.csv"
}


def read_cpu_time(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Total CPU Time"):
                return float(line.split(":")[1].strip())
    return None  # or raise an error if needed




cpu_times = []
cpu_dz_values = []

for dz, path in sorted(timing_files.items(), reverse=True):
    full_path = os.path.join(SCRIPT_DIR, path)
    cpu_time = read_cpu_time(full_path)
    if cpu_time is not None:
        cpu_dz_values.append(dz)
        cpu_times.append(cpu_time)
    else:
        print(f"Warning: CPU time not found for dz={dz}")

plt.figure()
plt.loglog(cpu_dz_values, cpu_times, marker='s', linestyle='-', linewidth=4, markersize=15, color='navy')
plt.xlabel("\Delta z (mm)")
plt.ylabel("CPU time (s)")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()



# -------------------------------------------------------------------
# Export CSV summary for LaTeX
# -------------------------------------------------------------------

csv_output_path = os.path.join(SCRIPT_DIR, 'convergence_summary.csv')

with open(csv_output_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['dz (mm)', 'Mean Error', 'CPU Time (s)'])
    for dz in sorted(dz_values+[0.0001] ):  # assure l’ordre croissant

        error = next((e for d, e in zip(dz_values, errors) if d == dz), None)
        cpu = next((c for d, c in zip(cpu_dz_values, cpu_times) if d == dz), None)
        writer.writerow([dz, f"{error:.2e}" if error else "--", f"{cpu:.2f}" if cpu else ""])





#%%


order_ok = 0.95 <= abs(slope) <= 1.05
fit_quality_ok = r_value**2 > 0.95

print("Order=",slope, "R²=",r_value**2)
if order_ok and fit_quality_ok:
    print("spatial convergence confirmed: order ≈ 1 and R² > 0.95.")
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("spatial convergence not confirmed.")
    print("Order=",slope)
    print("Validation FAILED.")
    sys.exit(1)

#%%
                   

