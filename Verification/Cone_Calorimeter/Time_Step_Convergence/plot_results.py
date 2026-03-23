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
results_file_path = os.path.join(SCRIPT_DIR, 'time_step_convergence.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'time_step_error.png')


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
    Gpyro01mmdt00005=pd.read_csv("0.1mm_dt0.0005s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt0005=pd.read_csv("0.1mm_dt0.005s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt005=pd.read_csv("0.1mm_dt0.05s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt05=pd.read_csv("0.1mm_dt0.5s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt1=pd.read_csv("0.1mm_dt1s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt2=pd.read_csv("0.1mm_dt2s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt5=pd.read_csv("0.1mm_dt5s/reference_cc_summary_01_0001.csv")
    Gpyro01mmdt10=pd.read_csv("0.1mm_dt10s/reference_cc_summary_01_0001.csv")

except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)


style_traces = {
    "dt10":    {'label':"dt=10s"    , "color": "black" , "linestyle": "-" , "linewidth": 3},
    "dt5":     {'label':"dt= 5s"    , "color": "blue"  , "linestyle": "-", "linewidth": 4},
    "dt2":     {'label':"dt= 2s"    , "color": "red"   , "linestyle": "-.", "linewidth": 4},
    "dt1":     {'label':"dt= 1s"    , "color": "green" , "linestyle": ":" , "linewidth": 5},
    "dt05":    {'label':"dt=0.5s"   , "color": "purple", "linestyle": "-" , "linewidth": 3},
    "dt005":   {'label':"dt=0.05s"  , "color": "orange", "linestyle": "--", "linewidth": 3},
    "dt0005":  {'label':"dt=0.005s" , "color": "gray"  , "linestyle": "-.", "linewidth": 3},
    "dt00005": {'label':"dt=0.0005s", "color": "brown" , "linestyle": ":" , "linewidth": 3}
}



# Premier graphique
plt.plot(Gpyro01mmdt10['t'].values, Gpyro01mmdt10['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt10"])
plt.plot(Gpyro01mmdt5['t'].values, Gpyro01mmdt5['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt5"])
plt.plot(Gpyro01mmdt2['t'].values, Gpyro01mmdt2['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values,  **style_traces["dt2"])
plt.plot(Gpyro01mmdt1['t'].values, Gpyro01mmdt1['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt1"])
plt.plot(Gpyro01mmdt05['t'].values, Gpyro01mmdt05['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt05"])
plt.plot(Gpyro01mmdt005['t'].values, Gpyro01mmdt005['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt005"])
plt.plot(Gpyro01mmdt0005['t'].values, Gpyro01mmdt0005['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt0005"])
plt.plot(Gpyro01mmdt00005['t'].values, Gpyro01mmdt00005['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values, **style_traces["dt00005"])



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
    10: Gpyro01mmdt10,
    5: Gpyro01mmdt5,
    2: Gpyro01mmdt2,
    1: Gpyro01mmdt1,
    0.5: Gpyro01mmdt05,
    0.05: Gpyro01mmdt005,
    0.005: Gpyro01mmdt0005,
}

errors = []
dt_values = []

t_ref = Gpyro01mmdt00005['t'].values
mlr_ref = Gpyro01mmdt00005['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values

for dt, df in sorted(cases.items(), reverse=True):
    t = df['t'].values
    mlr = df['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values
    error = compute_mean_absolute_error(t, mlr, t_ref, mlr_ref)
    dt_values.append(dt)
    errors.append(error)

# Conversion en log
log_dt = np.log10(dt_values)
log_err = np.log10(errors)




# Fit line in log-log space
log_dt = np.log10(dt_values)
log_err = np.log10(errors)
slope, intercept, r_value, p_value, std_err = linregress(log_dt, log_err)

# Plot fit line
fit_line = 10**(slope * log_dt + intercept)

plt.figure()
plt.loglog(dt_values, errors, marker='o', linestyle='-', linewidth=4, markersize=15, color='k', label='Simulation')
plt.loglog(dt_values, fit_line, '--', color='red', linewidth=4, label=f'Fit: slope={slope:.2f}, $R^2$={r_value**2:.3f}')
plt.xlabel("Time step (s)")
plt.ylabel("Mean error")
plt.legend(frameon=False, fontsize=30)

#plt.title("Temporal convergence of MLR")
plt.grid(True, which="both", ls="--", lw=0.5)

plt.tight_layout()

# Save figure
plt.savefig(results_file_path2)
#%%
timing_files = {
    10: '0.1mm_dt10s/reference_cc_timing.csv',
    5:  '0.1mm_dt5s/reference_cc_timing.csv',
    2:  '0.1mm_dt2s/reference_cc_timing.csv',
    1:  '0.1mm_dt1s/reference_cc_timing.csv',
    0.5: '0.1mm_dt0.5s/reference_cc_timing.csv',
    0.05: '0.1mm_dt0.05s/reference_cc_timing.csv',
    0.005: '0.1mm_dt0.005s/reference_cc_timing.csv',
    0.0005: '0.1mm_dt0.0005s/reference_cc_timing.csv'
}


def read_cpu_time(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Total CPU Time"):
                return float(line.split(":")[1].strip())
    return None  # or raise an error if needed




cpu_times = []
cpu_dt_values = []

for dt, path in sorted(timing_files.items(), reverse=True):
    full_path = os.path.join(SCRIPT_DIR, path)
    cpu_time = read_cpu_time(full_path)
    if cpu_time is not None:
        cpu_dt_values.append(dt)
        cpu_times.append(cpu_time)
    else:
        print(f"Warning: CPU time not found for dt={dt}")

plt.figure()
plt.loglog(cpu_dt_values, cpu_times, marker='s', linestyle='-', linewidth=4, markersize=15, color='navy')
plt.xlabel("Time step (s)")
plt.ylabel("CPU time (s)")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()


#%%
# -------------------------------------------------------------------
# Export CSV summary for LaTeX
# -------------------------------------------------------------------

csv_output_path = os.path.join(SCRIPT_DIR, 'convergence_summary.csv')

with open(csv_output_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Time step (s)', 'Mean Error', 'CPU Time (s)'])
    for dt in sorted(dt_values + [0.0005]):  # assure l’ordre croissant
        error = next((e for d, e in zip(dt_values, errors) if d == dt), None)
        cpu = next((c for d, c in zip(cpu_dt_values, cpu_times) if d == dt), None)
        writer.writerow([dt, f"{error:.2e}" if error else "--", f"{cpu:.2f}" if cpu else ""])





#%%


order_ok = 0.95 <= abs(slope) <= 1.05
fit_quality_ok = r_value**2 > 0.95

print("Order=",slope, "R²=",r_value**2)
if order_ok and fit_quality_ok:
    print("Temporal convergence confirmed: order ≈ 1 and R² > 0.95.")
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Temporal convergence not confirmed.")
    print("Order=",slope)
    print("Validation FAILED.")
    sys.exit(1)

#%%
                   

