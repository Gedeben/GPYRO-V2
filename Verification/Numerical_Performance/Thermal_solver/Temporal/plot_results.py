import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import linregress
import csv


import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))


results_file_path1 = os.path.join(SCRIPT_DIR, 'Temperature_Temporal_convergence.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'L2_error_Temporal_thermal_solver.png')
results_file_path3 = os.path.join(SCRIPT_DIR, 'CPU_Time_Temporal_thermal_solver.png')


# Load data
try:
    Gpyrodt100=pd.read_csv("dt_100s/radiation_1d_summary_01_0001.csv")
    #Gpyrodt50=pd.read_csv("dt_50s/radiation_1d_summary_01_0001.csv")
    Gpyrodt10=pd.read_csv("dt_10s/radiation_1d_summary_01_0001.csv")
    Gpyrodt1=pd.read_csv("dt_1s/radiation_1d_summary_01_0001.csv")
    Gpyrodt01=pd.read_csv("dt_01s/radiation_1d_summary_01_0001.csv")
    Gpyrodt001=pd.read_csv("dt_001s/radiation_1d_summary_01_0001.csv")
    Gpyrodt0001=pd.read_csv("dt_0001s/radiation_1d_summary_01_0001.csv")
    #Gpyrodt00001=pd.read_csv("dt_00001s/radiation_1d_summary_01_0001.csv")


except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)



#%% color and style parameters for the plots in black and white
style_traces = {
    "dt100":    {'label':"dt=100s"    , "color": "darkblue" , "linestyle": "-" , "linewidth": 3},
    "dt50":     {'label':"dt= 50s"    , "color": "blue"  , "linestyle": "-", "linewidth": 4},
    "dt10":     {'label':"dt= 10s"    , "color": "red"   , "linestyle": "-.", "linewidth": 4},
    "dt1":     {'label':"dt= 1s"    , "color": "green" , "linestyle": ":" , "linewidth": 5},
    "dt01":    {'label':"dt=0.1s"   , "color": "purple", "linestyle": "-" , "linewidth": 3},
    "dt001":   {'label':"dt=0.01s"  , "color": "orange", "linestyle": "--", "linewidth": 3},
    "dt0001":  {'label':"dt=0.001s" , "color": "gray"  , "linestyle": "-.", "linewidth": 3},
    "dt00005": {'label':"dt=0.0005s", "color": "brown" , "linestyle": ":" , "linewidth": 3},
    "an":      {'label':"analytique", "color": "black" , "linestyle": "-", 'marker':'+' , "linewidth": 3}
}



# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)


#%%

# Premier graphique
plt.plot(Gpyrodt100['t'].values, Gpyrodt100['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt100"])
#plt.plot(Gpyrodt50['t'].values, Gpyrodt50['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt50"])
plt.plot(Gpyrodt10['t'].values, Gpyrodt10['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt10"])
plt.plot(Gpyrodt1['t'].values, Gpyrodt1['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values,  **style_traces["dt1"])
plt.plot(Gpyrodt01['t'].values, Gpyrodt01['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt01"])
plt.plot(Gpyrodt001['t'].values, Gpyrodt001['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt001"])
#plt.plot(Gpyrodt0001['t'].values, Gpyrodt0001['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values, **style_traces["dt0001"])

# Mise en forme
plt.xlabel("Time (s)")
plt.ylabel(r"Temperature (°C)")
#plt.xlim([0, 5000])
#plt.ylim([0, 25])
plt.grid(False) #, linestyle='--', linewidth=0.5, alpha=0.7)  # Grille en pointillé
plt.legend(frameon=False)  # Légende sans cadre
plt.tight_layout()  # Ajuste automatiquement la disposition


# Save figure
plt.savefig(results_file_path1)




#%%


# Norme L1 : erreur moyenne absolue
def L1(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    # Grille commune sur l'intervalle d'intersection
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    # Interpolation des deux courbes
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    # Calcul de la norme L1
    return np.mean(np.abs(sim_interp - ref_interp))


# Norme L2 : erreur quadratique moyenne
def L2(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    return np.sqrt(np.mean((sim_interp - ref_interp)**2))


# Norme Linfini : erreur maximale
def Linf(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    return np.max(np.abs(sim_interp - ref_interp))


cases = {
    100: Gpyrodt100,
    #50: Gpyrodt50,
    10: Gpyrodt10,
    1: Gpyrodt1,
    0.1: Gpyrodt01,
    0.01: Gpyrodt001,
    #0.001: Gpyrodt0001

    }

errors = []
dt_values = []

an_t=Gpyrodt0001['t'].values
an_MLR= Gpyrodt0001['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values


for dt, df in sorted(cases.items(), reverse=True):
    t = df['t'].values
    mlr = df['003_TEMPERATURE( 0.0080_ 0.0000_ 0.0000)'].values
    error = L2(t, mlr, an_t, an_MLR)
    dt_values.append(dt)
    errors.append(error)


# Conversion en log
log_dt = np.log10(dt_values)
log_err = np.log10(errors)




# Fit line in log-log space
log_dt = np.log10(dt_values)
log_err = np.log10(errors)
slope, intercept, r_value, p_value, std_err = linregress(log_dt[1:], log_err[1:])

# Plot fit line
fit_line = 10**(slope * log_dt + intercept)

plt.figure()
plt.loglog(dt_values, errors, marker='o', linestyle='-', linewidth=4, markersize=15, color='k', label=r'Mesured $L^2 error$')
plt.loglog(dt_values, fit_line, '--', color='red', linewidth=4, label=f'Fit: slope={slope:.2f}, $R^2$={r_value**2:.3f}')
plt.xlabel("Time step (s)")
plt.ylabel("Error")
plt.legend(frameon=False, fontsize=30)

#plt.title("Temporal convergence of MLR")
plt.grid(True, which="both", ls="--", lw=0.5)

plt.tight_layout()
plt.savefig(results_file_path2)





#%%





timing_files = {
    100: 'dt_100s/radiation_1d_timing_01.csv',
    10:  'dt_10s/radiation_1d_timing_01.csv',
    1:  'dt_1s/radiation_1d_timing_01.csv',
    0.1: 'dt_01s/radiation_1d_timing_01.csv',
    0.01: 'dt_001s/radiation_1d_timing_01.csv',
    0.001: 'dt_0001s/radiation_1d_timing_01.csv',
    }

timing_files = {
    100: 'dt_100s/radiation_1d_timing.csv',
    10:  'dt_10s/radiation_1d_timing.csv',
    1:  'dt_1s/radiation_1d_timing.csv',
    0.1: 'dt_01s/radiation_1d_timing.csv',
    0.01: 'dt_001s/radiation_1d_timing.csv',
    0.001: 'dt_0001s/radiation_1d_timing.csv',
    }



def read_cpu_time(file_path):
    # --- Première tentative : lecture classique ---
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Total CPU Time"):
                try:
                    return float(line.split(":")[1].strip())
                except (IndexError, ValueError):
                    pass

    # --- Deuxième tentative : lecture CSV ---
    total = 0.0
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if len(row) >= 2:
                    try:
                        total += float(row[1])
                    except ValueError:
                        pass
        return total
    except Exception:
        return None



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


#data_old_Gpyro=pd.read_csv("convergence_summary_old.csv")
#cpu_dt_values_old=data_old_Gpyro["Time step (s)"].values
#cpu_time_old=data_old_Gpyro["CPU Time (s)"].values


# Fit line in log-log space
log_dt2 = np.log10(cpu_dt_values)
log_cpu = np.log10(cpu_times)
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(log_dt2[3:], log_cpu[3:])

# Plot fit line
fit_line2 = 10**(slope2 * log_dt2 + intercept2)

plt.figure()
#plt.loglog(cpu_dt_values_old, cpu_time_old, marker='s', linestyle='-', linewidth=4, markersize=15, color='green', label="Gpyro V0.8")

plt.loglog(cpu_dt_values, cpu_times, marker='s', linestyle='-', linewidth=4, markersize=15, color='navy', label="simulation")
plt.loglog(cpu_dt_values[1:], fit_line2[1:], '--', color='red', linewidth=4, label=f'Fit: slope={slope2:.2f}, $R^2$={r_value2**2:.3f}')

plt.xlabel("Time step (s)")
plt.ylabel("CPU time (s)")
plt.legend(frameon=False, fontsize=30)
plt.minorticks_on()

plt.grid(True, which="both", ls="--", lw=0.5)

plt.savefig(results_file_path3)


#%%
"""
# -------------------------------------------------------------------
# Export CSV summary 
# -------------------------------------------------------------------

csv_output_path = os.path.join(SCRIPT_DIR, 'convergence_summary.csv')

with open(csv_output_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Time step (s)', 'Mean Error', 'CPU Time (s)'])
    for dt in sorted(cpu_dt_values):  # assure l’ordre croissant
        error = next((e for d, e in zip(dt_values, errors) if d == dt), None)
        cpu = next((c for d, c in zip(cpu_dt_values, cpu_times) if d == dt), None)
        writer.writerow([dt, f"{error:.2e}" if error else "--", f"{cpu:.5f}" if cpu else ""])

"""



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