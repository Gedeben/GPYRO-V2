
import os
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
import sys



SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = SCRIPT_DIR  
# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

def read_cpu_time(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Total CPU Time"):
                return float(line.split(":")[1].strip())
    return None



style_traces = {
    1:  { "color": "black" , "linestyle": "-" ,  "linewidth": 3 },
    2:  { "color": "blue"  , "linestyle": "-" ,  "linewidth": 4 },
    3:  { "color": "red"   , "linestyle": "-.", "linewidth": 4 },
    4:  { "color": "green" , "linestyle": ":" ,  "linewidth": 5 },
    5:  { "color": "purple", "linestyle": "-" ,  "linewidth": 3 },
    6:  { "color": "orange", "linestyle": "--",  "linewidth": 3 },
    7:  { "color": "gray"  , "linestyle": "-.", "linewidth": 3 },
    8:  { "color": "brown" , "linestyle": ":" ,  "linewidth": 3 },
    9:  { "color": "cyan"  , "linestyle": "--",  "linewidth": 3 },
    10: { "color": "magenta", "linestyle": "-.", "linewidth": 3 },
}

omp_thread_counts = []
cpu_times = []

plt.figure(figsize=(12,8))  # Figure pour les MLR

# Boucle sur les sous-dossiers OMP_*
for folder in os.listdir(BASE_DIR):
    folder_path = os.path.join(BASE_DIR, folder)
    if os.path.isdir(folder_path) and folder.startswith("OMP_"):
        try:
            thread_count = int(folder.split("_")[1])
        except (IndexError, ValueError):
            print(f"Skipping folder '{folder}' (invalid name format)")
            continue

        # Lire le temps CPU
        for filename in os.listdir(folder_path):
            if filename.endswith("_timing.csv"):
                full_path = os.path.join(folder_path, filename)
                cpu_time = read_cpu_time(full_path)
                if cpu_time is not None:
                    omp_thread_counts.append(thread_count)
                    cpu_times.append(cpu_time)
                break

        # Lire et tracer le fichier summary
        summary_file = os.path.join(folder_path, "reference_cc_summary_01_0001.csv")
        if os.path.isfile(summary_file):
            try:
                Gpyro = pd.read_csv(summary_file)
                style = style_traces.get(thread_count, { "color": "black", "linestyle": "--", "linewidth": 2 })
                label = f"{thread_count} thread{'s' if thread_count > 1 else ''}"
                plt.plot(
                    Gpyro['t'].values,
                    Gpyro['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values,
                    label=label,
                    **style
                )
            except Exception as e:
                print(f"Error reading {summary_file}: {e}")
        else:
            print(f"Summary file not found for {folder}")

# Affichage du MLR
plt.xlabel("Time (s)")
plt.ylabel(r"MLR (g.m$^{-2}$.s$^{-1}$)")
plt.grid(False)
plt.legend(frameon=False)  # Légende sans cadre
plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_DIR, "omp_mlr_comparison.png"))
#%%


# Tris des données
omp_thread_counts, cpu_times = zip(*sorted(zip(omp_thread_counts, cpu_times)))
omp_thread_counts = np.array(omp_thread_counts)
cpu_times = np.array(cpu_times)

# Calculs
normalized_cpu_times = cpu_times / cpu_times[0]
speedups = cpu_times[0] / cpu_times

# Modèle d'Amdahl
def amdahl(N, f):
    return 1 / ((1 - f) + f / N)

# Ajustement pour estimer f
popt_amdahl, _ = curve_fit(amdahl, omp_thread_counts, speedups, bounds=(0, 1))
f_amdahl = popt_amdahl[0]
N_fine = np.linspace(min(omp_thread_counts), max(omp_thread_counts), 100)

# 📊 Création du graphique avec deux axes y
fig, ax1 = plt.subplots(figsize=(12,8))

# ➤ Axe Y1 : CPU Time normalisé
color1 = 'tab:red'
ax1.set_xlabel("Number of OpenMP Threads")
ax1.set_ylabel("Normalized CPU Time", color=color1)
ax1.plot(omp_thread_counts, normalized_cpu_times, marker='o', linestyle='-', linewidth=3, markersize=10, color=color1, label="Normalized CPU Time")
ax1.tick_params(axis='y', labelcolor=color1)

# ➤ Axe Y2 : Speedup + Amdahl
ax2 = ax1.twinx()
ax2.set_ylabel("Speedup", color='navy')
ax2.plot(omp_thread_counts, speedups, 'o', markersize=10,color='navy', label="Measured Speedup")
ax2.plot(N_fine, amdahl(N_fine, f_amdahl), '-', color='navy', label=f"Amdahl (f={f_amdahl:.2f})")
ax2.tick_params(axis='y', labelcolor='navy')

# ➤ Grille + légende combinée
fig.tight_layout()
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(
    lines_1 + lines_2, labels_1 + labels_2,
    borderaxespad=0,
    bbox_to_anchor=(0.14, 0.75),  # décale vers la droite (1.02 > 1)
    frameon=False,
    fontsize=24
)

#ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left',frameon=False, fontsize=24)  # Légende sans cadre

plt.grid(False)

# 💾 Sauvegarde
plt.savefig(os.path.join(SCRIPT_DIR, "omp_scaling.png"))


validation_results1 = f_amdahl > 0.65


if validation_results1:
    print(f"✅ Parallelization successful: speedup from 1 to 2 CPUs = {speedups[1]:.2f}")
else:
    print(f"❌ Parallelization not efficient enough: speedup from 1 to 2 CPUs = {speedups[1]:.2f}")


#%%
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



# ---------- VALIDATION DES COURBES MLR ----------

reference_thread = 1
validation_results2 = True
threshold = 1e-5
reference_thread = 1

# Collecte toutes les courbes pour comparaison
mlr_curves = {}

for folder in os.listdir(BASE_DIR):
    folder_path = os.path.join(BASE_DIR, folder)
    if os.path.isdir(folder_path) and folder.startswith("OMP_"):
        try:
            thread_count = int(folder.split("_")[1])
        except:
            continue

        summary_file = os.path.join(folder_path, "reference_cc_summary_01_0001.csv")
        if os.path.isfile(summary_file):
            df = pd.read_csv(summary_file)
            times = df["t"].values
            mlr = df["001_MLR( 0.0000_ 0.0000_ 0.0000)"].values
            mlr_curves[thread_count] = (times, mlr)

# Définir la courbe de référence (ex: 1 thread)
if reference_thread in mlr_curves:
    ref_times, ref_mlr = mlr_curves[reference_thread]

    for thread_count, (times, mlr) in mlr_curves.items():
        if thread_count == reference_thread:
            continue

        mae = compute_mean_absolute_error(times, mlr, ref_times, ref_mlr)
        if mae >=  threshold:
            validation_results2 = False
            break   


else:
    print("⚠️  Reference thread data not found for validation.")
    validation_results2 = False



#%%



if validation_results2 and validation_results1:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)

#%%
                   

