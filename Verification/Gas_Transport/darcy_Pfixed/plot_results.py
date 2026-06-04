import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import numpy as np
import os
import sys

# If running as a script, SCRIPT_DIR will be the script's directory.
# For interactive environments, we'll assume the data file is in the current working directory.
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
except NameError:
    SCRIPT_DIR = os.getcwd()


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))


results_file_path1 = os.path.join(SCRIPT_DIR, 'Viscosity_darcy_fixed.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'Pressure_darcy_fixed.png')
results_file_path3 = os.path.join(SCRIPT_DIR, 'flux_darcy_fixed.png')




# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

param_sets = [
    {'color': 'blue' , 'linewidth': 4, 'linestyle': '-'},
    {'color': 'blue' , 'marker': '+', 'linestyle': 'none', 'markersize': 15, 'markevery': 40,  'markeredgewidth': 3},
    {'color': 'green' , 'linewidth': 2, 'linestyle': '-.'},
    {'color': 'orange' , 'linewidth': 2, 'linestyle': ':'},
    {'color': 'blue' , 'linewidth': 2, 'linestyle': '', 'marker': 'o', 'markersize': 6, 'markevery': 10},
    {'color': 'red', 'linewidth': 2, 'linestyle': '', 'marker': 's', 'markersize': 6, 'markevery': 10},
    {'color': 'green' , 'linewidth': 2, 'linestyle': '', 'marker': '^', 'markersize': 6, 'markevery': 10},
    {'color': 'orange' , 'linewidth': 2, 'linestyle': '', 'marker': 'd', 'markersize': 6, 'markevery': 10}
 ]

#%%
# Load data
try:
    Pl=pd.read_csv("darcy_profile_01_01_PRESSURE.csv", header= None)
    Dl=pd.read_csv("darcy_profile_01_04_D12.csv", header= None)
    MFl=pd.read_csv("darcy_profile_01_05_MASS_FLUX_TOTAL_Z.csv", header= None)
    Tl=pd.read_csv("darcy_profile_01_03_TEMPERATURE.csv", header= None)
except FileNotFoundError:
    print("Error: The data file was not found ")
    sys.exit(1)
Pref=101300


Ntime=50
depths = Pl.iloc[0, 1:].astype(float).values
Time= Pl.iloc[Ntime, 0]
P=Pl.iloc[Ntime, 1:].values/Pref+1 # in Pa and absolute 
D=Dl.iloc[Ntime, 1:].values
MF=MFl.iloc[Ntime, 1:].values
T=Tl.iloc[Ntime, 1:].values

T0=300
M=29
sigma=3.617
EPSOK=97
Te=T0/EPSOK
a1=1.080794 
a2=-0.16033 
a3=0.605009
a4=-0.88524 
a5= 2.115672
a6=-2.98308
K=5E-10
Gamma=a1*Te**a2+a3*np.exp(a4*Te)+a5*np.exp(a6*Te)


alpha= 0.018829*(T0**3*2/M)**(1/2)/(sigma**2*Gamma)


P0=1
P1=2
L=0.10

n=2

P_an=(P0**n+(P1**n-P0**n)*depths/L)**(1/n)


# Calcul de la valeur constante
mf_value = K / (2 * alpha * L) * ((P0*Pref)**2 - (P1*Pref)**2)

# Liste/array de la même longueur que depths
MF_an = np.full_like(depths, mf_value, dtype=float)

#%%

D_an=alpha/(P_an*Pref)
plt.figure()
plt.plot(100*depths,D, label="Gpyro", **param_sets[0])
plt.plot(100*depths,D_an, label="Analytical", **param_sets[1])

plt.xlabel('Depth [cm]')
plt.ylabel('kinematic viscosity [m².s⁻¹]')

plt.legend(frameon=False, fontsize=30)
plt.tight_layout()

plt.savefig(results_file_path1)


#%%


plt.figure()

plt.plot(100*depths,P, label="Gpyro", **param_sets[0])
plt.plot(100*depths,P_an, label="Analytical", **param_sets[1])


plt.xlabel('Depth [cm]')
plt.ylabel(r'Pressure [bar]')

plt.legend(frameon=False, fontsize=30)


# --- Ajout des flèches rouges ---

# Point top surface (depth = 0)
plt.annotate(
    "top surface",
    xy=(0, P[0]),              # point à annoter
    xytext=(3, P[0] + 0.02*max(P)),  # position du texte
    arrowprops=dict(arrowstyle="->", color='red'),
    color='red',
    fontsize=24
)

# Point bottom surface (depth = max(depths))
plt.annotate(
    "bottom surface",
    xy=(100*depths[-1], P[-1]),
    xytext=(100*depths[-1] - 2.6, P[-1] - 0.15*max(P)),
    arrowprops=dict(arrowstyle="->", color='red'),
    color='red',
    fontsize=24
)
plt.tight_layout()
plt.savefig(results_file_path2)


#%%


plt.figure()


plt.plot(100*depths,MF/1000, label="Gpyro", **param_sets[0])
plt.plot(100*depths,MF_an, label="Analytical", **param_sets[1])


plt.ylim([0.90*mf_value, 1.10*mf_value])
plt.xlabel('Depth [cm]')
plt.ylabel(r'mass Flux [kg.m².s⁻¹]')

plt.legend(frameon=False, fontsize=30)
plt.tight_layout()

plt.savefig(results_file_path3)



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




P_error = L2(depths, P, depths, P_an)/np.mean(P_an)
D_error = L2(depths, D, depths, D_an)/ np.mean(D_an)
MF_error = L2(depths, MF/1000, depths, MF_an)/mf_value



error_value=np.max([P_error, D_error, MF_error])

print(f"Mass Flux relative error {np.round(MF_error,4)} %")
print(f"Viscosity relative error {np.round(D_error,4)} %")
print(f"Pressure relative error {np.round(P_error,4)} %")
# Validation check
threshold = 0.005
if (error_value <= threshold)  :
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)