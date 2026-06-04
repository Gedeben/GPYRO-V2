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


results_file_path1 = os.path.join(SCRIPT_DIR, 'Total_mass_flux_instantaneous.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'profile_mass_flux_instantaneous.png')


# File paths
file_path = os.path.join(SCRIPT_DIR, 'gas_transport_summary_01_0001.csv')

# --- NEW STYLING PARAMETERS ---
# Define colors for different quantities and styles for the data source
plot_colors = ['blue', 'red', 'black', 'green']
gpyro_style = {'linewidth': 6, 'linestyle': '-'}
analytical_style = {'linewidth': 6, 'linestyle': '--'}

# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

#%%

MFLUX_IN = 100 #g.m2.s-1
rho= 1000 #kg.m3
m_ini=1 #kg.m3
Z = 0.001 # s-1
theta=0.40
reaction_rate= Z*rho
L=0.1

Tend= rho/reaction_rate

omega_fg= theta*reaction_rate #Kg.m3.s-1

# Load data
try:
    data = pd.read_csv(file_path)
except FileNotFoundError as e:
    print(f"Error: The data file was not found at '{file_path}'. {e}")
    sys.exit(1)

# Time vector
time = data['t'].values

time_an= np.linspace(0, time[-1],50)

MF_an=[]

for t in time_an:
    if t<Tend:
        MF_an.append(-MFLUX_IN-1000*omega_fg*L)
    else:
        MF_an.append(-MFLUX_IN)
    


plt.figure()
#plt.plot(time, data["004_MLR( 0.0000_ 0.0000_ 0.0000)"].values, color='red', label="Gpyro MLR", lw=3)
plt.plot(time, data["005_MASS_FLUX_TOTAL_Z( 0.0000_ 0.0000_ 0.0000)"].values, lw=3, color='blue', label="Gpyro ")
plt.plot(time_an, MF_an, color='blue',markersize=15, marker = '+',markeredgewidth=3,linestyle=' ', label="Analytical")

plt.xlabel('Time [s]')
plt.ylabel(r'Mass Flux [g·m$^{-2}$·s$^{-1}$]')

plt.legend(frameon=False, fontsize=30)
ymin, ymax = plt.ylim()
#plt.yticks(np.linspace(-140, 40, 5))

plt.tight_layout()
plt.savefig(results_file_path1)


#%%

dflux=pd.read_csv("gas_transport_profile_01_01_MASS_FLUX_TOTAL_Z.csv", header= None)

Ntime=100

depths = dflux.iloc[0, 1:].astype(float).values

Time= dflux.iloc[Ntime, 0]

MF=dflux.iloc[Ntime, 1:].values




plt.figure()


MF_th= -(1000*(depths[-1]-depths[:])*omega_fg+MFLUX_IN)
        

plt.plot(100*depths[:],MF_th, color='blue', markersize=15, marker = '+',
         markeredgewidth=3, linestyle=' ',markevery=3, label="Analytical")

plt.plot(100*depths,MF, c='b',label='Gpyro', lw=3)

plt.xlabel('Depth [mm]')
plt.ylabel(r'Mass Flux [g·m$^{-2}$·s$^{-1}$]')

plt.legend(frameon=False, fontsize=30)

#plt.title("Mass flux through the sample")

# --- Ajout des flèches rouges ---

# Point top surface (depth = 0)
plt.annotate(
    "top surface",
    xy=(0, MF[0]),              # point à annoter
    xytext=(3, MF[0] - 0.02*max(MF)),  # position du texte
    arrowprops=dict(arrowstyle="->", color='red'),
    color='red',
    fontsize=24
)

# Point bottom surface (depth = max(depths))
plt.annotate(
    "bottom surface",
    xy=(100*depths[-1], MF[-1]),
    xytext=(100*depths[-1] - 2.6, MF[-1] + 0.15*max(MF)),
    arrowprops=dict(arrowstyle="->", color='red'),
    color='red',
    fontsize=24
)
plt.tight_layout()

plt.savefig(results_file_path2)



#%%



error_value=np.mean(abs(MF-MF_th))/MFLUX_IN

print(f"Mean relative error {np.round(error_value,3)} %")
# Validation check
threshold = 0.1
if error_value <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)