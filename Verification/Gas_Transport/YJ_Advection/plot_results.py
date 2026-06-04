import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys
import re
from matplotlib import cm



# If running as a script, SCRIPT_DIR will be the script's directory.
# For interactive environments, we'll assume the data file is in the current working directory.
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
except NameError:
    SCRIPT_DIR = os.getcwd()


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))


results_file_path1 = os.path.join(SCRIPT_DIR, 'O2_advection.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'N2_advection.png')


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 26  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (12, 8)

param_sets = [
    {'color': 'blue' , 'linewidth': 3, 'linestyle': '-'},
    {'color': 'red', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'blue' , 'linewidth': 2, 'linestyle': '', 'marker': '+','markersize': 15, 'markevery': 5,'markeredgewidth':3},
    {'color': 'red' , 'linewidth': 2, 'linestyle': '','marker': '+','markersize': 15,  'markevery': 5,'markeredgewidth':3},
    {'color': 'blue' , 'linewidth': 2, 'linestyle': '', 'marker': 'o', 'markersize': 6, 'markevery': 1},
    {'color': 'red', 'linewidth': 2, 'linestyle': '', 'marker': 's', 'markersize': 6, 'markevery': 10},
    {'color': 'green' , 'linewidth': 2, 'linestyle': '', 'marker': '^', 'markersize': 6, 'markevery': 10},
    {'color': 'orange' , 'linewidth': 2, 'linestyle': '', 'marker': 'd', 'markersize': 6, 'markevery': 10}
 ]




#%%


# Norme L2 : erreur quadratique moyenne
def L2(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    return np.sqrt(np.mean((sim_interp - ref_interp)**2))

def fix_fortran_number(s):
    if isinstance(s, float) or isinstance(s, int):
        return s
    s = s.strip()
    # si la chaîne est déjà correcte -> OK
    if "E" in s or "e" in s:
        return float(s)
    # sinon -> injecter un "E" avant le dernier signe +/-
    s = re.sub(r'([0-9])([+-][0-9]+)$', r'\1E\2', s)
    return float(s)


def Yj_advection_analytical(z, t, v, Y0, Yinf):
    """
    Analytical solution of a 1D pure advection problem
    with a moving step profile.

    Parameters
    ----------
    z : array_like
        Spatial coordinate
    t : float
        Time
    v : float
        Advection velocity
    Y0 : float
        Initial concentration
    Yinf : float
        Advected concentration

    Returns
    -------
    Y : ndarray
        Analytical solution
    """
    z_front = v * t
    Y = np.where(z[-1]-z < z_front, Yinf, Y0)
    return Y


def front_position(z, Y, Y0, Yinf):
    """
    Compute the front position defined by Y = (Y0+Yinf)/2
    using linear interpolation.
    """
    Ymid = 0.5 * (Y0 + Yinf)

    # différence par rapport au seuil
    dY = Y - Ymid

    # indices où le signe change
    idx = np.where(np.sign(dY[:-1]) != np.sign(dY[1:]))[0]

    if len(idx) == 0:
        return np.nan  # front hors domaine

    i = idx[0]

    # interpolation linéaire
    zf = z[i] + (z[i+1] - z[i]) * (Ymid - Y[i]) / (Y[i+1] - Y[i])
    return zf


def front_velocity(t, zf):
    """
    Compute front velocity from front position vs time.
    """
    t = np.array(t)
    zf = np.array(zf)

    v = np.diff(zf) / np.diff(t)
    t_mid = 0.5 * (t[1:] + t[:-1])

    return t_mid, v

#%%
# Load data
try:
    YJ1l=pd.read_csv("gas_advection_profile_01_01_YJ(01).csv", header= None)
    YJ2l=pd.read_csv("gas_advection_profile_01_02_YJ(02).csv", header= None)
    Gpyro=pd.read_csv("gas_advection_summary_01_0001.csv")
except FileNotFoundError:
    print("Error: The data file was not found ")
    sys.exit(1)
Pref=101300




T0=300

a1=1.080794 
a2=-0.16033 
a3=0.605009
a4=-0.88524 
a5= 2.115672
a6=-2.98308
Tref=300
R=8.314



M1=32
M2=32
sigma1=3433
sigma2=3681
EPSOK1=113
EPSOK2=91.5

sigma=0.5*(sigma1+sigma2)
EPSOK=(EPSOK1*EPSOK2)**0.5
Te=T0/EPSOK

C1=1/M1 +1/M2

Gamma=a1*Te**a2+a3*np.exp(a4*Te)+a5*np.exp(a6*Te)
alpha1= 0.018829*(T0**3*C1)**(1/2)/(sigma**2*Gamma)
D=alpha1/Pref

Ntime=50

rhogl=pd.read_csv("gas_advection_profile_01_06_GAS_DENSITY.csv", header= None)
rhog = rhogl.iloc[Ntime, 1:].values


#D=2.07016517e-05

#%%

Ntime=50
MF_list=pd.read_csv('gas_advection_profile_01_04_MASS_FLUX_TOTAL_Z.csv', header= None)
MF = MF_list.iloc[Ntime, 1:].values



#%%

PSI = 0.6 
v=-MF[3]/rhog[3]*10**(-3)/PSI

L=0.1 #m
Tc=L/v
print("v=",v," T=",Tc)


Yini_all = [0,1]
Ylim_all = [0.22,0.78]
NtimeL = [1,10,20,30,40, 50]

depths = YJ1l.iloc[0, 1:].astype(float).values

colors = cm.autumn(np.linspace(0, 0.8, len(NtimeL)))

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()
# faux tracés utilisés uniquement pour la légende
plt.plot([], [], color="black", linewidth=2.5, label="Gpyro")
plt.plot([], [],        marker=' ', markersize=15, linestyle='--', lw=2,
        markevery= 5, markeredgewidth=3, c="black",
         label="Analytical")


error02=[]
# -------- Tracé des courbes --------
for j, Ntime in enumerate(NtimeL):

    Time = YJ1l.iloc[Ntime, 0]
    YJ1 = YJ1l.iloc[Ntime, 1:].values
    YJ1 = np.array([fix_fortran_number(x) for x in YJ1])


    i=0
    Y0 = Yini_all[i]
    Yinf = Ylim_all[i]

    # solution analytique pure advection
    Yj_an = Yj_advection_analytical(depths, Time, v, Y0, Yinf)


    color = colors[j]

    #Analytical (marker only)
    plt.plot(
        100*depths, Yj_an,
        marker=' ', markersize=15, linestyle='--',
        markevery= 20, markeredgewidth=3,c=color,lw=3,
        markeredgecolor=color,
    )

    # Gpyro (line only)
    plt.plot(
        100*depths, YJ1,
        color=color, linewidth=2.5,
    )

    # Legende pour le temps (couleur uniquement)
    plt.plot([], [], color=color, linewidth=4, label=f"{Time:.1f} s")
    error=L2(depths, YJ1, depths, Yj_an, num_points=1000)
    error02.append(error)

# Flèche indiquant le sens du flux
ax = plt.gca()
ax.arrow(95*depths[-1], (Yinf+Y0)/2, -5, 0, head_width=0.02, head_length=0.5,lw=20, ec='black')
plt.text(80*depths[-1], (Yinf+Y0)*4/10, 'Flux direction', ha='center', va='center', color='black')

plt.xlabel('Depth [cm]')
plt.ylabel('Concentration [kg·kg⁻¹]')
plt.title("O₂ concentration")

#plt.ylim([0,1.2])
    
plt.legend(frameon=False, ncol=1)
plt.tight_layout()

plt.savefig(results_file_path1)



#%%

colors = cm.winter(np.linspace(0, 0.8, len(NtimeL)))

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()

errorN2=[]

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()
# faux tracés utilisés uniquement pour la légende
plt.plot([], [], color="black", linewidth=2.5, label="Gpyro")
plt.plot([], [],        marker=' ', markersize=15, linestyle='--', lw=2,
        markevery= 5, markeredgewidth=3, c="black",
         label="Analytical")

# -------- Tracé des courbes --------
for j, Ntime in enumerate(NtimeL):

    Time = YJ2l.iloc[Ntime, 0]
    YJ2 = YJ2l.iloc[Ntime, 1:].values
    YJ2 = np.array([fix_fortran_number(x) for x in YJ2])


    i=1
    Y0 = Yini_all[i]
    Yinf = Ylim_all[i]

    # solution analytique pure advection
    Yj_an = Yj_advection_analytical(depths, Time, v, Y0, Yinf)


    color = colors[j]

    #Analytical (marker only)
    plt.plot(
        100*depths, Yj_an,
        marker=' ', markersize=15, linestyle='--',
        markevery= 20, markeredgewidth=3,c=color,lw=3,
        markeredgecolor=color,
    )

    # Gpyro (line only)
    plt.plot(
        100*depths, YJ2,
        color=color, linewidth=2.5,
    )

    # Legende pour le temps (couleur uniquement)
    plt.plot([], [], color=color, linewidth=4, label=f"{Time:.1f} s")
    error=L2(depths, YJ1, depths, Yj_an, num_points=1000)
    error02.append(error)

# Flèche indiquant le sens du flux
ax = plt.gca()
ax.arrow(95*depths[-1], (Yinf+Y0)/2, -5, 0, head_width=0.02, head_length=0.5,lw=20, ec='black')
plt.text(80*depths[-1], (Yinf+Y0)*4.85/10, 'Flux direction', ha='center', va='center', color='black')

plt.xlabel('Depth [cm]')
plt.ylabel('Concentration [kg·kg⁻¹]')
plt.title("N₂ concentration")
    
plt.legend(frameon=False, ncol=1)
plt.tight_layout()
plt.savefig(results_file_path2)
#%%

front_pos=[]
time_list=[]

i=0
Y0 = Yini_all[i]
Yinf = Ylim_all[i]
for Ntime in range(len(YJ1l)):
    Time = YJ1l.iloc[Ntime, 0]
    YJ1 = YJ1l.iloc[Ntime, 1:].values
    YJ1 = np.array([fix_fortran_number(x) for x in YJ1])

    # --- position du front numérique ---
    zf = front_position(depths, YJ1, Y0, Yinf)
    front_pos.append(zf)
    time_list.append(Time)

# --- vitesse numérique ---
t_mid, v_num = front_velocity(time_list, front_pos)
mean_vnum=-np.mean(v_num[10:60])


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



# Norme Linfini : erreur maximale
def Linf(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    return np.max(np.abs(sim_interp - ref_interp))






error_value=abs(mean_vnum-v)/abs(v)
print(f"Abs error on the bvitesse of advection of Yj={np.round(error_value,6)}")
# Validation check
threshold = 0.002
if (error_value <= threshold)  :
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)