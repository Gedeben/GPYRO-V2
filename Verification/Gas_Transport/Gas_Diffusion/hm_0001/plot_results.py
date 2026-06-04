import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys
import re
from scipy.special import erfc
from matplotlib import cm



# If running as a script, SCRIPT_DIR will be the script's directory.
# For interactive environments, we'll assume the data file is in the current working directory.
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
except NameError:
    SCRIPT_DIR = os.getcwd()


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))


results_file_path1 = os.path.join(SCRIPT_DIR, 'O2_Diffusion_hm0001.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'N2_Diffusion_hm0001.png')


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 24  
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


# plotting the analytical Robin (mixed) semi-infinite solution for diffusion
def Yj_anlaytical(x, t, D, alpha,Ylim,Yini):
    z = x / (2.0 * np.sqrt(D * t))
    a = alpha * x + (alpha**2) * D * t
    arg = x / (2.0 * np.sqrt(D * t)) + alpha * np.sqrt(D * t)
    log_term = np.log(erfc(arg))
    term = np.exp(a + log_term)
    return Yini+(Ylim-Yini)*(erfc(z) - term)


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
#%%
# Load data
try:
    YJ1l=pd.read_csv("gas_diffusion_profile_01_01_YJ(01).csv", header= None)
    YJ2l=pd.read_csv("gas_diffusion_profile_01_02_YJ(02).csv", header= None)
    Gpyro=pd.read_csv("gas_diffusion_summary_01_0001.csv")
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
M=32*0.001
R=8.314

M1=32*1
M2=28.97
sigma1=5.061
sigma2=5.061
EPSOK1=254
EPSOK2=254

M1=32*1
M2=32
sigma1=3.433
sigma2=3.681
EPSOK1=113
EPSOK2=91.5

sigma=0.5*(sigma1+sigma2)
EPSOK=(EPSOK1*EPSOK2)**0.5
Te=T0/EPSOK

C1=1/M1 +1/M2




Gamma=a1*Te**a2+a3*np.exp(a4*Te)+a5*np.exp(a6*Te)
alpha1= 0.018829*(T0**3*C1)**(1/2)/(sigma**2*Gamma)


D=alpha1/Pref
#D=2.07016517e-05
#%%


hm  = 0.0001
K = 5E-10
PSI = 0.6

Yini_all = [0,1]
Ylim_all = [0.22,0.78]
M_all = np.array([32,32])*0.001
rhog_all = Pref*M_all/(R*Tref)
alpha_all = hm/(PSI*rhog_all*D)
NtimeL = [21, 401, 1001, 1801]
depths = YJ1l.iloc[0, 1:].astype(float).values

colors = cm.autumn(np.linspace(0, 0.8, len(NtimeL)))

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()
# faux tracés utilisés uniquement pour la légende
plt.plot([], [], color="black", linewidth=2.5, label="Gpyro (line)")
plt.plot([], [],        marker='+', markersize=15, linestyle='',
        markevery= 5, markeredgewidth=3, markeredgecolor="black",
         label="Analytical (marker)")


error02=[]
# -------- Tracé des courbes --------
for j, Ntime in enumerate(NtimeL):

    Time = YJ1l.iloc[Ntime, 0]
    YJ1 = YJ1l.iloc[Ntime, 1:].values
    YJ1 = np.array([fix_fortran_number(x) for x in YJ1])


    i=0
    Y0 = Yini_all[i]
    Ylim = Ylim_all[i]
    alpha = alpha_all[i]
    Yj_an02 = Yj_anlaytical(depths, Time, D, alpha, Ylim, Y0)

    color = colors[j]

    # Analytical (marker only)
    plt.plot(
        100*depths, Yj_an02,
        marker='+', markersize=15, linestyle='',
        markevery= 5, markeredgewidth=3,
        markeredgecolor=color,
    )

    # Gpyro (line only)
    plt.plot(
        100*depths, YJ1,
        color=color, linewidth=2.5,
    )

    # Legende pour le temps (couleur uniquement)
    plt.plot([], [], color=color, linewidth=4, label=f"{Time:.1f} s")
    error=L2(depths, YJ1, depths, Yj_an02, num_points=1000)
    error02.append(error)

plt.xlabel('Depth [cm]')
plt.ylabel('Concentration [kg·kg⁻¹]')
plt.title("O₂ concentration")
    
plt.legend(frameon=False, fontsize=26, ncol=2)
plt.tight_layout()

plt.savefig(results_file_path1)



#%%

colors = cm.winter(np.linspace(0, 0.8, len(NtimeL)))

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()
# faux tracés utilisés uniquement pour la légende
plt.plot([], [], color="black", linewidth=2.5, label="Gpyro (line)")
plt.plot([], [],        marker='+', markersize=15, linestyle='',
        markevery= 5, markeredgewidth=3, markeredgecolor="black",
         label="Analytical (marker)")

errorN2=[]
# -------- Tracé des courbes --------
for j, Ntime in enumerate(NtimeL):

    Time = YJ2l.iloc[Ntime, 0]
    YJ2 = YJ2l.iloc[Ntime, 1:].values
    YJ2 = np.array([fix_fortran_number(x) for x in YJ2])


    i=1
    Y0 = Yini_all[i]
    Ylim = Ylim_all[i]
    alpha = alpha_all[i]
    Yj_an = Yj_anlaytical(depths, Time, D, alpha, Ylim, Y0)

    color = colors[j]

    # Analytical (marker only)
    plt.plot(
        100*depths, Yj_an,
        marker='+', markersize=15, linestyle='',
        markevery= 5, markeredgewidth=3,
        markeredgecolor=color,
    )

    # Gpyro (line only)
    plt.plot(
        100*depths, YJ2,
        color=color, linewidth=2.5,
    )

    # Legende pour le temps (couleur uniquement)
    plt.plot([], [], color=color, linewidth=4, label=f"{Time:.1f} s")

    error=L2(depths, YJ1, depths, Yj_an02, num_points=1000)
    errorN2.append(error)

plt.xlabel('Depth [cm]')
plt.ylabel('Concentration [kg·kg⁻¹]')
plt.title("N₂ concentration")
    
plt.legend(frameon=False, fontsize=26, ncol=2)
plt.tight_layout()
plt.savefig(results_file_path2)



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






error_value=max(np.max(error02),np.max(errorN2))
print(f"Max L2 error on YJ={np.round(error_value,5)}")
# Validation check
threshold = 0.002
if (error_value <= threshold)  :
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)