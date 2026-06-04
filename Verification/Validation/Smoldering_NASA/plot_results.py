import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys
from matplotlib import cm
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))


results_file_path1 = os.path.join(SCRIPT_DIR, 'Thermocouple_Smoldering.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'MLR_Smoldering.png')
results_file_path3 = os.path.join(SCRIPT_DIR, 'Smoldering_front_propagation.png')

file_path = os.path.join(SCRIPT_DIR, '1d_smouldering_summary_01_0001.csv')

#%% color and style parameters for the plots in black and white
param_sets = [
    {'color': 'blue' , 'linewidth': 2, 'linestyle': '-'},
    {'color': 'red', 'linewidth': 2, 'linestyle': '-'},
    {'color': 'green' , 'linewidth': 2, 'linestyle': '-'},
    {'color': 'brown' , 'linewidth': 2, 'linestyle': '-'},
    {'color': 'black' , 'linewidth': 2, 'linestyle': '-'},
    {'color': 'pink', 'linewidth': 2, 'linestyle': '-'},
    {'color': 'green' , 'linewidth': 2, 'linestyle': '', 'marker': '^', 'markersize': 3, 'markevery': 100},
    {'color': 'orange' , 'linewidth': 2, 'linestyle': '', 'marker': 'd', 'markersize': 3, 'markevery': 40}
 ]

#%%

# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 25
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['legend.fontsize'] = 25
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['figure.figsize'] = (12, 8)

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
    data = pd.read_csv(file_path)
    # analytical_data = pd.read_csv(analytical_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)


indice=1
Time1 = data["t"].values[indice:]
MLR_TOTAL=data['001_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]
MLR_1=data['002_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]
MLR_2=data['003_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]
MLR_3=data['004_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]
MLR_4=data['005_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]
MLR_5=data['006_MLR( 0.0000_ 0.0000_ 0.0000)'].values[indice:]


plt.figure()
plt.plot(Time1, MLR_TOTAL, label='All gases', **param_sets[0])
plt.plot(Time1, MLR_1, label='T-pyrolysate', **param_sets[1])
plt.plot(Time1, MLR_4, label=r'O$_2$-pyrolysate', **param_sets[4])
plt.plot(Time1, MLR_5, label='product', **param_sets[5])
plt.plot(Time1, MLR_2, label=r'O$_2$', **param_sets[2])
plt.plot(Time1, MLR_3, label=r'N$_2$', **param_sets[3])


plt.xlabel('Time (s)')
plt.ylabel('Net Mass Flux '+r"[$\mathrm{g{\cdot}m^{-2}{\cdot}s^{-1}}$]")

plt.xlim([0,1100])
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(results_file_path2)


#%%




T_1=data['009_TEMPERATURE( 0.0050_ 0.0000_ 0.0000)'].values[indice:]
T_2=data['010_TEMPERATURE( 0.0250_ 0.0000_ 0.0000)'].values[indice:]
T_3=data['011_TEMPERATURE( 0.0450_ 0.0000_ 0.0000)'].values[indice:]
T_4=data['012_TEMPERATURE( 0.0650_ 0.0000_ 0.0000)'].values[indice:]
T_5=data['013_TEMPERATURE( 0.0850_ 0.0000_ 0.0000)'].values[indice:]
T_6=data['014_TEMPERATURE( 0.1050_ 0.0000_ 0.0000)'].values[indice:]
T_7=data['015_TEMPERATURE( 0.1250_ 0.0000_ 0.0000)'].values[indice:]
T_8=data['016_TEMPERATURE( 0.1400_ 0.0000_ 0.0000)'].values[indice:]


# Load experimental data data
try:
    exp_data = pd.read_csv("Experimental_data_smoldering.csv")
    exp_find = True
except FileNotFoundError :
    print("Error: The experimental files data was not found.")
    exp_find = False

# Définir les styles de tracé pour les données expérimentales et de simulation
param_sets_simulation = [
    {'color': 'grey', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'pink', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'darkred', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'turquoise', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'darkgreen', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'red', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'blue', 'linewidth': 3, 'linestyle': '-'},
    {'color': 'black', 'linewidth': 3, 'linestyle': '-'}
]


param_sets_experiment = [
    {'color': 'grey', 'marker': 's', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'pink', 'marker': 'D', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'darkred', 'marker': '^', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'turquoise', 'marker': 'x', 'markersize': 10, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'darkgreen', 'marker': 'v', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'red', 'marker': 's', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'blue', 'marker': 'D', 'markersize': 8, 'linestyle': ' ', 'markerfacecolor': 'none'},
    {'color': 'black', 'marker': '+', 'markersize': 10, 'linestyle': ' ', 'markerfacecolor': 'none'}
]

plt.figure(figsize=(16, 8))
plt.plot(Time1, T_1, label='  5 mm - Gpyro', **param_sets_simulation[0])
plt.plot(Time1, T_2, label=' 25 mm - Gpyro', **param_sets_simulation[1])
plt.plot(Time1, T_3, label=' 45 mm - Gpyro', **param_sets_simulation[2])
plt.plot(Time1, T_4, label=' 65 mm - Gpyro', **param_sets_simulation[3])
plt.plot(Time1, T_5, label=' 85 mm - Gpyro', **param_sets_simulation[4])
plt.plot(Time1, T_6, label='105 mm - Gpyro', **param_sets_simulation[5])
plt.plot(Time1, T_7, label='125 mm - Gpyro', **param_sets_simulation[6])
plt.plot(Time1, T_8, label='140 mm - Gpyro', **param_sets_simulation[7])

if exp_find:
    plt.plot(exp_data['Time'],exp_data['0_5cm'], label='  5 mm- Exp', **param_sets_experiment[0])
    plt.plot(exp_data['Time'],exp_data['2_5cm'], label=' 25 mm- Exp', **param_sets_experiment[1])
    plt.plot(exp_data['Time'],exp_data['4_5cm'], label=' 45 mm- Exp', **param_sets_experiment[2])
    plt.plot(exp_data['Time'],exp_data['6_5cm'], label=' 65 mm- Exp', **param_sets_experiment[3])
    plt.plot(exp_data['Time'],exp_data['8_5cm'], label=' 85 mm- Exp', **param_sets_experiment[4])
    plt.plot(exp_data['Time'],exp_data['10_5cm'], label='105 mm- Exp', **param_sets_experiment[5])
    plt.plot(exp_data['Time'],exp_data['12_5cm'], label='125 mm- Exp', **param_sets_experiment[6])
    plt.plot(exp_data['Time'],exp_data['14cm'], label='140 mm- Exp', **param_sets_experiment[7])


plt.xlabel('Time (s)')
plt.ylabel('Temperature [°C]')
plt.xlim([0,1000])
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False, ncol=1,fontsize=20 )
plt.tight_layout()

plt.savefig(results_file_path1)



#%%

YI1_list=pd.read_csv('1d_smouldering_profile_01_02_YI(01).csv', header= None)
YI2_list=pd.read_csv('1d_smouldering_profile_01_03_YI(02).csv', header= None)
YI3_list=pd.read_csv('1d_smouldering_profile_01_04_YI(03).csv', header= None)
YI4_list=pd.read_csv('1d_smouldering_profile_01_05_YI(04).csv', header= None)
YI5_list=pd.read_csv('1d_smouldering_profile_01_06_YI(05).csv', header= None)

depths = YI1_list.iloc[0, 1:].astype(float).values

NtimeL = [401,601,801,1001]# 400, 600, 800]

for Ntime in NtimeL:
        
    Time = YI1_list.iloc[Ntime, 0]
    
    
    YI1=YI1_list.iloc[Ntime, 1:].values
    YI2=YI2_list.iloc[Ntime, 1:].values
    YI3=YI3_list.iloc[Ntime, 1:].values
    YI4=YI4_list.iloc[Ntime, 1:].values
    YI5=YI5_list.iloc[Ntime, 1:].values
    
    
    plt.figure()
    plt.plot(100*depths, YI1, label='foam' , **param_sets[0])
    plt.plot(100*depths, YI2, label=r'$\beta$-foam', **param_sets[1])
    plt.plot(100*depths, YI3, label='T-char', **param_sets[2])
    plt.plot(100*depths, YI4, label='char', **param_sets[3])
    plt.plot(100*depths, YI5, label=r'$\alpha$-char', **param_sets[4])
    
    plt.xlabel('Depth [cm]')
    plt.ylabel('Mass fraction [kg/kg]')
    plt.legend(frameon=False)
    plt.xlim([0,14])
    plt.tight_layout()
    plt.savefig(os.path.join(SCRIPT_DIR, f"Smoldering_sample_composition_t_{Time:.0f}.png"))


#%%

YJ1_list=pd.read_csv('1d_smouldering_profile_01_07_YJ(01).csv', header= None)
YJ2_list=pd.read_csv('1d_smouldering_profile_01_08_YJ(02).csv', header= None)
YJ3_list=pd.read_csv('1d_smouldering_profile_01_09_YJ(03).csv', header= None)
YJ4_list=pd.read_csv('1d_smouldering_profile_01_10_YJ(04).csv', header= None)
YJ5_list=pd.read_csv('1d_smouldering_profile_01_11_YJ(05).csv', header= None)


for Ntime in NtimeL:
            
    Time = YI1_list.iloc[Ntime, 0]
    
    
    YJ1=YJ1_list.iloc[Ntime, 1:].values
    YJ2=YJ2_list.iloc[Ntime, 1:].values
    YJ3=YJ3_list.iloc[Ntime, 1:].values
    YJ4=YJ4_list.iloc[Ntime, 1:].values
    YJ5=YJ5_list.iloc[Ntime, 1:].values
    
    YJ1 = np.array([fix_fortran_number(x) for x in YJ1])
    YJ2 = np.array([fix_fortran_number(x) for x in YJ2])
    YJ3 = np.array([fix_fortran_number(x) for x in YJ3])
    YJ4 = np.array([fix_fortran_number(x) for x in YJ4])
    YJ5 = np.array([fix_fortran_number(x) for x in YJ5])
    
    YJsum=YJ1+YJ2+YJ3+YJ4+YJ5
    
    
    plt.figure()
    plt.plot(100*depths, YJ2, label=r'O$_2$', **param_sets[1])
    plt.plot(100*depths, YJ3, label=r'N$_2$', **param_sets[2])
    plt.plot(100*depths, YJ1, label='T-pyrolysate' , **param_sets[0])
    plt.plot(100*depths, YJ4, label=r'O$_2$-pyrolysate', **param_sets[3])
    plt.plot(100*depths, YJ5, label='product', **param_sets[4])
    
    plt.xlabel('Depth [cm]')
    plt.ylabel('Gas mass fraction [kg/kg]')
    plt.legend(frameon=False)
    plt.ylim((-0.02,1.02))
    plt.tight_layout()
    plt.savefig(os.path.join(SCRIPT_DIR, f"Smoldering_gas_mass_fraction_t_{Time:.0f}.png"))

#%%



NtimeL = [501,601,701,801,901]# 400, 600, 800]
colors1 = cm.autumn(np.linspace(0.8, 0, len(NtimeL)))
#colors2 = cm.winter(np.linspace(0.8, 0, len(NtimeL)))

# ---- Pour créer la légende "style" (Gpyro vs Analytical) ----
plt.figure()

# -------- Tracé des courbes --------
for j, Ntime in enumerate(NtimeL):

    Time = YI4_list.iloc[Ntime, 0]
    YI4 = YI4_list.iloc[Ntime, 1:].values
    #YI2 = YI2_list.iloc[Ntime, 1:].values

    color1 = colors1[j]
    #color2 = colors2[j]

    # Gpyro (line only)
    plt.plot(100*depths, YI4, color=color1, linewidth=2.5)
    #plt.plot(100*depths, YI2, color=color2, linewidth=2.5)


    # Legende pour le temps (couleur uniquement)
    plt.plot([], [], color=color1, linewidth=4, label=f"{Time:.0f} s")


# Flèche indiquant le sens du flux
ax = plt.gca()
ax.arrow(95*depths[-1], 0.4, -10, 0, head_width=0.03, head_length=0.2,lw=30, ec='black')
plt.text(60*depths[-1], 0.5, 'Smoldering front propagation', ha='center', va='center', color='black',fontsize=30)

plt.xlabel('Depth [cm]')
plt.ylabel('Char mass fraction [kg/kg]')
plt.xlim([0,14])
plt.ylim([-0.01,1])
plt.legend(frameon=False, ncol=3)
plt.tight_layout()
plt.savefig(results_file_path3)

#%%


t = YI4_list.iloc[:, 0].to_numpy() 
profiles = YI4_list.iloc[:, 1:].to_numpy()
x = np.asarray(depths)                      

def peak_pos_quadratic(x, y):
    """
    Retourne la position x_peak du maximum de y(x).
    1) prend l'indice du maximum global,
    2) si possible, ajuste une parabole locale sur (x[i-1], x[i], x[i+1]) pour sous-maille.
    Gère NaN et bords.
    """
    y = np.asarray(y)
    if np.all(~np.isfinite(y)):
        return np.nan

    i = np.nanargmax(y)
    # Si max en bord ou points voisins indisponibles, on prend la maille du max
    if i == 0 or i == len(y) - 1:
        return x[i]

    xs = np.array([x[i-1], x[i], x[i+1]], dtype=float)
    ys = np.array([y[i-1], y[i], y[i+1]], dtype=float)

    if np.any(~np.isfinite(ys)):
        return x[i]

    # Ajustement quadratique : y = a*x^2 + b*x + c
    a, b, c = np.polyfit(xs, ys, 2)
    if a == 0 or not np.isfinite(a):
        return x[i]

    x_peak = -b / (2*a)

    # Sécurité : ne pas extrapoler loin des 3 points
    if (x_peak < xs.min()) or (x_peak > xs.max()) or (not np.isfinite(x_peak)):
        return x[i]

    return x_peak

# ---------- Position du front pour chaque temps ----------
front_x = np.array([peak_pos_quadratic(x, profiles[k, :]) for k in range(profiles.shape[0])])

# ---------- Vitesse du front ----------
v = np.gradient(front_x[500:920], t[500:920])
V_front=-np.mean(v)*1000

print(f"Smoldering front speed ={V_front:.3f} mm/s")

