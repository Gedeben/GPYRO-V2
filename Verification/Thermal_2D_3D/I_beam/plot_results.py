import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path1 = os.path.join(SCRIPT_DIR, 'I_beam_point_T.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'I_beam_y_profile.png')
results_file_path3 = os.path.join(SCRIPT_DIR, 'I_beam_z_profile.png')


#%%


# Matplotlib styling
mpl.rcParams['font.family'] = 'serif'  
mpl.rcParams['font.size'] = 30  
mpl.rcParams['axes.labelsize'] = 30  
mpl.rcParams['legend.fontsize'] = 30  
mpl.rcParams['xtick.labelsize'] = 30 
mpl.rcParams['ytick.labelsize'] = 30  
mpl.rcParams['figure.figsize'] = (15, 8)


param= [
    {'color': 'grey', 'linewidth': 3},
    {'color': 'pink', 'linewidth': 3},
    {'color': 'darkred', 'linewidth': 3},
    {'color': 'turquoise', 'linewidth': 3},
    {'color': 'darkgreen', 'linewidth': 3},
    {'color': 'red', 'linewidth': 3},
    {'color': 'blue', 'linewidth': 3},
    {'color': 'black', 'linewidth': 3}
]



Gpyro_file="I_beam_summary_01_0001.csv"
data_Gpyro=pd.read_csv(Gpyro_file)

FDS_file="FDS/ht3d_ibeam_devc.csv"
FDS = pd.read_csv(FDS_file, skiprows=1)


G_t=data_Gpyro['t'].values
G_T1=data_Gpyro["001_TEMPERATURE( 0.0550_ 0.3950_ 0.3950)"].values
G_T2=data_Gpyro["002_TEMPERATURE( 0.0550_ 0.3450_ 0.2950)"].values
G_T3=data_Gpyro["003_TEMPERATURE( 0.0550_ 0.2950_ 0.1950)"].values
G_T4=data_Gpyro["004_TEMPERATURE( 0.0550_ 0.2250_ 0.3950)"].values
G_T5=data_Gpyro["005_TEMPERATURE( 0.0550_ 0.3950_ 0.0050)"].values
G_T6=data_Gpyro["006_TEMPERATURE( 0.0550_ 0.2250_ 0.0050)"].values


FDS_t=FDS["Time"].values
FDS_T1=FDS["TS_x195-40"].values
FDS_T2=FDS["TS_x145-30"].values
FDS_T3=FDS["TS_x095-20"].values
FDS_T4=FDS["TS_x025-40"].values
FDS_T5=FDS["TS_x195-01"].values
FDS_T6=FDS["TS_x025-01"].values


plt.figure()

plt.plot(G_t, G_T1,**param[0])
plt.plot(G_t, G_T2,**param[1])
plt.plot(G_t, G_T3,**param[2])
plt.plot(G_t, G_T4,**param[3])
plt.plot(G_t, G_T5,**param[4])
plt.plot(G_t, G_T6,**param[5])

plt.plot(FDS_t, FDS_T1, **param[0], linestyle='--', marker='+', markevery=10, markersize=15)
plt.plot(FDS_t, FDS_T2, **param[1], linestyle='--', marker='+', markevery=10, markersize=15)
plt.plot(FDS_t, FDS_T3, **param[2], linestyle='--', marker='+', markevery=10, markersize=15)
plt.plot(FDS_t, FDS_T4, **param[3], linestyle='--', marker='+', markevery=10, markersize=15)
plt.plot(FDS_t, FDS_T5, **param[4], linestyle='--', marker='+', markevery=10, markersize=15)
plt.plot(FDS_t, FDS_T6, **param[5], linestyle='--', marker='+', markevery=10, markersize=15)

plt.plot([], [], color='black', lw=3, linestyle='-', label="Gpyro")
plt.plot([], [], color='black', lw=3, linestyle='--',marker='+', markersize=15, label="FDS 3D")
plt.plot([], [], color='black', lw=3, linestyle=' ', label=" ")

plt.plot([], [], **param[0], label="Point 1")
plt.plot([], [], **param[1], label="Point 2")
plt.plot([], [], **param[2], label="Point 3")
plt.plot([], [], **param[3], label="Point 4")
plt.plot([], [], **param[4], label="Point 5")
plt.plot([], [], **param[5], label="Point 6")


# Ajout des labels et de la légende
plt.xlabel('Time [s]')
plt.ylabel('Temperature [°C]')
plt.legend(frameon=False, ncol=1,loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.tight_layout() 

plt.savefig(results_file_path1)


#%%


G_Tp1=pd.read_csv("I_beam_profile_01_01_TEMPERATURE.csv", header=None)
x_Gpyro = G_Tp1.iloc[0, 1:-1].tolist()


FDS_Tp1 = pd.read_csv("FDS/ht3d_ibeam_prof_1.csv", skiprows=2, header=None)


frac=[1,2,8]

plt.figure(figsize=(12, 8))



colors=['r', 'b', 'g', 'k']

i=0
for num in frac:
    indice=int((np.shape(G_Tp1)[0]-2)/num)+1
    Gpyro_time = G_Tp1.iloc[indice, 0]
    Gpyro_T1= G_Tp1.iloc[indice, 1:-1].tolist()
    
    #indice=np.shape(df)[0]-1
    indice=int((np.shape(FDS_Tp1)[0]-1)/num)-0
    FDS_Time=FDS_Tp1.iloc[indice, 0]
    n = int(FDS_Tp1.iloc[indice, 1])
    x_FDS = FDS_Tp1.iloc[indice,2:n+2].tolist()
    FDS_Tl1 = FDS_Tp1.iloc[indice, n+2:2*n+2].tolist()
    
    color=colors[i]
    i+=1    
    # Tracer la quantité en fonction de la profondeur pour chaque temps
    plt.plot(x_FDS, FDS_Tl1, color=color, label=f'FDS t={round(FDS_Time,)} s', lw=3)
    plt.plot(x_Gpyro, Gpyro_T1[::-1], color=color, linestyle='--', label=f'Gpyro t={round(Gpyro_time,)} s',lw=3)


plt.legend(frameon=False,fontsize=25)
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
#plt.title("Profile temperature in x direction (y,z=0.3,0.03)")
plt.tight_layout() 

plt.savefig(results_file_path2)


#%%



G_Tp2=pd.read_csv("I_beam_profile_01_02_TEMPERATURE.csv", header=None)
x_Gpyro = G_Tp2.iloc[0, 1:-1].tolist()


FDS_Tp2 = pd.read_csv("FDS/ht3d_ibeam_prof_2.csv", skiprows=2, header=None)


frac=[1,2,8]

plt.figure(figsize=(12, 8))



colors=['r', 'b', 'g', 'k']

i=0
for num in frac:
    indice=int((np.shape(G_Tp2)[0]-2)/num)+1
    Gpyro_time = G_Tp2.iloc[indice, 0]
    Gpyro_T2= G_Tp2.iloc[indice, 1:-1].tolist()
    
    indice=int((np.shape(FDS_Tp2)[0]-1)/num)-0
    FDS_Time=FDS_Tp2.iloc[indice, 0]
    n = int(FDS_Tp2.iloc[indice, 1])
    x_FDS = FDS_Tp2.iloc[indice,2:n+2].tolist()
    FDS_Tl2 = FDS_Tp2.iloc[indice, n+2:2*n+2].tolist()
    
    color=colors[i]
    i+=1    
    # Tracer la quantité en fonction de la profondeur pour chaque temps
    plt.plot(x_FDS, FDS_Tl2, color=color, label=f'FDS t={round(FDS_Time,)} s', lw=3)
    plt.plot(x_Gpyro, Gpyro_T2[::-1], color=color, linestyle='--', label=f'Gpyro t={round(Gpyro_time,)} s',lw=3)


plt.legend(frameon=False,fontsize=25)
plt.xlabel("y position [m]")
plt.ylabel("Temperature [°C]")
#plt.title(""Profile temperature in y direction (x,z=0.3,0.03)")
plt.tight_layout() 

#%%


G_Tp2=pd.read_csv("I_beam_profile_01_03_TEMPERATURE.csv", header=None)
x_Gpyro = G_Tp2.iloc[0, 1:-1].tolist()


FDS_Tp2 = pd.read_csv("FDS/ht3d_ibeam_prof_3.csv", skiprows=2, header=None)


frac=[1,2,8]

plt.figure(figsize=(12, 8))



colors=['r', 'b', 'g', 'k']

i=0
for num in frac:
    indice=int((np.shape(G_Tp2)[0]-2)/num)+1
    Gpyro_time = G_Tp2.iloc[indice, 0]
    Gpyro_T2= G_Tp2.iloc[indice, 1:-1].tolist()
    
    indice=int((np.shape(FDS_Tp2)[0]-1)/num)-0
    FDS_Time=FDS_Tp2.iloc[indice, 0]
    n = int(FDS_Tp2.iloc[indice, 1])
    x_FDS = FDS_Tp2.iloc[indice,2:n+2].tolist()
    FDS_Tl2 = FDS_Tp2.iloc[indice, n+2:2*n+2].tolist()
    
    color=colors[i]
    i+=1    
    # Tracer la quantité en fonction de la profondeur pour chaque temps
    plt.plot(x_FDS, FDS_Tl2, color=color, label=f'FDS t={round(FDS_Time,)} s', lw=3)
    plt.plot(x_Gpyro, Gpyro_T2[::-1], color=color, linestyle='--', label=f'Gpyro t={round(Gpyro_time,)} s',lw=3)


plt.legend(frameon=False,fontsize=25)
plt.xlabel("z position [m]")
plt.ylabel("Temperature [°C]")
#plt.title("Profile temperature in z direction (x,y=0.2,0.3)")
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



# Norme Linfini : erreur maximale
def Linf(sim_x, sim_y, ref_x, ref_y, num_points=1000):
    min_x = max(min(sim_x), min(ref_x))
    max_x = min(max(sim_x), max(ref_x))
    uniform_x = np.linspace(min_x, max_x, num_points)
    
    sim_interp = np.interp(uniform_x, sim_x, sim_y)
    ref_interp = np.interp(uniform_x, ref_x, ref_y)
    
    return np.max(np.abs(sim_interp - ref_interp))



err1=L1(G_t, G_T1, FDS_t, FDS_T1)
err2=L1(G_t, G_T2, FDS_t, FDS_T2)
err3=L1(G_t, G_T3, FDS_t, FDS_T3)
err4=L1(G_t, G_T4, FDS_t, FDS_T4)
err5=L1(G_t, G_T5, FDS_t, FDS_T5)
err6=L1(G_t, G_T6, FDS_t, FDS_T6)

err=[err1, err2, err3, err4, err5, err6]
Error= np.max(err)

#%%
G_t=data_Gpyro['t'].values
G_T1=data_Gpyro["001_TEMPERATURE( 0.0550_ 0.3950_ 0.3950)"].values
G_T2=data_Gpyro["002_TEMPERATURE( 0.0550_ 0.3450_ 0.2950)"].values
G_T3=data_Gpyro["003_TEMPERATURE( 0.0550_ 0.2950_ 0.1950)"].values
G_T4=data_Gpyro["004_TEMPERATURE( 0.0550_ 0.2250_ 0.3950)"].values
G_T5=data_Gpyro["005_TEMPERATURE( 0.0550_ 0.3950_ 0.0050)"].values
G_T6=data_Gpyro["006_TEMPERATURE( 0.0550_ 0.2250_ 0.0050)"].values


FDS_t=FDS["Time"].values
FDS_T1=FDS["TS_x195-40"].values
FDS_T2=FDS["TS_x145-30"].values
FDS_T3=FDS["TS_x095-20"].values
FDS_T4=FDS["TS_x025-40"].values
FDS_T5=FDS["TS_x195-01"].values
FDS_T6=FDS["TS_x025-01"].values




print(f"Abs error on the point temperature ={np.round(Error,1)} °C")
# Validation check
threshold = 20
if (Error <= threshold)  :
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
