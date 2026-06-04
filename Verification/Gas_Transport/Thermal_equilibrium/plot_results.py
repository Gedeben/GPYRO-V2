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



style_traces = {
    "Gpyro": {"linewidth": 3},
    "Analytical": {'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 3, 'markeredgewidth': 4},
}

case_styles = [
{'color': 'blue'},
{'color': 'black'},
{'color': 'red'}
]

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
    T= pd.read_csv("Thermal_equilibrium_profile_01_03_TEMPERATURE.csv", header= None)
    T2= pd.read_csv("Thermal_equilibrium_profile_02_03_TEMPERATURE.csv", header= None)
    dQSC= pd.read_csv("Thermal_equilibrium_profile_01_02_QSG.csv", header= None)

except FileNotFoundError as e:
    print("Error: The data file was not found ")
    sys.exit(1)

# Time vector
results_file_path1 = os.path.join(SCRIPT_DIR, 'thermal_exchange_temperature.png')
results_file_path2 = os.path.join(SCRIPT_DIR, 'thermal_exchange_QSG.png')

#%%


Ntime=-1

depths = T.iloc[0, 1:].astype(float).values
Time= T.iloc[Ntime, 0]

Temp=T.iloc[Ntime, 1:].values
Temp2=T2.iloc[Ntime, 1:].values

plt.figure()
plt.plot(100*depths,Temp2, **style_traces["Gpyro"], **case_styles[1], label="Gpyro without thermal-solid\n heat exchange")
plt.plot(100*depths,Temp, **style_traces["Gpyro"], **case_styles[0], label="Gpyro with thermal-solid\n heat exchange")

plt.xlabel('depth [cm]')
plt.ylabel('Temperature [°C]')
plt.legend( frameon=False)
plt.tight_layout()
plt.savefig(results_file_path1)



#%%


CPG=1000
MFlux=-100
TinG=300-273.15

QSG_lim= CPG*MFlux*(1/1000)* (Temp[-1]-TinG)/(depths[-2]-depths[-1])

QSC_an= np.array(CPG*MFlux*(1/1000)* (Temp[1:]-Temp[:-1])/(depths[1:]-depths[:-1]))
depths = dQSC.iloc[0, 1:].astype(float).values
Time= dQSC.iloc[Ntime, 0]

QSC_an = np.append(QSC_an, QSG_lim)

QSC=dQSC.iloc[Ntime, 1:].values


plt.figure()

# Traces
plt.plot(100*depths[:], QSC, label="Gpyro", **style_traces["Gpyro"], **case_styles[0])
plt.plot(100*depths[:], QSC_an, label="Analytical", **style_traces["Analytical"], **case_styles[0])

# Axes et titre
plt.xlabel('Depth [mm]')
plt.ylabel('Heat Source [W·m⁻³]')  # ancien: Heat Flux [W.m-3]
plt.legend(frameon=False)

x_last = 100 * depths[-1]
y_last = QSC[-1]            

plt.annotate(
    "Impact of the gas entering \nat a different temperature \nthan the solid.",
    xy=(x_last, y_last),
    xytext=(x_last - 0.5*(100*(depths.max()-depths.min())), y_last + 0.2*(QSC.max()-QSC.min())),
    color='red',
    fontsize=24,
    arrowprops=dict(
        arrowstyle="->",
        color="red",
        lw=1.5,
        shrinkA=0, shrinkB=4
    ))

plt.tight_layout()
plt.savefig(results_file_path2)





#%%



error_value=np.max(QSC_an-QSC)/QSG_lim

error2=np.min(Temp2-Temp)
# Validation check
threshold = 0.005
if error_value <= threshold and error2 >=150:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)