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


# File paths
file_path = os.path.join(SCRIPT_DIR, 'Mixture_Law_summary_01_0001.csv')
# Load data
try:
    data = pd.read_csv(file_path)
except FileNotFoundError as e:
    print(f"Error: The data file was not found at '{file_path}'. {e}")
    sys.exit(1)


    
compare_to_old_gpyro = False

# Execute only if the specific argument is provided
if len(sys.argv) > 1 and sys.argv[1] == "compare_old_gpyro":
    Gpyro_old_path1 = os.path.join(SCRIPT_DIR, 'Gpyro_V08_Mixture_Law_summary_01_0001.csv')
    data_old = pd.read_csv(Gpyro_old_path1)
    compare_to_old_gpyro= True

# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 25
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['legend.fontsize'] = 22
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['figure.figsize'] = (25, 15) # Adjusted for single plots


#%%
r=0.001
K_A = 1
R_A = 1000
RS_A = 2000
Cp_A = 1200
EMIS_A = 1
PSI_A=(1-R_A/RS_A)

K_B = 0.5
R_B = 600
RS_B = 800
Cp_B = 1900
EMIS_B = 0
PSI_B=(1-R_B/RS_B)


# Analytical solutions

t = data['t'].values

Y_A = np.maximum(1 - r * t, 0)
Y_B = 1 - Y_A


RHO_BULK = 1/(Y_A/R_A + Y_B/R_B)

X_A= RHO_BULK*(Y_A/R_A)
X_B= RHO_BULK*(Y_B/R_B)

RHO_BULK = X_A*R_A + X_B*R_B
RHO_SOLID = 1/(Y_A/RS_A + Y_B/RS_B)
#RHO_SOLID = X_A*RS_A + X_B*RS_B  #False
#RHO_SOLID = Y_A*RS_A + Y_B*RS_B  #False

POROSITY =X_A*PSI_A+X_B*PSI_B
CP = Y_A * Cp_A + Y_B * Cp_B
KZ = X_A * K_A + X_B * K_B


#%%
style_traces = {
    "Gpyro": {"linewidth": 3},
    "GpyroV08": {"linewidth": 5,"linestyle":'--'},
    "Analytical":{'linestyle': 'none', 'markersize': 15,'marker':'+', 'markevery': 10, 'markeredgewidth': 6}
}

case_styles = [
 {'color': 'blue'},
{'color': 'black'},
{'color': 'red'},
{'color': 'green'},
{'color': 'purple'}
]

# --------------------
# Combined subplot figure (2x3 grid)
# --------------------
fig, axes = plt.subplots(2, 3, sharex=True)
axes = axes.flatten()


if compare_to_old_gpyro:    
    ax=axes[0]
    ax.plot(t, data_old["002_YI( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[2],label=r"Gpyro V0.8 $Y_A$")
    ax.plot(t, data_old["003_YI( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[3],label=r"Gpyro V0.8 $Y_B$")
    ax=axes[1]
    ax.plot(t, data_old["010_XI( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[2],label=r"Gpyro V0.8 $X_A$")
    ax.plot(t, data_old["011_XI( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[3],label=r"Gpyro V0.8 $X_B$")
    ax=axes[2]
    ax.plot(t, data_old["008_BULK_DENSITY( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[2],label=r"Gpyro V0.8 $\bar{\rho}$")
    ax.plot(t, data_old["007_SOLID_DENSITY( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[3],label=r"Gpyro V0.8 $\bar{\rho}_s$")
    ax=axes[3]
    ax.plot(t, data_old["005_THERMAL_CONDUCTIVITY_Z( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[2],label=r"Gpyro V0.8 $\bar{k}$")
    ax=axes[4]
    ax.plot(t, data_old["006_SPECIFIC_HEAT_CAPACITY( 0.0000_ 0.0000_ 0.0000)"],
            **style_traces["GpyroV08"],**case_styles[2],label=r"Gpyro V0.8 $\bar{c}_p$")


ax=axes[0]
ax.plot(t, data["002_YI( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $Y_A$")
ax.plot(t, data["003_YI( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[1],label=r"Gpyro $Y_B$")
ax.plot(t,Y_A, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $Y_A$")
ax.plot(t,Y_B, **style_traces["Analytical"],**case_styles[1],label=r"Analytical $Y_B$")
ax.set_title("Species Mass Fraction")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Mass Fraction (kg/kg)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)


ax=axes[1]
ax.plot(t, data["010_XI( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $X_A$")
ax.plot(t, data["011_XI( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[1],label=r"Gpyro $X_B$")
ax.plot(t,X_A, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $X_A$")
ax.plot(t,X_B, **style_traces["Analytical"],**case_styles[1],label=r"Analytical $X_B$")
ax.set_title("Species Volumetric Fraction")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Volume Fraction (m³/m³)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)

ax=axes[2]
ax.plot(t, data["008_BULK_DENSITY( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $\bar{\rho}$")
ax.plot(t, data["007_SOLID_DENSITY( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[1],label=r"Gpyro $\bar{\rho}_s$")
ax.plot(t,RHO_BULK, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $\bar{\rho}$")
ax.plot(t,RHO_SOLID, **style_traces["Analytical"],**case_styles[1],label=r"Analytical $\bar{\rho}_s$")
ax.set_title("Densities")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Density (kg/m³)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)

ax=axes[3]
ax.plot(t, data["005_THERMAL_CONDUCTIVITY_Z( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $\bar{k}$")
ax.plot(t,KZ, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $\bar{k}$")
ax.set_title("Thermal conductivity")
ax.set_xlabel("Time (s)")
ax.set_ylabel("k (W.m⁻¹·K⁻¹)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)


ax=axes[4]
ax.plot(t, data["006_SPECIFIC_HEAT_CAPACITY( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $\bar{c}_p$")
ax.plot(t,CP, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $\bar{c}_p$")
ax.set_title("Heat capacity")
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"$c_p$ (J.kg⁻¹·K⁻¹)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)

ax=axes[5]
ax.plot(t, data["009_POROSITY( 0.0000_ 0.0000_ 0.0000)"], **style_traces["Gpyro"],**case_styles[0],label=r"Gpyro $\bar{\psi}$")
ax.plot(t,POROSITY, **style_traces["Analytical"],**case_styles[0],label=r"Analytical $\bar{\psi}$")
ax.set_title("Porosity")
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"$\psi$ (-)")
ax.ticklabel_format(style="plain", axis="both")
ax.legend(frameon=False)
ax.grid(False)



fig.tight_layout()
combined_file = Path(SCRIPT_DIR) / "All_Variables_Subplots_ML.png"
plt.savefig(combined_file, dpi=300)

#%%


error1=np.max(np.abs(data["002_YI( 0.0000_ 0.0000_ 0.0000)"] - Y_A))
error2=np.max(np.abs(data["003_YI( 0.0000_ 0.0000_ 0.0000)"] - Y_B))
error3=np.max(np.abs(data["010_XI( 0.0000_ 0.0000_ 0.0000)"] - X_A))
error4=np.max(np.abs(data["011_XI( 0.0000_ 0.0000_ 0.0000)"] - X_B))
error5=np.max(np.abs(data["008_BULK_DENSITY( 0.0000_ 0.0000_ 0.0000)"] - RHO_BULK))/1000
error6=np.max(np.abs(data["007_SOLID_DENSITY( 0.0000_ 0.0000_ 0.0000)"] - RHO_SOLID))/1000
error7=np.max(np.abs(data["005_THERMAL_CONDUCTIVITY_Z( 0.0000_ 0.0000_ 0.0000)"] - KZ))
error8=np.max(np.abs(data["006_SPECIFIC_HEAT_CAPACITY( 0.0000_ 0.0000_ 0.0000)"] - CP))/1000
error9=np.max(np.abs(data["009_POROSITY( 0.0000_ 0.0000_ 0.0000)"] - POROSITY))

error=np.max([error1,error2,error3,error4,error5,error6,error7,error8,error9])
print(f"error for Y_A={error1} %")
print(f"error for Y_B={error2} %")
print(f"error for X_A={error3} %")
print(f"error for X_B={error4} %")
print(f"error for bulk density ={error5} %")
print(f"error for solid density={error6} %")
print(f"error for thermal conductivity={error7} %")
print(f"error for heat capacity={error8} %")
print(f"error for porosity={error9} %")

#%%


# Validation check
threshold = 0.00001
if error <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
