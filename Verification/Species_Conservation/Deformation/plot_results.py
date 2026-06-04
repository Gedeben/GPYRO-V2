import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from matplotlib.lines import Line2D
import os
import sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_path1 = os.path.join(SCRIPT_DIR, 'deformation_pure_condesed/deformation_pure_condensed_summary_01_0001.csv')
file_path2 = os.path.join(SCRIPT_DIR, 'no_deformation/no_deformation_summary_01_0001.csv')
file_path3 = os.path.join(SCRIPT_DIR, 'swelling/swelling_summary_01_0001.csv')
file_path4 = os.path.join(SCRIPT_DIR, 'shrinking/shrinking_summary_01_0001.csv')
file_path5 = os.path.join(SCRIPT_DIR, 'non_charring/non_charring_summary_01_0001.csv')


#%% color and style parameters for the plots in black and white

style_traces = {
    "Gpyro": {"linewidth": 3},
    "Analytical": {'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 1, 'markeredgewidth': 4},
}

case_styles = [
 {'color': 'blue'},
{'color': 'black'},
{'color': 'red'},
{'color': 'green'},
{'color': 'purple'}
]

# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 30
mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['legend.fontsize'] = 24
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['figure.figsize'] = (15,8)

#%%
# Load data
try:
    data1 = pd.read_csv(file_path1)
    data2 = pd.read_csv(file_path2)
    data3 = pd.read_csv(file_path3)
    data4 = pd.read_csv(file_path4)
    data5 = pd.read_csv(file_path5)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)


#%% Analytical

# Case 1: deformation_pure_condesed

rho_0 = 1000  # Reactant density (kg/m³)
rho_1 = 2000  # Product density (kg/m³)
chi=0
zdim_0 = 0.01 # Initial thickness (m)
rho_ratio=(rho_0 / rho_1)
theta = 1+((rho_1 / rho_0)-1)*chi
D=rho_ratio *theta
zdim_1 = zdim_0 *D  # Final thickness (m)


# Case 2: no_deformation

rho_0 = 1000  # Reactant density (kg/m³)
rho_1 = 630  # Product density (kg/m³)
chi=1
zdim_0 = 0.01 # Initial thickness (m)
rho_ratio=(rho_0 / rho_1)
theta = 1+((rho_1 / rho_0)-1)*chi

D=rho_ratio *theta
zdim_2 = zdim_0 *D  # Final thickness (m)


# Case 3: swelling
rho_0 = 1000  # Reactant density (kg/m³)
rho_1 = 400  # Product density (kg/m³)
chi=0.3333
zdim_0 = 0.01 # Initial thickness (m)
rho_ratio=(rho_0 / rho_1)
theta = 1+((rho_1 / rho_0)-1)*chi
D=rho_ratio *theta
zdim_3 = zdim_0 *D  # Final thickness (m)


# Case 4: shrinking

rho_0 = 1000  # Reactant density (kg/m³)
rho_1 = 600  # Product density (kg/m³)
chi=1.5
zdim_0 = 0.01 # Initial thickness (m)
rho_ratio=(rho_0 / rho_1)
theta = 1+((rho_1 / rho_0)-1)*chi
D=rho_ratio *theta
zdim_4 = zdim_0 *D  # Final thickness (m)

# Case 5: non_charring
zdim_5=0


#%%

# Time vector
t1 = data1['t'].values
t2 = data2['t'].values
t3 = data3['t'].values
t4 = data4['t'].values
t5 = data5['t'].values

zdim_gpyro1 = data1['002_THICKNESS( 0.0000_ 0.0000_ 0.0000)'].values  # in m
zdim_gpyro2 = data2['002_THICKNESS( 0.0000_ 0.0000_ 0.0000)'].values  # in m
zdim_gpyro3 = data3['002_THICKNESS( 0.0000_ 0.0000_ 0.0000)'].values  # in m
zdim_gpyro4 = data4['002_THICKNESS( 0.0000_ 0.0000_ 0.0000)'].values  # in m
zdim_gpyro5 = data5['002_THICKNESS( 0.0000_ 0.0000_ 0.0000)'].values  # in m


#%%


# ================== FIGURE & AXES ==================
fig = plt.figure()

# Axe principal (gauche)
ax = fig.add_axes([0.08, 0.12, 0.7, 0.8])  # [left, bottom, width, height]

# Axe fantôme pour les légendes (droite)
ax_leg = fig.add_axes([0.78, 0.12, 0.30, 0.8])
ax_leg.axis("off")


n=9
tp =[t1[int(np.size(zdim_gpyro1)*i/20)-1]/60 for i in range(12,21)]

# ================== PLOTS ==================
ax.plot(t1/60, zdim_gpyro1/zdim_0, **style_traces["Gpyro"], **case_styles[0])
ax.plot(t2/60, zdim_gpyro2/zdim_0, **style_traces["Gpyro"], **case_styles[1])
ax.plot(t3/60, zdim_gpyro3/zdim_0, **style_traces["Gpyro"], **case_styles[2])
ax.plot(t4/60, zdim_gpyro4/zdim_0, **style_traces["Gpyro"], **case_styles[3])
ax.plot(t5/60, zdim_gpyro5/zdim_0, **style_traces["Gpyro"], **case_styles[4])


ax.plot(tp, zdim_1/zdim_0*np.ones(n), **style_traces["Analytical"], **case_styles[0])
ax.plot(tp, zdim_2/zdim_0*np.ones(n), **style_traces["Analytical"], **case_styles[1])
ax.plot(tp, zdim_3/zdim_0*np.ones(n), **style_traces["Analytical"], **case_styles[2])
ax.plot(tp, zdim_4/zdim_0*np.ones(n), **style_traces["Analytical"], **case_styles[3])
ax.plot(tp, zdim_5/zdim_0*np.ones(n), **style_traces["Analytical"], **case_styles[4])

ax.set_xlabel("Time (min)")
ax.set_ylabel("Deforamation (cm/cm)")
ax.set_xlim(0, t1[-1]/60)
ax.grid(False)

# ================== LEGENDES ==================
legend_codes = [
    Line2D([0], [0], **style_traces["Gpyro"], label="Gpyro", color="gray"),
    Line2D([0], [0], **style_traces["Analytical"], label="Analytical \nfinal thickness", color="gray"),
]

case_labels = [
    "Pure solid \nphase reaction",
    "No deformation",
    "Swelling",
    "Shrinking",
    "Non-charring",
]

legend_cases = [
    Line2D(
        [0], [0],
        lw=3,
        label=case_labels[k],
        **case_styles[k],
        marker="s",
        markersize=12,
        linestyle="none"
    )
    for k in range(len(case_styles))
]

leg1 = ax_leg.legend(
    handles=legend_codes,
    loc="upper left",
    frameon=False,
)

leg2 = ax_leg.legend(
    handles=legend_cases,
    loc="lower left",
    frameon=False,
)

ax_leg.add_artist(leg1)

# ================== SAVE ==================
output_file = Path(SCRIPT_DIR) / "deformation.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')



#%%
# Function to compute absolute error between simulation and analytical solution

erreur1 = abs(zdim_gpyro1[-1] - zdim_1)/zdim_0
erreur2 = abs(zdim_gpyro2[-1] - zdim_2)/zdim_0
erreur3 = abs(zdim_gpyro3[-1] - zdim_3)/zdim_0
erreur4 = abs(zdim_gpyro4[-1] - zdim_4)/zdim_0
erreur5 = abs(zdim_gpyro5[-1] - zdim_5)/zdim_0

print(" Error in final thickness:")
print(f"# Case 1- deformation_pure_condesed : {erreur1:.4f} %")
print(f"# Case 2- no_deformation: {erreur2:.4f} %")
print(f"# Case 3- swelling: {erreur3:.4f} %")
print(f"# Case 4- shrinking: {erreur4:.4f} %")
print(f"# Case 5- non-charring: {erreur5:.4f} %")

erreur=np.max([erreur1,erreur2, erreur3, erreur4,erreur5])

#%%

# Validation threshold
threshold = 0.001
if erreur <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
