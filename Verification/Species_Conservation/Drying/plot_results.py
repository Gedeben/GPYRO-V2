import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import numpy as np
from matplotlib.lines import Line2D



import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

# File paths
file_path = os.path.join(SCRIPT_DIR, 'drying_summary_01_0001.csv')

#%%


style_traces = {
    "Gpyro": {"linewidth": 3},
    "Analytical": {'marker': '+',  'linestyle': 'none', 'markersize': 16, 'markevery': 1, 'markeredgewidth': 4},

}

case_styles = [
 {'color': 'blue'},
{'color': 'black'},
{'color': 'red'},
{'color': 'green'},
{'color': 'purple'},
{'color': 'orange'}
]

# Matplotlib style settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 25
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['legend.fontsize'] = 25
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['figure.figsize'] = (12,8)



def custom_cumtrapz(y, x, initial=0):
    """
    Calcul de l'intégrale cumulative via la règle des trapèzes.
    Équivalent à scipy.integrate.cumtrapz pour des données non uniformément espacées.
    
    Parameters:
    - y: array-like, valeurs à intégrer
    - x: array-like, abscisses (doivent être de même longueur que y)
    - initial: float, valeur initiale de l'intégrale (par défaut 0)
    
    Returns:
    - result: array, intégrale cumulative
    """
    if len(x) != len(y):
        raise ValueError("x et y doivent avoir la même longueur.")
    
    dx = np.diff(x)  # Différences entre les abscisses
    area = dx * (y[:-1] + y[1:]) / 2  # Aire des trapèzes entre chaque paire de points
    cumulative = np.cumsum(area)  # Cumul des aires
    result = np.concatenate([[initial], cumulative])  # Ajout de la valeur initiale
    return result


#%%
# Load data
try:
    data = pd.read_csv(file_path)
except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the path. {e}")
    sys.exit(1)

# Time vector
t = data['t'].values

T = data["002_TEMPERATURE( 0.0000_ 0.0000_ 0.0000)"].values


theoretical_water_content = 1-0.8
theoretical_emitted_gas = 0.8-0.48
theoretical_char_content = 1 - theoretical_water_content - theoretical_emitted_gas


gpyro_water_content = custom_cumtrapz(data["007_MLR( 0.0000_ 0.0000_ 0.0000)"].values/1000, t, initial=0)
gpyro_emitted_gas = custom_cumtrapz(data["008_MLR( 0.0000_ 0.0000_ 0.0000)"].values/1000, t, initial=0)
mass_wet = data["003_YI( 0.0000_ 0.0000_ 0.0000)"].values*(data["009_TOTAL_MASS( 0.0000_ 0.0000_ 0.0000)"].values/1000)
mass_dray = data["004_YI( 0.0000_ 0.0000_ 0.0000)"].values*(data["009_TOTAL_MASS( 0.0000_ 0.0000_ 0.0000)"].values/1000)
mass_char = data["005_YI( 0.0000_ 0.0000_ 0.0000)"].values*(data["009_TOTAL_MASS( 0.0000_ 0.0000_ 0.0000)"].values/1000)


theoretical_total = 1


Total_mass=mass_wet+mass_dray+mass_char+gpyro_emitted_gas+gpyro_water_content


n=5
tp =[t[int(np.size(mass_wet)*i/20)-1]/60 for i in range(16,21)]
fig = plt.figure()

# Axe principal (gauche)
ax = fig.add_axes([0.08, 0.12, 0.75, 0.8])  # [left, bottom, width, height]

# Axe fantôme pour les légendes (droite)
ax_leg = fig.add_axes([0.85, 0.12, 0.30, 0.8])
ax_leg.axis("off")

ax.plot(t/60, mass_wet,**style_traces["Gpyro"],**case_styles[0])
ax.plot(t/60, mass_dray, **style_traces["Gpyro"],**case_styles[1])
ax.plot(t/60, mass_char, **style_traces["Gpyro"],**case_styles[2])
ax.plot(t/60, gpyro_water_content,**style_traces["Gpyro"],**case_styles[3])
ax.plot(t/60, gpyro_emitted_gas,**style_traces["Gpyro"],**case_styles[4])
ax.plot(t/60, Total_mass,**style_traces["Gpyro"],**case_styles[5])

ax.plot(tp,0*np.ones(n), **style_traces["Analytical"],**case_styles[0])
ax.plot(tp,0*np.ones(n), **style_traces["Analytical"],**case_styles[1])
ax.plot(tp,theoretical_char_content*np.ones(n), **style_traces["Analytical"],**case_styles[2])

ax.plot(tp,theoretical_water_content*np.ones(n), **style_traces["Analytical"],**case_styles[3])
ax.plot(tp,theoretical_emitted_gas*np.ones(n), **style_traces["Analytical"],**case_styles[4])
ax.plot(tp,1*np.ones(n), **style_traces["Analytical"],**case_styles[5])

ax.set_xlabel("Time (min)")
ax.set_ylabel("Mass (kg)")
ax.grid(False)

# ================== LEGENDES ==================
legend_codes = [
    Line2D([0], [0], **style_traces["Gpyro"], label="Gpyro", color="gray"),
    Line2D([0], [0], **style_traces["Analytical"], label="Analytical \nfinal mass", color="gray"),
]

case_labels = [
    r"Solid${_{wet}}$",
    r"Solid${_{dry}}$",
    r"$Char$",
    "Emitted "+r"$H_2O$",
    "Emitted "+r"$pyrolysate$",
    "Total mass"
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
output_file = Path(SCRIPT_DIR) / "Mass_drying.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')


#%%
# Function to compute absolute error between simulation and analytical solution

erreur_water = np.abs(gpyro_water_content - theoretical_water_content)/theoretical_water_content * 100
erreur_gas = np.abs(gpyro_emitted_gas - theoretical_emitted_gas)/theoretical_emitted_gas * 100
erreur_char = np.abs(mass_char - theoretical_char_content)/theoretical_char_content * 100
erreur_total = np.abs(Total_mass - 1)/theoretical_total * 100

erreur = max(erreur_water[-1], erreur_gas[-1], erreur_char[-1], erreur_total[-1])
print(f"\nMaximum relative error at equilibrium: {erreur} %")

# Validation threshold
threshold = 5  # 5%
if erreur <= threshold:
    print("Validation PASSED.")
    sys.exit(0)
else:
    print("Validation FAILED.")
    sys.exit(1)
