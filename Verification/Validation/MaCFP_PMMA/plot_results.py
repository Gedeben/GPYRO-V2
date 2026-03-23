#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 16:33:36 2025

@author: d90183
"""


#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from scipy.signal import savgol_filter

# Configuration des paramètres de style pour correspondre aux normes d’un article scientifique
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 30
mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['legend.fontsize'] = 30
mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30

# Chargement des données
Gpyro = pd.read_csv("PMMA_summary_01_0001.csv")

Exp = pd.read_csv("MaCFP-PMMA_Gasification_q50_Mass_R3.csv", skiprows=1)

S=1*3.1415*(0.035**2)
Exp['MLR'] = -np.gradient(Exp['[g]'].values, Exp['[s]'].values)/S


Exp['MLR'] = savgol_filter(Exp['MLR'].values, window_length=11, polyorder=2)


# Définition des styles de tracé pour chaque courbe
style_traces = {
    "GpyroV08": {"color": "blue", "linestyle": "--", "linewidth": 4},
    "Gpyro": {"color": "red", "linestyle": "-.", "linewidth": 4},
    "Exp": {"color": "black", "linestyle": " ", "marker":'+', "markersize": 15, 'markeredgewidth':2}
}

# Création de la figure
plt.figure(figsize=(12, 7.5))


plt.plot(Gpyro['t'].values, Gpyro['007_MLR( 0.0000_ 0.0000_ 0.0000)'].values,
         label="Gpyro", **style_traces["Gpyro"])



plt.plot(Exp['[s]'].values, 1*Exp['MLR'].values,
         label="Exp", **style_traces["Exp"], markevery=2)

# Mise en forme finale
plt.title("")
plt.ylabel("MLR "+r"$\mathrm{(g{\cdot}m^{-2}{\cdot}s^{-1})}$")

plt.xlabel("Time (s)")
plt.xlim([0, 450])
plt.grid(False)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("result_pmma_macfp.png")



#%%

# Interpolation de la courbe Exp sur les temps de Gpyro
exp_interp = np.interp(Gpyro['t'].values, Exp['[s]'].values, Exp['MLR'].values)



def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


rmse_value = rmse(Gpyro['007_MLR( 0.0000_ 0.0000_ 0.0000)'].values, exp_interp)



# Calcul du pourcentage de similarité (plus RMSE est petit, plus les courbes sont proches)
max_val = max(Gpyro['007_MLR( 0.0000_ 0.0000_ 0.0000)'].max(), exp_interp.max())
error_value = min(100, 100 * ( rmse_value / max_val))


print("error_value =",round(error_value,2),"%")