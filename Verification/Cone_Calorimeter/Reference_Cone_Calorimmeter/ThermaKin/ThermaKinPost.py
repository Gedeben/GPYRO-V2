# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re





filename = 'reference_cc.txt'  # <-- Ton fichier ici






# --- Fonction de parsing ---
def parse_thermakin_output(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    time_blocks = {}
    current_time = None
    reading_profiles = False
    read_species_list = True
    profiles = []
    boundary_info = {'TOP': {}, 'BOTTOM': {}}
    total_thickness = None
    total_mass = None

    for line in lines:
        line = line.strip()

        if line.startswith("Time [s]"):
            if current_time is not None:
                time_blocks[current_time] = {
                    "profiles": np.array(profiles),
                    "boundary": boundary_info,
                    "thickness": total_thickness,
                    "mass": total_mass
                }
            current_time = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", line)[0])
            profiles = []
            boundary_info = {'TOP': {}, 'BOTTOM': {}}
            total_thickness = None
            total_mass = None
            reading_profiles = False

        elif line.startswith("FROM TOP"):
            reading_profiles = True

        elif line.startswith("BOUNDARY"):
            reading_profiles = False
            
            if read_species_list :
                # Extraire les espèces depuis cette ligne
                header_parts = line.strip().split(':')
                if len(header_parts) > 1:
                    species_list = header_parts[1].strip().split()
                    read_species_list= False
                print (species_list)
    

        elif line.startswith("Total thickness"):
            reading_profiles = False
            total_thickness = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", line)[0])

        elif line.startswith("Total mass"):
            reading_profiles = False
            total_mass = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", line)[0])

        elif line.startswith("***"):
            reading_profiles = False

        elif reading_profiles:
            if line:
                try:
                    values = list(map(float, line.split()))
                    profiles.append(values)
                except ValueError:
                    reading_profiles = False

        elif line.startswith("TOP") or line.startswith("BOTTOM"):
            parts = line.split()
            side = parts[0]
            area = float(parts[1])
            heat_flow_in = float(parts[2])
            mass_flows = list(map(float, parts[3:]))
            boundary_info[side] = {
                'area': area,
                'heat_flow_in': heat_flow_in,
                'mass_flow_out': dict(zip(species_list, mass_flows))
            }

    if current_time is not None and current_time not in time_blocks:
        time_blocks[current_time] = {
            "profiles": np.array(profiles),
            "boundary": boundary_info,
            "thickness": total_thickness,
            "mass": total_mass
        }

    return time_blocks

# --- Fonction pour construire les (x, y) à partir des profils ---
def extract_profile_xy(time_blocks, time_to_plot, quantity='Temperature', species=None):
    if time_to_plot not in time_blocks:
        print(f"Time {time_to_plot} not found.")
        return None, None

    block = time_blocks[time_to_plot]
    profiles = block["profiles"]

    if profiles is None or len(profiles) == 0:
        print(f"No profile data at {time_to_plot}s.")
        return None, None

    z = profiles[:, 0]  # Depth from top

    if quantity.lower() == 'temperature':
        y = profiles[:, 1]
    elif quantity.lower() == 'concentration' and species is not None:
        species_list = list(time_blocks[0]["boundary"]["TOP"]["mass_flow_out"].keys())
        if species not in species_list:
            print(f"Species '{species}' not found. Available: {species_list}")
            return None, None
        idx = species_list.index(species) + 2  # Columns: z, T, then concentrations
        y = profiles[:, idx]
    else:
        print(f"Unsupported quantity '{quantity}'. Available : Temperature, concentration(species)")
        return None, None

    return y, z

# --- Fonction pour construire les (x, y) scalaires ---
def extract_scalar_xy(time_blocks, scalar='thickness', boundary='TOP', species=None):
    times = sorted(time_blocks.keys())
    values = []
    if scalar == 'mass_flow_out':
        species_list = list(time_blocks[0]["boundary"]["TOP"]["mass_flow_out"].keys())
        if (species not in species_list):
            scalar
            print(f"Species '{species}' not found. Available: {species_list}")
            return None, None
    for t in times:
        block = time_blocks[t]
        if scalar == 'thickness':
            values.append(block["thickness"])
        elif scalar == 'mass':
            values.append(block["mass"])
        elif scalar == 'heat_flow_in':
            values.append(block["boundary"][boundary]['heat_flow_in'])
        elif scalar == 'mass_flow_out' and species is not None:
            values.append(block["boundary"][boundary]['mass_flow_out'].get(species, 0.0))
        else:
            raise ValueError(f"Unknown scalar {scalar}, available : thickness, mass, heat_flow_in, mass_flow_out(species)")

    return times, values




def extract_point_quantity_over_time(time_blocks, z_point, quantity='Temperature', species=None):
    times = sorted(time_blocks.keys())
    values = []

    # On suppose que toutes les grilles z sont identiques (à vérifier si besoin)
    for t in times:
        block = time_blocks[t]
        profiles = block["profiles"]

        if profiles is None or len(profiles) == 0:
            values.append(np.nan)
            continue

        z = profiles[:, 0]

        # Sélection de la bonne colonne
        if quantity.lower() == 'temperature':
            y = profiles[:, 1]
        elif quantity.lower() == 'concentration' and species is not None:
            species_list = list(block["boundary"]["TOP"]["mass_flow_out"].keys())
            if species not in species_list:
                print(f"[t={t}s] Species '{species}' not found. Available: {species_list}")
                values.append(np.nan)
                continue
            idx = species_list.index(species) + 2
            y = profiles[:, idx]
        else:
            raise ValueError(f"Unsupported quantity '{quantity}'. Available: Temperature, concentration(species)")

        # Interpolation linéaire en z_point
        if z_point < np.min(z) or z_point > np.max(z):
            values.append(np.nan)  # z hors domaine
        else:
            y_interp = np.interp(z_point, z, y)
            values.append(y_interp)

    return times, values




# --- Fonction de plot scientifique ---
def scientific_plot(x, y, xlabel='X', ylabel='Y', title='', xlim=None, ylim=None, color='k', label=None):
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['legend.fontsize'] = 18
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color=color, linewidth=2, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(False)
    if label:
        plt.legend(frameon=False)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.tight_layout()
    plt.show()

# =====================
# === Exemple d'utilisation ===
# =====================

species_to_plot='GAS_1'
input_time_to_plot= 300 #s
time_blocks = parse_thermakin_output(filename)

#%%

# Exemple : tracer l'évolution de la température à z=0.002 m (2 mm)
z_target = 2.500000e-04  # [m]
times, tz0 = extract_point_quantity_over_time(
    time_blocks,
    z_point=z_target,
    quantity='Temperature'
)
z_target = 0.003  # [m]
times, tz3 = extract_point_quantity_over_time(
    time_blocks,
    z_point=z_target,
    quantity='Temperature'
)

z_target = 0.006  # [m]
times, tz6 = extract_point_quantity_over_time(
    time_blocks,
    z_point=z_target,
    quantity='Temperature'
)

plt.plot(times,tz0)
plt.plot(times,tz3)
plt.plot(times,tz6)
#%%
scientific_plot(times, tz6,
                xlabel='Time (s)', ylabel='Temperature (K)',
                title=f'Temperature at z={z_target*1000:.1f} mm')

#%%

# Exemple : tracer l'évolution de la concentration de GAS_1 à z = 0.005 m
times, c_A = extract_point_quantity_over_time(
    time_blocks,
    z_point=0.003,
    quantity='concentration',
    species='MAT_A'
)

times, c_B = extract_point_quantity_over_time(
    time_blocks,
    z_point=0.003,
    quantity='concentration',
    species='MAT_B'
)

times, c_C = extract_point_quantity_over_time(
    time_blocks,
    z_point=0.003,
    quantity='concentration',
    species='MAT_C'
)

plt.plot(times,c_A)
plt.plot(times,c_B)
plt.plot(times,c_C)

#%%
# Exemples :
# Plot profil de température
y, z = extract_profile_xy(time_blocks, time_to_plot=input_time_to_plot, quantity='Temperature')
scientific_plot(z, y, xlabel='Depth from top (m)', ylabel='Temperature (K)', title=f'Temperature Profile at t={input_time_to_plot} s')
#%%
# Plot profil de concentration de PMMA
y, z = extract_profile_xy(time_blocks, time_to_plot=input_time_to_plot, quantity='Concentration', species=species_to_plot)
scientific_plot(z, y, xlabel='Depth from top (m)', ylabel='Concentration (kg/m³)', title=f'{species_to_plot} Concentration Profile at t={input_time_to_plot} s')
#%%
# Plot scalaires : masse vs temps
times, thikness = extract_scalar_xy(time_blocks, scalar='thickness')
scientific_plot(times, thikness, xlabel='Time (s)', ylabel='thikness (m)', title='thikness Evolution')

"""
m_pmma=(np.array(thikness)- thikness[-1])*120
scientific_plot(times, m_pmma, xlabel='Time (s)', ylabel='mass (kg)', title='mass Evolution')
mlr= -(m_pmma[1:]-m_pmma[:-1])
scientific_plot(times[:-1], mlr, xlabel='Time (s)', ylabel='mlr', title='mass Evolution')
"""

#%%
# Plot scalaires : flux de masse PMMA_g
times, MLR1 = extract_scalar_xy(time_blocks, scalar='mass_flow_out', boundary='TOP', species=species_to_plot)
scientific_plot(times, MLR1, ylim=None, xlabel='Time (s)', ylabel='MLR  (kg/s)', title=f'Mass flow out {species_to_plot}')

times, MLR2 = extract_scalar_xy(time_blocks, scalar='mass_flow_out', boundary='TOP', species='GAS_2')
scientific_plot(times, MLR2, ylim=None, xlabel='Time (s)', ylabel='MLR  (kg/s)', title=f'Mass flow out {species_to_plot}')
MLR=np.array(MLR1)+np.array(MLR2)

scientific_plot(times, MLR,ylim=None, xlabel='Time (s)', ylabel='MLR  (kg/s)', title=f'Mass flow out {species_to_plot}')

#%%
Temp=[]
for time in times:
    y, T = extract_profile_xy(time_blocks, time_to_plot=time, quantity='Temperature')
    Temp.append(T[0])
    
plt.plot(times, Temp)




#%%
import csv

# Génération du fichier CSV
csv_filename = 'output_ThermaKin.csv'
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['time (s)', 'MLR (kg/s)', 'thikness (m)','T(z=0)','T(z=3)','T(z=6)','C_A (kg/m^3)','C_B (kg/m^3)', 'C_C (kg/m^3)'])  # En-têtes

    # Supposer que times, MLR et thikness ont la même longueur
    for t, m, h, z0, z3, z6, A, B, C in zip(times, MLR, thikness, tz0, tz3, tz6, c_A, c_B, c_C):
        writer.writerow([t, m, h, z0, z3, z6, A, B, C])

