import numpy as np
import pandas as pd

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
results_file_path = os.path.join(SCRIPT_DIR, 'analytical_results.csv')

def calculate_temperature_profile_plate(L_steel, L_ins, k_steel, k_ins, h, T_hot_air_C, T_cold_air_C):
    """
    Calculates the steady-state temperature profile through an insulated steel plate.

    Args:
        L_steel (float): Thickness of the steel plate in meters.
        L_ins (float): Thickness of the insulation layer in meters.
        k_steel (float): Thermal conductivity of steel in W/(m K).
        k_ins (float): Thermal conductivity of insulation in W/(m K).
        h (float): Heat transfer coefficient in W/(m^2 K).
        T_hot_air_C (float): Temperature of the hot air in Celsius.
        T_cold_air_C (float): Temperature of the cold air in Celsius.

    Returns:
        pandas.DataFrame: A DataFrame with 'depth' (in meters) and 'temp' (temperature in Celsius).
    """

    # Convert temperatures to Kelvin for calculations
    T_hot_air_K = T_hot_air_C + 273.15
    T_cold_air_K = T_cold_air_C + 273.15

    # Calculate total thermal resistance
    R_total = (1 / h) + (L_ins / k_ins) + (L_steel / k_steel) + (L_ins / k_ins) + (1 / h)

    # Calculate steady-state heat flux (q'')
    q_double_prime = (T_hot_air_K - T_cold_air_K) / R_total

    #print(f"Calculated Heat Flux (q''): {q_double_prime:.2f} W/m^2")

    # Calculate interface temperatures
    # Temperature at the surface of the first insulation layer (hot side)
    T0_K = T_hot_air_K - q_double_prime / h

    # Temperature at the interface between first insulation and steel
    T1_K = T0_K - q_double_prime * (L_ins / k_ins)

    # Temperature at the interface between steel and second insulation
    T2_K = T1_K - q_double_prime * (L_steel / k_steel)

    # Temperature at the surface of the second insulation layer (cold side)
    T3_K = T2_K - q_double_prime * (L_ins / k_ins)

    # Generate x-coordinates and temperatures for the profile
    num_points_per_layer = 5 # Number of points to sample within each solid layer
    x_coords = []
    temperatures_K = []

    # # Hot air and hot air-insulation interface
    # x_coords.append(-0.001) # Small offset for hot air point (for visualization)
    # temperatures_K.append(T_hot_air_K)
    # x_coords.append(0)
    # temperatures_K.append(T0_K)

    # First insulation layer
    x_ins1 = np.linspace(0, L_ins, num_points_per_layer, endpoint=False) # Exclude endpoint to avoid duplicates at interfaces
    T_ins1 = T0_K - q_double_prime * (x_ins1 / k_ins)
    x_coords.extend(x_ins1.tolist())
    temperatures_K.extend(T_ins1.tolist())
    # Add the interface point explicitly
    x_coords.append(L_ins)
    temperatures_K.append(T1_K)


    # Steel layer
    x_steel = np.linspace(L_ins, L_ins + L_steel, num_points_per_layer, endpoint=False)
    T_steel = T1_K - q_double_prime * ((x_steel - L_ins) / k_steel)
    x_coords.extend(x_steel.tolist())
    temperatures_K.extend(T_steel.tolist())
    # Add the interface point explicitly
    x_coords.append(L_ins + L_steel)
    temperatures_K.append(T2_K)

    # Second insulation layer
    x_ins2 = np.linspace(L_ins + L_steel, L_ins + L_steel + L_ins, num_points_per_layer, endpoint=False)
    T_ins2 = T2_K - q_double_prime * ((x_ins2 - (L_ins + L_steel)) / k_ins)
    x_coords.extend(x_ins2.tolist())
    temperatures_K.extend(T_ins2.tolist())
    # Add the interface point explicitly
    x_coords.append(L_ins + L_steel + L_ins)
    temperatures_K.append(T3_K)


    # # Cold air-insulation interface and cold air
    # x_coords.append(x_total_length + 0.001) # Small offset for cold air point (for visualization)
    # temperatures_K.append(T_cold_air_K)

    # Convert temperatures back to Celsius
    temperatures_C = [temp - 273.15 for temp in temperatures_K]

    # Create a Pandas DataFrame
    df = pd.DataFrame({
        'depth': x_coords,
        'temp': temperatures_C
    })

    # Sort by depth to ensure the profile is in order
    df = df.sort_values(by='depth').reset_index(drop=True)

    return df

# Given parameters
L_steel = 0.01  # 1 cm = 0.01 m
L_ins = 0.02    # 2 cm = 0.02 m
k_steel = 50    # W/(m K)
k_ins = 0.2     # W/(m K)
h = 10          # W/(m^2 K)
T_hot_air_C = 480 # °C
T_cold_air_C = 20 # °C

# Calculate the temperature profile and get it as a DataFrame
temperature_df = calculate_temperature_profile_plate(L_steel, L_ins, k_steel, k_ins, h, T_hot_air_C, T_cold_air_C)

# --- Save to CSV ---
output_filename = results_file_path
temperature_df.to_csv(output_filename, index = False)
