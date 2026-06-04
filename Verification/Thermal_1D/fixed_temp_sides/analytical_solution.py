import numpy as np
import pandas as pd
import os
import sys

# Define the directory where the script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

def calculate_temperature(x, t, L, k, rho, cp, T1, T2, T_i, n_terms=150):
    """
    Calculates transient temperature in a slab using thermal conductivity.

    Args:
        x (np.ndarray): Position array along the slab's thickness [m].
        t (float): Time at which to calculate the temperature [s].
        L (float): Thickness of the slab [m].
        k (float): Thermal conductivity of the material [W/m·K].
        rho (float): Density of the material [kg/m³].
        cp (float): Specific heat capacity of the material [J/kg·K].
        T1 (float): Fixed surface temperature at x=0 [K or °C].
        T2 (float): Fixed surface temperature at x=L [K or °C].
        T_i (float): Uniform initial temperature of the slab [K or °C].
        n_terms (int): Number of terms to use in the series approximation.

    Returns:
        np.ndarray: Temperature profile T(x) at the given time t.
    """
    # 1. Calculate thermal diffusivity (alpha) from the given properties
    alpha = k / (rho * cp)
    
    # 2. Calculate the steady-state linear temperature profile
    # The linear profile is calculated for each position in the x array.
    T_steady_state = T1 + (T2 - T1) * x / L
    
    # 3. Calculate the transient part of the solution
    transient_sum = np.zeros_like(x, dtype=float)
    
    # Sum the first n_terms of the series
    for n in range(1, n_terms + 1):
        Cn = (2 / (n * np.pi)) * ((T_i - T1) + (T_i - T2) * (-1)**(n+1))
        term = (
            Cn * np.sin(n * np.pi * x / L) * np.exp(-alpha * (n * np.pi / L)**2 * t)
        )
        transient_sum += term
        
    # 4. The final solution is the sum of the steady-state and transient parts
    return T_steady_state + transient_sum

# --- Main script to run the simulation and save to CSV ---
if __name__ == '__main__':
    # 1. Define physical and simulation parameters
    L = 0.1             # Slab thickness [m]
    T_i = 0          # Initial temperature [°C]
    T2 = 400         # Surface temperature at x=0 [°C]
    T1 = 20           # Surface temperature at x=L [°C]
    
    # Material properties (e.g., for Carbon Steel)
    k = 0.1             # Thermal conductivity [W/m·K]
    rho = 100.0         # Density [kg/m³]
    cp = 1000.0          # Specific heat capacity [J/kg·K]
    
    num_terms = 150    # Number of terms for the series

    # 2. Define spatial and temporal points for the calculation
    x_points = np.array([0, 0.04, 0.1])
    t_points = np.linspace(0, 14000, 100) # in seconds

    # 3. Prepare data for the CSV file in a new format (time-based rows)
    data_rows = []
    
    # Calculate the temperature profile for each time point and store it
    for t in t_points:
        T_profile = calculate_temperature(x_points, t, L, k, rho, cp, T1, T2, T_i, num_terms)
        
        # Create a dictionary for the current time step
        row_dict = {'Time (s)': t}
        for i, x in enumerate(x_points):
            # Create a column for each position
            row_dict[f'Temperature at x={x}m (°C)'] = T_profile[i]
        
        data_rows.append(row_dict)

    # 4. Create a pandas DataFrame from the list of dictionaries and save to CSV
    df = pd.DataFrame(data_rows)
    output_filename = os.path.join(SCRIPT_DIR, 'analytical_results.csv')
    df.to_csv(output_filename, index=False)

    print(f"Solution successfully saved to '{output_filename}'")
