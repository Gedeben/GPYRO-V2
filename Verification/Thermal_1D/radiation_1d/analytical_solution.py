import numpy as np
from scipy.special import erfc
import pandas as pd
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

results_file_path = os.path.join(SCRIPT_DIR, "analytical_results.csv")

def calculate_temperature_heat_flux(q0, k, rho, cp, Ti, x, t):
    """
    Calculates the temperature in a semi-infinite solid with constant surface heat flux.

    Args:
        q0 (float): Constant heat flux at the surface (W/m^2).
        k (float): Thermal conductivity of the material (W/(m*K)).
        rho (float): Density of the material [kg/m³].
        cp (float): Specific heat capacity of the material [J/kg·K].
        Ti (float): Initial uniform temperature of the solid (°C or K).
        x (float): Position (depth) from the surface (m).
        t (float): Time elapsed (s).

    Returns:
        float: The temperature T(x, t) at the specified position and time.
    """
    # The equation is not defined for t=0, so we handle this case.
    if t <= 0:
        return Ti

    # Calculate the two terms of the temperature rise equation
    alpha = k / (rho * cp)
    
    term1 = (2 * q0 / k) * np.sqrt(alpha * t / np.pi) * np.exp(-x**2 / (4 * alpha * t))
    term2 = (q0 * x / k) * erfc(x / (2 * np.sqrt(alpha * t)))

    # The temperature rise is delta_T = T(x, t) - Ti
    delta_T = term1 - term2

    # The final temperature is T = Ti + delta_T
    T_xt = Ti + delta_T
    
    return T_xt

if __name__ == "__main__":
    thermal_conductivity = 0.186   # k in W/(m*K)
    density = 380		   # rho in kg/m³
    heat_capacity = 1764	   # cp in J/kg.K
    initial_temp = 27              # Ti in °C

    # Heating conditions
    heat_flux = 2000              # q0 in W/m^2 (50 kW/m^2)
    time_elapsed = 300            # t in seconds (5 minutes)
    num_steps = 150		  # number of time steps
    time_points = np.linspace(0, time_elapsed, num_steps)
       
    # Positions where we want to find the temperature
    depths = [0.0, 0.008, 0.016, 0.04]                   # x in meters

    # Calculate the temperature
    results = []
    for t in time_points:
        current_row = {'Time (s)': t}
        for depth in depths:
            temperature = calculate_temperature_heat_flux(
                q0=heat_flux,
                k=thermal_conductivity,
                rho=density,
                cp=heat_capacity,
                Ti=initial_temp,
                x=depth,
                t=t
            )
            column_name = f"Depth {depth * 100:.1f} cm"
            current_row[column_name] = temperature
        results.append(current_row)
    results_df = pd.DataFrame(results)
    results_df.to_csv(results_file_path, index=False)

    print(f"Results saved to {results_file_path}")