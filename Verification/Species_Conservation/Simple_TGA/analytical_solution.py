# main_script.py

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from pathlib import Path
from typing import Dict, Any

# ---------------------------
# CONSTANTS
# ---------------------------
# These parameters define the material's properties and the experimental conditions.
PYROLYSIS_PARAMS: Dict[str, Any] = {
    # Kinetic parameters (Arrhenius equation)
    "A": 8.7E+12,      # Pre-exponential factor [1/s]
    "E": 163e3,    # Activation energy [J/mol]
    "R": 8.314,    # Universal gas constant [J/mol/K]
    
    # Initial and boundary conditions
    "m0": 1.0,     # Initial mass [kg]
    "T0": 300,     # Initial temperature [K]
    "Tmax": 1000,  # Maximum temperature [K]
    "beta": 5.0/60,   # Constant heating rate [K/s]
    
    # Numerical parameters
    "num_points": 1000 # Number of calculation points
}

# ---------------------------
# CORE CALCULATION
# ---------------------------

def calculate_pyrolysis(params: Dict[str, Any]) -> pd.DataFrame:
    """
    Calculates the analytical solution for single-step, first-order pyrolysis
    under a constant heating rate.

    This function models how a material's mass changes as it's heated over time.

    Args:
        params (Dict[str, Any]): A dictionary containing all simulation parameters.

    Returns:
        pd.DataFrame: A DataFrame containing the time-resolved results of the simulation.
    """
    # --- Unpack parameters for easier access ---
    A, E, R = params["A"], params["E"], params["R"]
    m0, T0, Tmax = params["m0"], params["T0"], params["Tmax"]
    beta, num_points = params["beta"], params["num_points"]

    # --- Set up independent variables: Temperature and Time ---
    # Create a linear array of temperatures from start to max
    T = np.linspace(T0, Tmax, num_points)
    # Calculate the corresponding time array based on the constant heating rate
    t = (T - T0) / beta

    # --- Calculate key physical quantities ---
    # 1. Rate Constant (k): Governed by the Arrhenius equation.
    # It describes how fast the reaction occurs at a given temperature.
    k = A * np.exp(-E / (R * T))

    # 2. Integral of k/beta: This term appears in the analytical solution for mass.
    # We integrate numerically using the trapezoidal rule.
    # The 'initial=0' adds the starting point of the integration.
    integral = cumulative_trapezoid(k / beta, T, initial=0)

    # 3. Mass (m): The remaining mass as a function of temperature.
    m = m0 * np.exp(-integral)

    # 4. Mass Fraction (Y): The normalized mass (m/m0).
    Y = m / m0

    # 5. Mass Loss Rate (MLR): The positive rate at which mass is lost.
    # MLR = -dm/dt. For a first-order reaction, MLR = k * m.
    MLR = k * m * 1000

    # --- Assemble results into a DataFrame ---
    results_df = pd.DataFrame({
        "Time (s)": t,
        "Temperature (K)": T,
        "Rate Constant k (1/s)": k,
        "Mass Fraction Y (-)": Y,
        "Mass m (kg)": m,
        "MLR = -dm/dt (g/s)": MLR
    })
    
    return results_df

# ---------------------------
# DATA OUTPUT
# ---------------------------

def save_data_to_csv(df: pd.DataFrame, filename: str) -> None:
    """
    Saves the provided DataFrame to a CSV file in the same directory as the script.

    Args:
        df (pd.DataFrame): The DataFrame to be saved.
        filename (str): The desired name for the output file.
    """
    # Create a file path relative to the current script location
    output_path = Path(__file__).resolve().parent / filename
    
    # Save the DataFrame to a CSV file, without the index column
    df.to_csv(output_path, index=False, float_format='%.6g')
    
    print(f"Results successfully saved to: {output_path}")

# ---------------------------
# MAIN EXECUTION BLOCK
# ---------------------------

if __name__ == "__main__":
    """
    This block runs when the script is executed directly.
    """
    # Perform the pyrolysis calculation
    pyrolysis_results = calculate_pyrolysis(PYROLYSIS_PARAMS)
    
    # Save the results to a file
    save_data_to_csv(pyrolysis_results, "analytical_results.csv")