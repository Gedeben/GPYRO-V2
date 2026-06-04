import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from pathlib import Path
from typing import Dict, Any

# ---------------------------
# CONSTANTS
# ---------------------------
OXIDATION_PARAMS: Dict[str, Any] = {
    # Kinetic parameters (Arrhenius equation)
    "A": 7E6,      # Pre-exponential factor [1/s]
    "E": 190e3,        # Activation energy [J/mol]
    "R": 8.314,        # Universal gas constant [J/mol/K]

    # Initial and boundary conditions
    "rho_a": 1000,     # Initial sample density [kg/m³]
    "rho_b": 800,      # Final residue density [kg/m³]
    "T0": 300,         # Initial temperature [K]
    "tmax": 18000,     # Maximum time [s]
    "beta": 5.0/60,    # Constant heating rate [K/s]
    "YO2": [0, 0.043, 0.082, 0.205],   # Ambient oxygen mass fraction [-]
    # Numerical parameters
    "num_points": 1000 # Number of calculation points
}

# ---------------------------
# CORE CALCULATION
# ---------------------------

def calculate_oxidation(params: Dict[str, Any]) -> pd.DataFrame:
    """
    Calculates the analytical solution for single-step, first-order oxidation
    under a constant heating rate, accounting for a non-reacting residue.

    Args:
        params (Dict[str, Any]): A dictionary containing all simulation parameters.

    Returns:
        pd.DataFrame: A DataFrame containing the time-resolved results.
    """
    # --- Unpack parameters for easier access ---
    A, E, R = params["A"], params["E"], params["R"]
    rho_a, rho_b, T0, tmax, YO2 = params["rho_a"], params["rho_b"], params["T0"], params["tmax"], params["YO2"]
    beta, num_points = params["beta"], params["num_points"]

    # --- Correctly define reacting and residual mass fractions ---
    # The non-reacting portion that remains as char/residue (e.g., 0.8 or 80%)
    SF_char = rho_b / rho_a
    # The reacting portion of the initial mass (e.g., 0.2 or 20%)
    SF_gas = 1.0 - SF_char

    # --- Set up independent variables: Temperature and Time ---
    Tmax = T0 + beta * tmax # Final temperature based on heating rate and time
    T = np.linspace(T0, Tmax, num_points)
    t = (T - T0) / beta      # Time corresponding to each temperature point

    # --- Calculate key physical quantities ---
    # 1. Reaction rate constant k(T), including oxygen dependency.
    # The model assumes the rate is first-order with respect to O2.
    # For YO2=0, k=0, so no reaction occurs (pyrolysis case).
    k = [y * A * np.exp(-E / (R * T)) for y in YO2]

    # 2. Integral term for the analytical solution.
    # We integrate k/beta with respect to temperature T.
    integral_k_beta = [cumulative_trapezoid(k_i / beta, T, initial=0) for k_i in k]

    # 3. Decay of the Reacting Mass: The decay of the volatile fraction.
    m_reactant_decay = [np.exp(-integral_i) for integral_i in integral_k_beta]

    # 4. Total Mass (Corrected): The sum of the remaining reacting mass and the CONSTANT residue.
    # m_total = (SF_gas * m_reactant_decay) + SF_char
    m_total = [SF_gas * decay + SF_char for decay in m_reactant_decay]

    # 5. Mass Loss Rate (MLR) (Corrected): The positive rate at which mass is lost.
    # MLR = -dm/dt = k * m_reacting. Only the reacting mass contributes to MLR.
    # The current amount of reacting mass is (SF_gas * m_reactant_decay)
    m_reacting_current = [SF_gas * decay for decay in m_reactant_decay]
    MLR = [k_i * m_reacting_i for k_i, m_reacting_i in zip(k, m_reacting_current)]

    # --- Assemble results into a properly structured DataFrame (Corrected) ---
    results = {"Time (s)": t, "Temperature (K)": T}
    yo2_labels = [f"{y*100:.1f}%" for y in YO2] # e.g., "0.0%", "4.3%"

    for i, label in enumerate(yo2_labels):
        results[f'Total Mass {label} O2 (norm)'] = m_total[i]
        results[f'Reacting Mass Decay {label} O2 (norm)'] = m_reactant_decay[i]
        results[f'MLR {label} O2 (g/s)'] = MLR[i]*1000  # Convert from kg/s to g/s assuming initial mass is 1g

    results_df = pd.DataFrame(results)
    return results_df

# ---------------------------
# DATA OUTPUT
# ---------------------------

def save_data_to_csv(df: pd.DataFrame, filename: str) -> None:
    """
    Saves the provided DataFrame to a CSV file in the same directory as the script.
    """
    # Create the output directory if it doesn't exist
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    
    df.to_csv(output_path, index=False, float_format='%.6g')
    print(f"Results successfully saved to: {output_path}")

# ---------------------------
# MAIN EXECUTION BLOCK
# ---------------------------

if __name__ == "__main__":
    # Perform the oxidation calculation
    oxidation_results = calculate_oxidation(OXIDATION_PARAMS)
    
    # Save the results to a file
    save_data_to_csv(oxidation_results, "analytical_results.csv")