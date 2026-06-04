import numpy as np
from scipy.optimize import brentq
import math
import csv

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

results_file_path = os.path.join(SCRIPT_DIR, "analytical_results.csv")

def solve_transcendental_equation(Bi, num_roots=200):
    """
    Finds the positive roots of the equation beta * tan(beta) = Bi.
    """
    roots = []
    # Search for roots in intervals ((n-1)*pi, (n-1)*pi + pi/2)
    # For n=1, interval is (eps, pi/2 - eps)
    # Correct intervals for beta * tan(beta) = Bi (where Bi > 0):
    # beta must be in ( (m-1)*pi, (m-1)*pi + pi/2 ) for m=1, 2, ...
    
    # Initial small epsilon to avoid division by zero or exact boundaries
    eps = 1e-9 

    # First root in (eps, pi/2 - eps)
    try:
        root = brentq(lambda beta: beta * np.tan(beta) - Bi, eps, np.pi/2 - eps)
        roots.append(root)
    except ValueError:
        # This can happen if Bi is very small, root is close to 0
        # or if Bi is large and root is very close to pi/2
        if Bi < 1e-6 and Bi >= 0:
             roots.append(np.sqrt(Bi) if Bi > 0 else eps) # Approximate for very small Bi
        elif Bi > 1e5: # if Bi is very large, root is very close to pi/2
            try:
                # Try a tighter interval if default fails for large Bi
                root = brentq(lambda beta: beta * np.tan(beta) - Bi, np.pi/2 - 10*eps, np.pi/2 - eps)
                roots.append(root)
            except ValueError:
                print(f"Warning: Could not find the first root for Bi={Bi}. Approximating as pi/2.")
                roots.append(np.pi/2 - eps) # Fallback, might not be accurate
        else:
            print(f"Warning: Could not find the first root for Bi={Bi} in (eps, pi/2-eps).")


    # Subsequent roots
    for n in range(2, num_roots + 1):
        lower_bound = (n - 1) * np.pi + eps
        upper_bound = (n - 1) * np.pi + np.pi/2 - eps
        if lower_bound >= upper_bound: 
            continue
        try:
            root = brentq(lambda beta: beta * np.tan(beta) - Bi, lower_bound, upper_bound)
            roots.append(root)
        except ValueError:
            print(f"Warning: Could not find root {n} for Bi={Bi} in interval ({lower_bound:.4f}, {upper_bound:.4f}). Skipping.")
            pass 
            
    return np.array(roots)

def calculate_temperature(x, t, L, alpha, Tg, T_initial, Bi_val, beta_roots):
    """
    Calculates temperature T(x,t) using the series solution.
    """
    if t == 0:
        return T_initial

    theta_sum = 0
    dimensionless_pos = x / L
    fourier_number = alpha * t / (L**2)

    for beta_n in beta_roots:
        if beta_n == 0: 
             continue

        Cn_num = 4 * np.sin(beta_n)
        Cn_den = 2 * beta_n + np.sin(2 * beta_n)

        if abs(Cn_den) < 1e-9: 
            C_n = 0 
        else:
            C_n = Cn_num / Cn_den
        
        term = C_n * np.cos(beta_n * dimensionless_pos) * np.exp(-beta_n**2 * fourier_number)
        theta_sum += term
    
    temperature = Tg + (T_initial - Tg) * theta_sum
    return temperature

def main():
    # Case A Parameters
    L = 1  # m, thickness of the slab
    k = 1  # W/(mK), thermal conductivity
    rho = 1000  # kg/m^3, density
    c = 1  # J/(kgK), specific heat capacity (1 kJ/kgK = 1000 J/kgK)
    h = 1  # W/(m^2K), convective heat transfer coefficient
    
    T_initial = 1000  # °C, initial temperature of the slab
    Tg = 0  # °C, ambient gas temperature
    
    # Calculate thermal diffusivity (alpha) and Biot number (Bi)
    alpha = k / (rho * c)
    Bi = (h * L) / k
    
    #print(f"Calculated parameters:")
    #print(f"Alpha (thermal diffusivity): {alpha:.2e} m^2/s")
    #print(f"Biot number (Bi): {Bi:.2f}")
    #print("-" * 30)

    # Number of roots (terms in series)
    num_terms = 200 
    beta_roots = solve_transcendental_equation(Bi, num_roots=num_terms)
    
    if len(beta_roots) < 5: 
        print(f"Warning: Only {len(beta_roots)} roots found. Results might be inaccurate.")
        if len(beta_roots) == 0 and Bi > 1e-7: 
             print("CRITICAL: No roots found for Bi > 0. Check root finding logic or Bi value.")
             return

    # Time points for evaluation (seconds)
    time_points = np.linspace(0,1800,30)
                   
    x_coords_m = {
        "Back": 0.0,
        "50 cm": 0.5,
        "Front": L 
    }
    output_column_names = ["Back", "50 cm", "Front"]
    output_x_values = [x_coords_m[name] for name in output_column_names]

    # --- CSV Output Setup ---
    csv_filename = results_file_path
    header_row_csv = ["Time"] + output_column_names
    
    with open(csv_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(header_row_csv)

        # Calculate and print/write temperatures
        for t_sec in time_points:
            results_row_data = [f"{t_sec}"] # Start with time for this row
            for x_m in output_x_values:
                temp = calculate_temperature(x_m, t_sec, L, alpha, Tg, T_initial, Bi, beta_roots)
                results_row_data.append(f"{temp:.2f}")
            csv_writer.writerow(results_row_data) # Write to CSV

    #print(f"Results saved to {csv_filename}")

if __name__ == "__main__":
    main()
