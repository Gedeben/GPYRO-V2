import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(sys.argv[0] if '__file__' in locals() else os.getcwd()))
output_filename = os.path.join(SCRIPT_DIR, "analytical_results.csv")


def find_eigenvalues(num_eigenvalues, delta, k, hc):
    """
    Find the first 'num_eigenvalues' eigenvalues λ_m from:
        cot(λ_m δ) = (k / hc) λ_m
    """
    eigenvalues = []
    equation = lambda lam: np.cos(lam * delta) / np.sin(lam * delta) - (k / hc) * lam

    for m in range(1, num_eigenvalues + 1):
        bracket_low = (m - 1) * np.pi / delta + 1e-6
        bracket_high = m * np.pi / delta - 1e-6

        try:
            sol = root_scalar(equation, bracket=[bracket_low, bracket_high])
            eigenvalues.append(sol.root)
        except ValueError:
            print(f"Warning: Could not find eigenvalue for m={m}")

    return np.array(eigenvalues)


def calculate_temperature(times, depths, T0, q_e, k, rho, cp, kappa, delta, hc, num_terms=100):
    """
    Implements Eq. (3.136a) exactly.
    """
    alpha = k / (rho * cp)
    depths_np = np.array(depths)

    # Find eigenvalues λ_m
    lambdas = find_eigenvalues(num_terms, delta, k, hc)

    # Set up result array
    T_result = np.zeros((len(times), len(depths)))

    # Loop over time and depth
    for it, t in enumerate(times):
        for iz, z in enumerate(depths_np):
            series_sum = 0.0
            for lm in lambdas:
                # 1 / [ λ_m (k^2 + λ_m^2) ]
                pref_series = 1.0 / (lm * (kappa**2 + lm**2))

                # numerator: cos(λ δ) + (λ/κ) sin(λ δ) - exp(-κ δ)
                num_first = np.cos(lm * delta) + (lm / kappa) * np.sin(lm * delta) - np.exp(-kappa * delta)

                # denominator: λ δ + 0.5 sin(2 δ λ)
                denom = (lm * delta) + (0.5 * np.sin(2 * delta * lm))

                # multiply by cos(λ (δ - z))
                spatial_part = np.cos(lm * (delta - z))

                # decay term
                time_part = (1.0 - np.exp(-alpha * lm**2 * t))

                # full term
                term = pref_series * (num_first / denom) * spatial_part * time_part
                series_sum += term

            # assemble temperature
            T_result[it, iz] = T0 + (2 * q_e * kappa**2 / k) * series_sum

    return T_result


if __name__ == '__main__':
    # Parameters
    T0 = 20.0       # Initial temp (°C)
    q_e = 5000.0    # Incident heat flux (W/m^2)
    k = 0.1         # Thermal conductivity (W/(m·K))
    rho = 100       # kg/m³
    cp = 1000       # J/(kg·K)
    hc = 10         # Convective coefficient (W/(m²·K))
    kappa = 24.0    # Radiation absorption coefficient (1/m)
    delta = 0.1     # Slab thickness (m)
    num_terms = 200

    depth_points = [0.0, 0.04, 0.1]
    time_points = np.linspace(0, 2000, 500)

    print("Calculating temperature evolution...")
    T_evolution = calculate_temperature(
        times=time_points,
        depths=depth_points,
        T0=T0,
        q_e=q_e,
        k=k,
        rho=rho,
        cp=cp,
        kappa=kappa,
        delta=delta,
        hc=hc,
        num_terms=num_terms
    )
    print("Calculation complete.")

    # Save to CSV
    column_headers = ['Time (s)'] + [f'Temp_at_{d:.4f}m' for d in depth_points]
    df = pd.DataFrame(T_evolution, columns=column_headers[1:])
    df.insert(0, 'Time (s)', time_points)
    df.to_csv(output_filename, index=False, float_format='%.4f')

    print(f"Results saved to '{output_filename}'")
