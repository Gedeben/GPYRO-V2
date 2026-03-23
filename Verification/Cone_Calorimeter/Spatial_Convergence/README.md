# Spatial Convergence Verification Case

## Overview

This validation case is designed to verify the **spatial convergence** behavior of the Gpyro solvers. It reproduces a one-dimensional pyrolysis scenario under cone calorimeter conditions using increasingly refined spatial grids while keeping the time step fixed at \(\Delta t = 0.05\,\text{s}\).

The objective is to assess how the accuracy of the numerical solution improves with finer spatial resolution. The spatial discretizations range from \(1\,\text{mm}\) down to \(0.0001\,\text{mm}\).

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable [--full]
```

- If the `--full` option is provided, **all** simulations are executed, including the **two most computationally expensive cases**: `0.001mm_dt0.05s` and `0.0001mm_dt0.05s`.
- Without `--full`, only the moderate and fast simulations are run.

This command will:
- Launch the appropriate Gpyro simulations.
- Post-process the results to extract convergence metrics.
- Compile a LaTeX report into a PDF format summarizing the setup and convergence analysis.

### 2. Run a single case manually

```bash
cd 0.01mm_dt0.05s
/path/to/gpyro/executable Reference_CC.data
```

This runs the simulation for a single mesh resolution without post-processing or report generation.

## Folder Structure

| File/Folder               | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| `run.sh`                  | Bash script to run the simulations and optionally include heavy cases.      |
| `plot_results.py`         | Python script for computing errors and generating plots.                    |
| `spatial_convergence.tex`| LaTeX source file for the convergence report.                               |
| `spatial_convergence.pdf`| Precompiled version of the report.                                          |
| `*/Reference_CC.data`     | Input file for each spatial discretization test case.                       |
| `0.0001mm_dt0.05s/`       | Simulation case with finest spatial resolution (already run due to cost).   |
| `0.001mm_dt0.05s/`        | High-resolution simulation (already run due to cost).                       |
| `0.01mm_dt0.05s/`, ..., `1mm_dt0.05s/` | Other mesh refinement cases with moderate computational cost. |


