# Time Step Convergence Verification Case

## Overview

This verification case evaluates the **temporal convergence** of the Gpyro solver. A one-dimensional pyrolysis scenario is simulated using a fixed spatial resolution (\( \Delta x = 0.1\,\text{mm} \)) while varying the time step across several orders of magnitude — from \( 0.0005\,\text{s} \) to \( 10\,\text{s} \).

The objective is to assess how the numerical error evolves as the time step is refined and to estimate the convergence order of the time integration scheme.

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable [--full]
```

- If the `--full` option is provided, **all** simulations are executed, including the **he most computationally expensive cases**: `0.1mm_dt0.0005s`.
- Without `--full`, only the moderate and fast simulations are run.

This command will:
- Launch the appropriate Gpyro simulations.
- Post-process the results to extract convergence metrics.
- Compile a LaTeX report into a PDF format summarizing the setup and convergence analysis.

### 2. Run a single case manually

```bash
cd 0.1mm_dt0.05s
/path/to/gpyro/executable Reference_CC.data
```


## Folder Structure

| File/Folder                      | Description                                                              |
|----------------------------------|--------------------------------------------------------------------------|
| `run.sh`                         | Script to run all simulations and generate the summary report.           |
| `plot_results.py`                | Python script for computing error norms and plotting convergence.        |
| `time_step_convergence.tex`      | LaTeX source of the convergence report.                                  |
| `time_step_convergence.pdf`      | Compiled report with convergence plots and numerical analysis.           |
| `*/Reference_CC.data`            | Gpyro input files for each case with a specific time step.               |
| `0.1mm_dt0.0005s/`               | Simulation case with finest temporal resolution (already run due to cost).|
| `0.1mm_dt10s/`, ..., `0.1mm_dt0.005s/` | Other time step refinement cases with moderate computational cost. |








