# 1D Heat Conduction Verification Case

## Overview
Four distinct cases of one-dimensional transient heat conduction through a flat plate with a thickness of 0.1 meters are analyzed. In each case, the plate is initially at a uniform temperature of 20°C and is exposed on one side to ambient air at 120°C, while the other side is perfectly insulated. Radiative heat transfer is considered negligible in all scenarios.

## How to Run

### 1. Run the full pipeline (simulation + analytical solution + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```

This command will:
- Launch the Gpyro simulation using the provided input.
- Calculate the analytical solution.
- Perform automatic post-processing of the results.
- Generate a LaTeX report in PDF format summarizing the setup and results.

### 2. Run only the Gpyro simulation

```bash
/path/to/gpyro/executable radiation_1d.data
```

This executes the Gpyro model using the provided input file without any post-processing or reporting.

## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `heat_conduction_x/heat_conduction_x.data`     | Input file for Gpyro simulation.                                            |
| `heat_conduction.tex`      | LaTeX source file for the report.                                           |
| `heat_conduction.pdf`      | Precompiled report summarizing the case and results.                        |
| `heat_conduction_x/plot_results.py`       | Python script for post-processing and validation checks.                    |
| `heat_conduction_x/analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

