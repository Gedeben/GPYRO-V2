# 1D Convective Cooling Verification Case

## Overview
This case examines one-dimensional heat transfer through a slab that's 1 meter thick. The initial temperature of the plate is 1000°C, and it's exposed to the air at 0°C. The back is insulated and there's no radiation from the surface.

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
| `convective_cooling.data`     | Input file for Gpyro simulation.                                            |
| `convective_cooling.tex`      | LaTeX source file for the report.                                           |
| `convective_cooling.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

