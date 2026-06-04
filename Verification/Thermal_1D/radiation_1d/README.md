# 1D Radiation Verification Case

## Overview

This verification case examines the one-dimensional (1D) heating of a slab, with radiation as the only heat source and no heat losses. The simulation results are validated against the analytical solution for a semi-infinite slab, as presented in Incropera et al.

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
| `radiation_1d.data`     | Input file for Gpyro simulation.                                            |
| `radiation_1d.tex`      | LaTeX source file for the report.                                           |
| `radiation_1d.bib`      | Bibliography file for the report.                                           |
| `radiation_1d.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

