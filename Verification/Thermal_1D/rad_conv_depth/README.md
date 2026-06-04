# 1D radiation in depth Validation Case

## Overview

This validation scenario represents a one-dimensional (1D) radiation heat transfer case, where multiple sub-cases are taken into consideration. In one case, we consider a slab's heating by radiation, with absorption into the depth and cooling by convection at the surface. We look at the temperature evolution and the results were validated by comparison with an analytical solution from Lautenberger's thesis dissertation.

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

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
/path/to/gpyro/executable 1d_conv_rad.data
```

This executes the Gpyro model using the provided input file without any post-processing or reporting.

## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `rad.data`     | Input file for Gpyro simulation.                                            |
| `rad.tex`      | LaTeX source file for the report.                                           |
| `rad.bib`      | Bibliography file for the report.                                           |
| `rad.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |
