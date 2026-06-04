# Fixed Boundary Temperatures 1D Verification Case

## Overview

This verification case represents one-dimensional (1D) heat transfer inside a slab with fixed temperatures at boundary conditions. The results are compared with an analytical solution from Incropera et al. 

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
/path/to/gpyro/executable Reference_CC.data
```

This executes the Gpyro model using the provided input file without any post-processing or reporting.

## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `fixed_temp.data`     | Input file for Gpyro simulation.                                            |
| `fixed_temp.tex`      | LaTeX source file for the report.                                           |
| `fixed_temp.bib`      | Bibliography file for the report.                                           |
| `fixed_temp.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

