# Cone Calorimeter 1D Validation Case

## Overview

This validation scenario represents a one-dimensional (1D) mixed heat transfer case, where one side of the slab gains heat by radiation while simultaneously losing heat through re-radiation and convection, with the opposite side being adiabatic. It is designed to verify heat balance calculations under complex conditions. The results were validated by comparison with **FDS**.

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```

This command will:
- Launch the Gpyro simulation using the provided input.
- Perform automatic post-processing of the results and comparing them to the FDS results.
- Generate a LaTeX report in PDF format summarizing the setup and results.

### 2. Run only the Gpyro simulation

```bash
/path/to/gpyro/executable 1d_conv_rad.data
```

This executes the Gpyro model using the provided input file without any post-processing or reporting.

## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `1d_conv_rad.data`     | Input file for Gpyro simulation.                                            |
| `rad_conv.tex`      | LaTeX source file for the report.                                           |
| `rad_conv.bib`      | Bibliography file for the report.                                           |
| `rad_conv.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |
| `FDS/`                  | Contains input and output files for the FDS simulation.                     |
