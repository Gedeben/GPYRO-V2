# Cone Calorimeter 2D/3D Extension Case

## Overview

This validation scenario extends the reference 1D cone calorimeter configuration to two- and three-dimensional (2D and 3D) domains in order to assess the consistency of the numerical solver in multi-dimensional setups.

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```

This command will:
- Run Gpyro simulations for the 1D, 2D, and 3D configurations.
- Perform post-processing to compare Mass Loss Rate (MLR) results.
- Compile the LaTeX report summarizing the methodology and results.

### 2. Run only the Gpyro simulation

```bash
/path/to/gpyro/executable <input_file.data>
```

Run this manually inside the `1D/`, `2D/`, or `3D/` folder, depending on the configuration to test.

## Folder Content

| File/Folder              | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| `1D/`                    | Input files for the 1D cone calorimeter simulation.                         |
| `2D/`                    | Input files for the 2D simulation.                                          |
| `3D/`                    | Input files for the 3D simulation.                                          |
| `plot_results.py`        | Python script for plotting and comparing MLR between 1D, 2D, and 3D cases.  |
| `ref_cc_2D_3D_extention.tex` | LaTeX source file describing the simulation setup and results.          |
| `ref_cc_2D_3D_extention.pdf` | Precompiled version of the LaTeX report.                                |
| `run.sh`                 | Bash script to automate simulation, post-processing, and report generation. |


