# 3D Lateral Convection Validation Case

## Overview

This validation scenario extends the 1D reference cone calorimeter setup into three dimensions by enabling natural convection on lateral surfaces. The goal is to evaluate Gpyro's ability to simulate multidimensional pyrolysis front propagation.

---

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```

This command will:
- Launch the Gpyro simulation using the provided input file.
- Optionally use Smokeview to generate visual outputs.
- Automatically post-process the results with the Python script.
- Generate a PDF report describing the setup and key results.

> **Note:** You will be prompted to define the path to **Smokeview** if it is not already set. You can choose to skip this step if you do not wish to generate graphical outputs.

### 2. Run only the Gpyro simulation

```bash
/path/to/gpyro/executable Reference_CC_3D_Conv.data
```
---

## Folder Content

| File/Folder                     | Description                                                                 |
|--------------------------------|-----------------------------------------------------------------------------|
| `Reference_CC_3D_Conv.data`    | Input file for Gpyro simulation.                                            |
| `reference_cc_3D_conv_01.ini`  | Smokeview initialization file (viewpoints, coloring, etc.).                |
| `reference_cc_3D_conv_01.ssf`  | Smokeview script for automating image generation.                          |
| `plot_results.py`              | Python script for extracting and visualizing simulation results.            |
| `ref_cc_3D_conv.tex`           | LaTeX source for the PDF report.                                           |
| `ref_cc_3D_conv.pdf`           | Precompiled report summarizing the simulation.                              |
| `run.sh`                       | Bash script to execute the complete simulation and post-processing pipeline.|


---



