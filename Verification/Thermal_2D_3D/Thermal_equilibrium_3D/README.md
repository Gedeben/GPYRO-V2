# 3D Thermal Equilibrium Verification Case

## Overview

This verification case simulates a 3D thermal equilibrium scenario to validate the conservative behaviour of the 3D thermal solver.

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
/path/to/gpyro/executable Thermal_equilibrium.data
```
---
## Folder Content
| File/Folder                  | Description                                                                  |
| ---------------------------- | ---------------------------------------------------------------------------- |
| `Thermal_equilibrium.data`   | Input file for Gpyro simulation.                                             |
| `thermal_equilibrium_01.ini` | Smokeview initialization file (viewpoints, coloring, etc.).                  |
| `thermal_equilibrium_01.ssf` | Smokeview script for automated image generation.                             |
| `plot_results.py`            | Python script for extracting and visualizing simulation results.             |
| `thermal_equilibrium_3D.tex` | LaTeX source for the PDF report.                                             |
| `thermal_equilibrium_3D.pdf` | Precompiled report summarizing the simulation.                               |
| `run.sh`                     | Bash script to execute the complete simulation and post-processing pipeline. |


---

