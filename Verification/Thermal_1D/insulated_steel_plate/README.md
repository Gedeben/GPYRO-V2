# 1D Insulated Steel Plate Verification Case

## Overview
This case examines one-dimensional heat transfer through a multi-layered, insulated. The composite plate consists of 1 cm of steel covered by 2 cm of insulation on each side. The initial temperature of the plate is 20°C. The air's temperature is 480°C at the front of the slab, and 20°C at its back. Thermal radiation is neglected. The slab is heated for 10h to obtain the steady-state temperature profile. the convective heat transfer coefficient (h) is 10 W/m²K.


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
| `insulated_steel_plate.data`     | Input file for Gpyro simulation.                                            |
| `insulated_steel_plate.tex`      | LaTeX source file for the report.                                           |
| `insulated_steel_plate.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `analytical_solution.py`       | Python script to calculate the analytical solution.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

