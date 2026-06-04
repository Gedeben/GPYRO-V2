# Heat Conduction K and C variables 1D Verification Case

## Overview

This case examines one-dimensional heat transfer through a slab that's 0.01 meters thick. The initial temperature of the plate is 0°C, and it's exposed at its front surface to gas at 700°C. The thermal conductivity and specific heat capacity vary with the temperature, in the form of a ramp. The back is exposed to gas at 0°C. Thermal radiation is considered negligible.

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
| `heat_conduction_kc.data`     | Input file for Gpyro simulation.                                            |
| `heat_conduction_kc.tex`      | LaTeX source file for the report.                                           |
| `heat_conduction_kc.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `theoretical_results.csv`       | CSV file containing theoretical results obtained from FDS.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

