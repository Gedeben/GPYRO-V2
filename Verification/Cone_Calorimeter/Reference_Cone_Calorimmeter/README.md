# Cone Calorimeter 1D Validation Case

## Overview

This validation scenario reproduces a one-dimensional (1D) version of a cone calorimeter experiment. It is designed to compare the results of three fire pyrolysis codes: **ThermaKin**, **Gpyro**, and **FDS**. The input conditions are chosen such that the three models should yield theoretically identical results under ideal assumptions.

## How to Run

### 1. Run the full pipeline (simulation + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```

This command will:
- Launch the Gpyro simulation using the provided input.
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
| `Reference_CC.data`     | Input file for Gpyro simulation.                                            |
| `reference_cc.tex`      | LaTeX source file for the report.                                           |
| `reference_cc.bib`      | Bibliography file for the report.                                           |
| `reference_cc.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |
| `FDS/`                  | Contains input and output files for the FDS simulation.                     |
| `ThermaKin/`            | Contains input and output files for the ThermaKin simulation.               |

