# MaCFP PMMA Validation Case

## Overview

This validation case reproduces the MaCFP benchmark experiment for **poly(methyl methacrylate) (PMMA)** performed in the NIST Gasification Apparatus under an external radiant heat flux of **50 kW/m²**.

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
/path/to/gpyro/executable PMMA.data
```
This executes the Gpyro model using the provided input file without any post-processing or reporting.


## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `PMMA.data`             | Input file for Gpyro simulation.                                            |
| `PMMA.tex`              | LaTeX source file for the report.                                           |
| `PMMA.pdf`              | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

