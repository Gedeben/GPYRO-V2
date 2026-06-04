# Radiation Loss 1D Verification Case

## Overview

This verification scenario represents a one-dimensional (1D) case of heat loss by radiation. It is compared to results obtained by running the same case using FDS. Two subcases are ran, each with a different surrounding temperature.

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
| `Tgas0K/radiation_loss.data`     | Input file for Gpyro simulation of sub-case 1.                                            |
| `Tgas0K/FDS/`                  | Contains input and output files for the FDS simulation of sub-case 1.                     |
| `Tgas273K/radiation_loss.data`     | Input file for Gpyro simulation of sub-case 2.                                            |
| `Tgas273K/FDS/`                  | Contains input and output files for the FDS simulation of sub-case 2.                     |
| `radiation_loss.tex`      | LaTeX source file for the report.                                           |
| `radiation_loss.pdf`      | Precompiled report summarizing the case and results.                        |
| `plot_results.py`       | Python script for post-processing and validation checks.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

