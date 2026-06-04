# Reaction Enthalpy Verification Case

## Overview
These three cases are used to verify the calculation of the reaction enthalpy in Gpyro. Each case verifies a type of reaction : endothermic, exothermic and sensible entalpy variation.


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
/path/to/gpyro/executable x_reaction.data
```

This executes the Gpyro model using the provided input file without any post-processing or reporting.

## Folder Content

| File/Folder             | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `x_reaction/x_reaction.data`     | Input file for Gpyro simulation.                                            |
| `Reaction_Enthalpy.tex`      | LaTeX source file for the report.                                           |
| `Reaction_Enthalpy.pdf`      | Precompiled report summarizing the case and results.                        |
| `x_reaction/plot_results.py`       | Python script for post-processing and validation checks.                    |
| `run.sh`                | Bash script to run the entire pipeline automatically.                       |

