# OMP Parallelization Verification

## Overview

This case evaluates the parallel performance of the cone calorimeter simulation by varying the number of OpenMP threads. 

---

## How to Run

### Full pipeline (simulations + post-processing + report generation)

```bash
bash run.sh /path/to/gpyro/executable
```
This command will:

- Run multiple Gpyro simulations with increasing numbers of OpenMP threads.
- Automatically post-process the results with the Python script.
- Generate a PDF report that summarizes the parallelization performance .

> **Note:**  The tested thread counts follow the predefined list (1, 2, 3, 4, 6, 8, 12, 14, 22, 30). Simulations are executed sequentially for each value in this list, until the number exceeds the maximum number of logical CPU cores detected on the system.

### Manual execution with a specific number of threads
```bash
export OMP_NUM_THREADS=4   # for exemple 4 Threads
/path/to/gpyro/executable ref_CC_NZ5000.data
```
---
## Folder Content

| File/Folder               | Description                                                            |
| ------------------------- | ---------------------------------------------------------------------- |
| `ref_CC_NZ5000.data`      | Input file for Gpyro simulation with 5000 cellls.                      |
| `plot_results.py`         | Python script to extract, process, and visualize performance data.     |
| `OMP_parallelization.tex` | LaTeX source file for the PDF report.                                  |
| `OMP_parallelization.pdf` | Precompiled report containing results and figures.                     |
| `run.sh`                  | Bash script to automate simulation and post-processing across threads. |



---


