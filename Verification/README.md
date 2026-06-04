# Validation and Verification Cases for GPYRO

This folder contains all **Validation** and **Verification** cases presented in the validation manual.
It provides the input files, automation scripts, and post-processing utilities.

---

## How to Run the Cases

### Run All Verification Cases
A bash script is provided to automatically run all verification cases, post-process results, and generate a PDF report. 

From this folder, run:
```bash
run_verification.sh /path/to/gpyro_executable
```
Here, /path/to/gpyro_executable is the executable version of Gpyro that you want to use.

### Run a Single Case
To run only one case, navigate to its folder and execute:
```
run.sh /path/to/gpyro_executable
```
Each folder may also include additional remarks in its local README.md.


## Visualization with Smokeview

The **3D reference cone calorimeter case with lateral convection** requires **Smokeview** for visualization.
When running the automated script, you will be asked whether you want to launch Smokeview:

- If you select **No**,the post processing will not be performed.
- If you select **Yes**, ensure that Smokeview is correctly installed and accessible in your environment. The post-processing and will then be executed.

### Installing Smokeview
You can download Smokeview from the official NIST page:
[https://pages.nist.gov/fds-smv/downloads.html]



## Repository Content

| File/Folder                | Description                                               |
|----------------------------|-----------------------------------------------------------|
| `clean.sh`                 | Automatic cleanup of simulation results                   |
| `run_verification.sh`      | Master script to run all verification cases               |
| `run_utils.sh`             | bash functions used in run.sh scripts                     |
| `verification_report.tex`  | LaTeX source of the verification report                   |
| `verification_report.pdf`  | Precompiled verifcation report                            |
| `Validation`               | Folder with validation cases                              |
| `Cone_Calorimeter`         | Main folder for cone calorimeter-based verification cases |


    




