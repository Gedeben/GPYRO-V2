# Validation and Verification Cases for GPYRO

This folder contains all **Validation** and **Verification** cases presented in the *SoftwareX* article.
It provides the input files, automation scripts, and post-processing utilities.

---

## Available Cases

| Folder                                          | Description                                           |
|-------------------------------------------------|-------------------------------------------------------|
| `./Cone_Calorimeter/Reference_Cone_Calorimeter` | 1D reference cone calorimeter case                    |
| `./Cone_Calorimeter/Time_Step_Convergence`      | Time-step convergence study                           |
| `./Cone_Calorimeter/Spatial_Convergence`        | Spatial mesh convergence study                        |
| `./Cone_Calorimeter/Extention_Cone_Calo_2D_3D`  | Extension of the reference case to 2D/3D              |
| `./Cone_Calorimeter/3D_Lateral_convection`      | 3D cone calorimeter case including lateral convection |
| `./Validation/MaCFP_PMMA`                       | MaCFP PMMA cone calorimeter validation case           |

Each case folder contains its own `README.md` file with additional notes and specific instructions.

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


## Comparing GPYRO Versions

To compare the legacy version of GPYRO (V0.8) with the new version (v2) presented in this repository, simply run the cases twice, specifying different executables, and compare the result:

```
run_verification.sh /path/to/legacy/gpyro      # v0.8
run_verification.sh /path/to/new/gpyro         # v2
```

For that you need the executable of both version 

### Get Gpyro V0.8 executable 

The legacy version can be obtained directly from GitHub:
https://github.com/lautenberger/gpyro

### Get Gpyro V2 executable 

You need to compile the code on a **Linux** system. Two build scripts are provided in build/linux/:
make_gnu.sh → GNU Fortran compiler
make_intel.sh → Intel Fortran compiler
Both scripts compile release and debug versions of gpyro and gpyro_propest.

Example (GNU Fortran):
```
cd build/linux
bash make_gnu.sh
```
Executables will be located in:

build/linux/gpyro/
build/linux/gpyro_propest/

*If compilation fails, check your compiler installation and update paths inside the scripts.*


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


    




