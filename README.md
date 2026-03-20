
# Gpyro V2 – A Generalized Pyrolysis Solver (0D → 3D)

                      (      (   '  )      (      (         )(    
                    ()/)   (()     ((       )    )/(      ()/(    
                 /█(▌▌))  (█(███)  (█)    (▌) (████(▖)  (/ ██\\)  
                |█(  __   |█|__)█|  (██\_/█)  |█(__)▌|  |█|  |█|  
                |█( |██|  |█ ███/    \███/    |█████/   |█|  |█|  
                |█(__|█|  |█|         |█|     |█| \█\   |█|__|█|  
                 \█████|  |█|          █|     |█|  \█\   \████/                                                  


**Gpyro V2** is a general-purpose pyrolysis modeling tool, capable of simulating thermal degradation of materials under external heat flux, from **0D lumped models to fully 3D geometries**.  

This repository is a fork of the original [Gpyro](https://github.com/lautenberger/gpyro) code developed by Lautenberger, with major refactoring, performance gains, and improved physics.

---

## What’s New in Gpyro V2?

### User Experience
- Complete **User Guide**.  
- Updated **Technical Reference Guide**.  
- Simplified input structure and clearer error messages.  

### Performance
- Full refactorization of the codebase for better maintainability.  
- Massive **OpenMP parallelization**.  
- New **non-recursive solid species solver and thermal solver**.  
- Speed-up ranging from 4× to 200×, particularly important in 3D cases.  
- Inverse simulation tool (**gpyro_propest**) is now practical on desktop machines.  

### Debugging and Physics
- Fixed **pore radiation** model contribution.  
- New implementation of in-depth radiation in 3D.  
- Stabilization of species conservation in near-extinction conditions.  
- Fixed inverse modeling tool for ATG analysis.  

### Verification and Validation
- New systematic verification procedure included.  
- Dedicated **Verification Guide**.  

---

## Repository Structure

**Language**: Gpyro is written in modern **Fortran 90**.  
**Executables**: The project compiles into three main applications:  

- **`gpyro`** → standalone pyrolysis solver.  
- **`gpyro_propest`** → inverse modeling tool for parameter estimation.  
- **`gpyro_fds`** → coupling of Gpyro with FDS v6.2.0 (release soon).  

**Directories**:  

- `/src` → Fortran source code.  
- `/Documentation` → User Guide, Technical Guide, Verification Guide.  
- `/build` → build instructions and precompiled binaries.  
- `/Verification` → automated verification procedures and benchmark cases.  

---

## Current Limitations
- Gas-phase solvers (pressure, diffusion, energy) are not fully functional.  
- Coupling with **FDS** is under finalization and will be released soon.  

---

## Disclaimer
This repository is a fork of [Gpyro](https://github.com/lautenberger/gpyro).
The FDS coupling is debugged but not included here.  

---



