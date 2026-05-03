# 50 ps Restrained MD Pilot R2 Summary

## Purpose

- Run a 50 ps restrained MD pilot for only 14.2 / L002 / 20chol from the successful 10 ps R2 endpoint.
- Continue the R2 PBC-aware restraint strategy and evaluate whether the ligand pose remains retained before any longer restrained pilot.
- This is not production MD.

## Inputs

- PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\targeted_clash_cleanup\outputs\branch_B_remove_POPC94\branch_B_cleaned.psf`
- Input PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\outputs\TRKB_20chol_L002_14.2_shortMD_R2_final.pdb`
- GRO box: `C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro`
- CHARMM toppar: `C:\TRKB_WP2\TRKB_20CHOL\toppar`
- Ligand stream: `C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str`

## R2 Restraint Strategy

- Protein backbone: medium PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.
- Ligand heavy atoms: weak PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.
- Lipid/cholesterol heavy atoms: extremely weak PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.
- Water/ions: unrestrained.
- Force constants: protein backbone 500, ligand heavy 25, lipid/chol heavy 5 kJ/mol/nm^2.
- The old absolute Cartesian lipid/cholesterol restraint expression was not used.

## Settings

- OpenMM platform used: OpenCL
- CUDA used: NO
- Nonbonded method: PME
- Nonbonded cutoff: 1.2 nm
- constraints: HBonds
- rigidWater: True
- Timestep: 0.500 fs
- Integrator: LangevinMiddleIntegrator
- Temperature: 310 K
- Total planned time: 50.000 ps
- Monitor interval: 500 steps
- DCD interval: 1000 steps

## Outcome

- Completed 50 ps: YES
- Last completed step: 100000 / 100000
- Last completed time: 50.000 ps
- Stop reason: none
- NaN/nonfinite detected: NO
- Temperature stable below 500 K: YES
- Temperature range / average / final: 307.18 to 313.589 K / 310.313 K / 309.929 K
- Potential energy initial/final/min/max/average: -451823 / -453981 / -454269 / -449971 / -452255 kJ/mol
- Potential energy finite through monitor rows: YES
- Potential energy final-minus-initial: -2157.97 kJ/mol
- Ligand min distance to TYR13/VAL17/SER20 initial/final/min/max/average: 3.82059 / 3.70177 / 2.90834 / 4.50632 / 3.62435 A
- Ligand near chains always at least 2: YES
- Ligand heavy atom RMSD vs 50 ps input pose initial/final/max: 5.60855e-06 / 1.26568 / 1.7428 A
- Severe clashes appeared: NO
- Close contacts appeared: NO
- Final severe/close contacts: 0 / 0
- Recommend next stage: YES
- Recommended next stage: The monitored criteria support a non-production 100 ps restrained MD pilot with the same R2 PBC-aware restraints; 40chol / other-ligand expansion can be considered after that pilot remains stable.

## Output Files

- Final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2_final.pdb`
- DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2.dcd`
- State XML: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2_state.xml`
- Checkpoint: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2.chk`
- Monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\reports\md50ps_R2_monitor.csv`
- Log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\logs\md50ps_R2.log`

This stage is a 50 ps restrained MD pilot only, not production MD.

## Next Step

- The monitored criteria support a non-production 100 ps restrained MD pilot with the same R2 PBC-aware restraints; 40chol / other-ligand expansion can be considered after that pilot remains stable.
