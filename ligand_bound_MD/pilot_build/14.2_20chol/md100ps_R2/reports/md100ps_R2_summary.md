# 100 ps Restrained MD Pilot R2 Summary

## Purpose

- Run a 100 ps restrained MD pilot continuation for only 14.2 / L002 / 20chol from the successful 50 ps R2 endpoint.
- Continue the R2 PBC-aware restraint strategy and evaluate whether the ligand pose remains retained before any system expansion.
- This is not production MD.

## Inputs

- PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\targeted_clash_cleanup\outputs\branch_B_remove_POPC94\branch_B_cleaned.psf`
- 50 ps final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2_final.pdb`
- 50 ps checkpoint: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2.chk`
- 50 ps state XML: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md50ps_R2\outputs\TRKB_20chol_L002_14.2_50ps_R2_state.xml`
- Checkpoint-compatible restraint anchor PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\outputs\TRKB_20chol_L002_14.2_shortMD_R2_final.pdb`
- Original R2 input PDB for RMSD reference: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\restraint_strategy_fix\outputs\R2_pbcaware_minimized.pdb`
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
- Total planned time: 100.000 ps
- Monitor interval: 2000 steps
- DCD interval: 4000 steps
- Started from checkpoint: YES
- Start mode: checkpoint
- Start fallback reason: none

## Outcome

- Completed 100 ps: YES
- Last completed step: 200000 / 200000
- Last completed time: 100.000 ps
- Stop reason: none
- NaN/nonfinite detected: NO
- Temperature stable below 500 K: YES
- Temperature range / average / final: 306.937 to 313.268 K / 310.414 K / 310.952 K
- Potential energy initial/final/min/max/average: -453981 / -453165 / -454421 / -450740 / -452422 kJ/mol
- Potential energy finite through monitor rows: YES
- Potential energy final-minus-initial: 815.945 kJ/mol
- Ligand min distance to TYR13/VAL17/SER20 initial/final/min/max/average: 3.70177 / 3.6668 / 3.24725 / 4.53247 / 3.74823 A
- Ligand near chains always at least 2: YES
- Ligand heavy atom RMSD vs 100 ps input pose initial/final/max: 0.000476316 / 1.19815 / 1.52269 A
- Ligand heavy atom RMSD vs original R2 input pose initial/final/max: 1.17115 / 1.56563 / 2.00164 A
- Severe clashes appeared: NO
- Close contacts appeared: NO
- Final severe/close contacts: 0 / 0
- Recommend next stage: YES
- Recommended next stage: The monitored criteria support considering a separate non-production 14.2 / 40chol pilot with the same Branch-B/R2 strategy; do not treat this as production MD.

## Output Files

- Final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2_final.pdb`
- DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2.dcd`
- State XML: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2_state.xml`
- Checkpoint: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2.chk`
- Monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\reports\md100ps_R2_monitor.csv`
- Log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\logs\md100ps_R2.log`

This stage is a 100 ps restrained MD pilot only, not production MD.

## Next Step

- The monitored criteria support considering a separate non-production 14.2 / 40chol pilot with the same Branch-B/R2 strategy; do not treat this as production MD.
