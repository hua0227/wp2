# Short Restrained MD Retry R2 Summary

## Inputs

- PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\targeted_clash_cleanup\outputs\branch_B_remove_POPC94\branch_B_cleaned.psf`
- R2 minimized PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\restraint_strategy_fix\outputs\R2_pbcaware_minimized.pdb`
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
- Schedule: 50 K 2 ps, 100 K 2 ps, 200 K 2 ps, 310 K 4 ps.

## Outcome

- Completed 10 ps: YES
- Last completed step: 20000 / 20000
- Stop reason: none
- NaN/nonfinite detected: NO
- Temperature stable below 500 K: YES
- Maximum sampled temperature: 314.215 K
- Final potential energy: -451125.874752 kJ/mol
- Ligand min distance to TYR13/VAL17/SER20: 3.82108 A
- Ligand near two chains: YES
- Ligand heavy atom RMSD vs R2 input pose: 0.914187 A
- Final severe/close contacts: 0 / 0
- Recommend next 50 ps restrained MD pilot: YES

## Output Files

- Final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\outputs\TRKB_20chol_L002_14.2_shortMD_R2_final.pdb`
- DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\outputs\TRKB_20chol_L002_14.2_shortMD_R2.dcd`
- Monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\reports\short_md_retry_R2_monitor.csv`
- Log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\short_md_retry_R2\logs\short_md_retry_R2.log`

This stage is a short restrained MD pilot retry only, not production MD.

## Next Step

- The monitored criteria support a 50 ps restrained MD pilot reattempt using the same R2 PBC-aware restraints; this should remain non-production.
