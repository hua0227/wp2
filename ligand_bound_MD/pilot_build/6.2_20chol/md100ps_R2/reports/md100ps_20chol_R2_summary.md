# 20chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 6.2 / L008 / 20chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_20chol\md50ps_R2\outputs\TRKB_20chol_L008_6.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_20chol\md50ps_R2\outputs\TRKB_20chol_L008_6.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_20chol\md50ps_R2\outputs\TRKB_20chol_L008_6.2_50ps_R2_final.pdb`
- Fallback used: False
- Fallback reason: none

## Method

- OpenCL used: True
- Platform: OpenCL
- Timestep: 0.0005 ps (0.5 fs)
- Integrator: LangevinMiddleIntegrator
- Nonbonded method: PME
- Constraints: HBonds
- rigidWater: True
- Temperature: 310 K
- R2 restraints: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.
- PBC-aware restraint expression: `0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2`
- Monitoring cadence: 1 ps
- DCD cadence: 2 ps
- This is restrained MD pilot work, not production MD.

## Outcome

- Completed 100 ps: True
- Last step: 200000/200000
- Stop reason: none
- NaN detected: False
- Temperature K initial/final/min/max/mean: 312.913/309.611/307.972/313.095/310.344
- Potential energy kJ/mol initial/final/mean: -455924.522993/-455392.786177/-456347.711055

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 4.14219/3.57923/3.27508/5.03961/3.88429
- VAL17 distance A initial/final/min/max/mean: 5.94911/6.27969/3.77848/6.43541/4.84524
- SER20 distance A initial/final/min/max/mean: 4.03959/4.12865/3.20376/4.96546/3.75247
- Key-residue min distance A initial/final/min/max/mean: 4.03959/3.57923/3.20376/4.30441/3.60377
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.654458/1.54022/2.27282
- Ligand RMSD vs original R2 minimized input A initial/final/max: 3.56056/4.50161/5.18175

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
