# 20chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 17.2 / L010 / 20chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_20chol\md50ps_R2\outputs\TRKB_20chol_L010_17.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_20chol\md50ps_R2\outputs\TRKB_20chol_L010_17.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_20chol\md50ps_R2\outputs\TRKB_20chol_L010_17.2_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 310.215/308.375/307.219/313.529/310.141
- Potential energy kJ/mol initial/final/mean: -454911.741019/-456596.755668/-456290.068143

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.68113/4.24509/3.24076/5.70841/4.08645
- VAL17 distance A initial/final/min/max/mean: 3.52566/3.70658/2.87317/4.74454/3.87052
- SER20 distance A initial/final/min/max/mean: 3.69715/3.98709/3.18933/4.52134/3.73501
- Key-residue min distance A initial/final/min/max/mean: 3.52566/3.70658/2.87317/4.14936/3.60992
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.737751/0.75166/2.13821
- Ligand RMSD vs original R2 minimized input A initial/final/max: 0.580888/0.745881/2.66303

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
