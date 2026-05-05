# 20chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 9.2 / L004 / 20chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md50ps_R2\outputs\TRKB_20chol_L004_9.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md50ps_R2\outputs\TRKB_20chol_L004_9.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md50ps_R2\outputs\TRKB_20chol_L004_9.2_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 309.05/309.277/306.879/313.012/310.176
- Potential energy kJ/mol initial/final/mean: -457283.561736/-457908.861541/-456892.107107

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 4.62984/4.08762/3.22645/5.2729/3.90034
- VAL17 distance A initial/final/min/max/mean: 3.68632/3.69957/3.248/4.9946/3.95614
- SER20 distance A initial/final/min/max/mean: 3.46316/4.17884/3.1388/6.48965/5.12799
- Key-residue min distance A initial/final/min/max/mean: 3.46316/3.69957/3.1388/4.43629/3.6944
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 1.0577/1.00458/2.63139
- Ligand RMSD vs original R2 minimized input A initial/final/max: 4.39515/3.40547/4.7255

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
