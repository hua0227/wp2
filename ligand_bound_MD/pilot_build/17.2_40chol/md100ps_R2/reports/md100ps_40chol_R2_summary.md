# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 17.2 / L010 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_40chol\md50ps_R2\outputs\TRKB_40chol_L010_17.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_40chol\md50ps_R2\outputs\TRKB_40chol_L010_17.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_40chol\md50ps_R2\outputs\TRKB_40chol_L010_17.2_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 311.464/309.123/308.192/313.012/310.417
- Potential energy kJ/mol initial/final/mean: -446616.420979/-449384.705159/-447306.884328

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.95627/3.51278/3.26222/4.32584/3.68892
- VAL17 distance A initial/final/min/max/mean: 3.7613/4.31575/3.48308/5.27464/4.16168
- SER20 distance A initial/final/min/max/mean: 3.75924/3.6854/3.03222/6.18016/4.66625
- Key-residue min distance A initial/final/min/max/mean: 3.75924/3.51278/3.03222/4.25953/3.62771
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.554685/0.963409/2.13564
- Ligand RMSD vs original R2 minimized input A initial/final/max: 1.35668/1.5255/3.17244

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
