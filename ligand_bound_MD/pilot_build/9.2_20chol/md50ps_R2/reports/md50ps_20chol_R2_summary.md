# 20chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 9.2 / L004 / 20chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\outputs\TRKB_20chol_L004_9.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\outputs\TRKB_20chol_L004_9.2_10ps_R2_state.xml`
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
- This is restrained MD pilot work, not production MD.

## Outcome

- Completed 50 ps: True
- Last step: 100000/100000
- Stop reason: none
- NaN detected: False
- Temperature K initial/final/min/max/mean: 311.054/310.307/306.756/313.483/310.221
- Potential energy kJ/mol initial/final/mean: -453148.317107/-455190.30783/-454567.123675

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.33938/3.93064/3.17503/4.88694/3.82478
- VAL17 distance A initial/final/min/max/mean: 4.21601/3.7739/3.34419/4.90448/4.05175
- SER20 distance A initial/final/min/max/mean: 3.92381/3.85862/3.30532/6.53237/4.41549
- Key-residue min distance A initial/final/min/max/mean: 3.33938/3.7739/3.17503/4.85873/3.70768
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.620255/1.94916/3.30601

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
