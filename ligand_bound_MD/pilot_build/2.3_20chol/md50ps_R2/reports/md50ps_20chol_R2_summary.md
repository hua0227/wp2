# 20chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 2.3 / L006 / 20chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\outputs\TRKB_20chol_L006_2.3_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\outputs\TRKB_20chol_L006_2.3_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 308.334/313.008/307.332/313.542/310.313
- Potential energy kJ/mol initial/final/mean: -453472.976237/-453917.650554/-453688.023977

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.63026/3.57358/3.25437/4.03414/3.63188
- VAL17 distance A initial/final/min/max/mean: 3.63255/4.02495/3.38729/4.96114/4.02767
- SER20 distance A initial/final/min/max/mean: 6.08891/5.1281/4.17543/7.85034/6.50516
- Key-residue min distance A initial/final/min/max/mean: 3.63026/3.57358/3.25437/4.03414/3.61433
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.518204/1.9349/2.12378

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
