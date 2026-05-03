# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 6.2 / L008 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_40chol\outputs\TRKB_40chol_L008_6.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_40chol\outputs\TRKB_40chol_L008_6.2_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 307.58/310.336/307.562/313.874/310.141
- Potential energy kJ/mol initial/final/mean: -443740.8343/-443907.888988/-444981.417777

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.5946/3.68427/3.03612/4.44415/3.56691
- VAL17 distance A initial/final/min/max/mean: 4.10637/5.01022/3.55331/6.79328/5.12374
- SER20 distance A initial/final/min/max/mean: 5.32724/5.33035/4.1925/7.86207/6.05212
- Key-residue min distance A initial/final/min/max/mean: 3.5946/3.68427/3.03612/4.44415/3.56365
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.519507/1.7245/2.73923

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
