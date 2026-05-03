# 20chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 19.1 / L003 / 20chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\outputs\TRKB_20chol_L003_19.1_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\outputs\TRKB_20chol_L003_19.1_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 310.701/312.637/306.66/314.33/310.274
- Potential energy kJ/mol initial/final/mean: -454268.815949/-454017.376007/-453597.244787

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.83019/4.43781/3.30224/4.46045/3.76664
- VAL17 distance A initial/final/min/max/mean: 3.20948/3.28734/2.90202/4.25877/3.60637
- SER20 distance A initial/final/min/max/mean: 3.68595/3.21435/3.00849/4.24534/3.56861
- Key-residue min distance A initial/final/min/max/mean: 3.20948/3.21435/2.90202/3.86764/3.39745
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.576986/1.6911/1.74851

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
