# 20chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 17.2 / L010 / 20chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_20chol\outputs\TRKB_20chol_L010_17.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_20chol\outputs\TRKB_20chol_L010_17.2_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 312.187/311.485/307.926/314.185/310.456
- Potential energy kJ/mol initial/final/mean: -453940.711723/-454211.894828/-454011.470795

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.61295/3.72805/3.29197/6.93926/4.43896
- VAL17 distance A initial/final/min/max/mean: 3.85471/3.94055/3.15972/5.02023/3.97336
- SER20 distance A initial/final/min/max/mean: 3.92816/3.77915/2.8514/4.31877/3.59703
- Key-residue min distance A initial/final/min/max/mean: 3.61295/3.72805/2.8514/4.05715/3.53447
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.503547/0.875924/1.39511

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
