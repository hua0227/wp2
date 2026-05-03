# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 17.2 / L010 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_40chol\outputs\TRKB_40chol_L010_17.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\17.2_40chol\outputs\TRKB_40chol_L010_17.2_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 310.797/310.128/307.772/312.881/310.327
- Potential energy kJ/mol initial/final/mean: -445148.302815/-444594.019612/-445337.451907

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.89904/3.99387/3.2524/5.2813/3.89904
- VAL17 distance A initial/final/min/max/mean: 3.8814/4.09961/3.30463/4.98426/4.08147
- SER20 distance A initial/final/min/max/mean: 4.74028/4.03725/3.07043/6.25363/4.36623
- Key-residue min distance A initial/final/min/max/mean: 3.8814/3.99387/3.07043/4.3124/3.75782
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.435182/0.75275/2.28052

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
