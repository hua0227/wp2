# 20chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 6.2 / L008 / 20chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_20chol\outputs\TRKB_20chol_L008_6.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_20chol\outputs\TRKB_20chol_L008_6.2_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 311.074/312.71/306.567/313.617/310.178
- Potential energy kJ/mol initial/final/mean: -454356.174849/-453765.623091/-453868.335571

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.73986/4.44023/2.9692/4.71992/3.74823
- VAL17 distance A initial/final/min/max/mean: 4.9063/5.57951/3.64153/5.92294/4.57334
- SER20 distance A initial/final/min/max/mean: 4.07472/3.91111/3.08011/4.3788/3.6928
- Key-residue min distance A initial/final/min/max/mean: 3.73986/3.91111/2.9692/4.11489/3.57428
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.64098/1.67825/2.27392

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
