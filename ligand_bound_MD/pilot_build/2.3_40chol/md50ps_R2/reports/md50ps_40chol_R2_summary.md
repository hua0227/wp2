# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 2.3 / L006 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\outputs\TRKB_40chol_L006_2.3_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\outputs\TRKB_40chol_L006_2.3_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 309.814/310.207/305.423/314.035/310.211
- Potential energy kJ/mol initial/final/mean: -445616.532357/-444901.319466/-445185.699999

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.52963/4.24819/3.03458/4.24819/3.50776
- VAL17 distance A initial/final/min/max/mean: 3.75926/3.87577/3.24215/4.79085/3.94395
- SER20 distance A initial/final/min/max/mean: 4.11084/5.59426/3.2239/5.59426/4.35305
- Key-residue min distance A initial/final/min/max/mean: 3.52963/3.87577/3.03458/4.18195/3.47633
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.627241/1.94715/1.94715

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
