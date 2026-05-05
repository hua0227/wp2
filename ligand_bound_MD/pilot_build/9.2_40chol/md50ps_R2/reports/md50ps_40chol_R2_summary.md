# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 9.2 / L004 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\outputs\TRKB_40chol_L004_9.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\outputs\TRKB_40chol_L004_9.2_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 312.811/309.114/306.994/312.903/310.541
- Potential energy kJ/mol initial/final/mean: -446011.336422/-447030.062985/-445250.117438

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.34044/3.09142/2.91965/4.02204/3.50982
- VAL17 distance A initial/final/min/max/mean: 4.10029/4.37571/3.49752/6.5078/4.71055
- SER20 distance A initial/final/min/max/mean: 5.20166/4.63644/3.59968/6.11746/5.05184
- Key-residue min distance A initial/final/min/max/mean: 3.34044/3.09142/2.91965/4.0182/3.49868
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.442889/0.627492/2.04181

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
