# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 19.1 / L003 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\outputs\TRKB_40chol_L003_19.1_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\outputs\TRKB_40chol_L003_19.1_10ps_R2_state.xml`
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
- Temperature K initial/final/min/max/mean: 311.508/311.674/307.033/313.449/310.366
- Potential energy kJ/mol initial/final/mean: -444766.907041/-446305.203428/-444830.001738

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 4.93712/4.0159/3.6724/6.00497/4.56184
- VAL17 distance A initial/final/min/max/mean: 4.07357/3.78361/3.37564/4.6959/3.92164
- SER20 distance A initial/final/min/max/mean: 3.63606/3.76404/2.93497/4.29236/3.55351
- Key-residue min distance A initial/final/min/max/mean: 3.63606/3.76404/2.93497/4.18992/3.52274
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.62992/2.51241/2.51241

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
