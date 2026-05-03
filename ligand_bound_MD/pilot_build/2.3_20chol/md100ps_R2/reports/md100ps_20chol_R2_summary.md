# 20chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 2.3 / L006 / 20chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md50ps_R2\outputs\TRKB_20chol_L006_2.3_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md50ps_R2\outputs\TRKB_20chol_L006_2.3_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md50ps_R2\outputs\TRKB_20chol_L006_2.3_50ps_R2_final.pdb`
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
- Monitoring cadence: 1 ps
- DCD cadence: 2 ps
- This is restrained MD pilot work, not production MD.

## Outcome

- Completed 100 ps: True
- Last step: 200000/200000
- Stop reason: none
- NaN detected: False
- Temperature K initial/final/min/max/mean: 310.536/310.976/307.158/314.423/310.341
- Potential energy kJ/mol initial/final/mean: -453681.960612/-457109.896159/-456174.70427

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.95805/3.63059/3.14915/4.5542/3.72437
- VAL17 distance A initial/final/min/max/mean: 4.19474/3.97804/3.4699/4.88547/4.08797
- SER20 distance A initial/final/min/max/mean: 5.23463/4.69417/3.77329/7.47164/5.80177
- Key-residue min distance A initial/final/min/max/mean: 3.95805/3.63059/3.14915/4.24144/3.69646
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.675177/1.2911/2.20578
- Ligand RMSD vs original R2 minimized input A initial/final/max: 1.623/1.67641/2.62674

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
