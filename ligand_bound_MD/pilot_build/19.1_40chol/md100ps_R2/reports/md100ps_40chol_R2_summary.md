# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 19.1 / L003 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md50ps_R2\outputs\TRKB_40chol_L003_19.1_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md50ps_R2\outputs\TRKB_40chol_L003_19.1_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md50ps_R2\outputs\TRKB_40chol_L003_19.1_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 308.366/311.199/307.428/313.78/310.107
- Potential energy kJ/mol initial/final/mean: -446846.796689/-447778.342099/-446917.410835

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 4.3977/3.7453/3.25597/5.44908/4.02567
- VAL17 distance A initial/final/min/max/mean: 4.30323/4.085/3.41715/5.56912/3.94914
- SER20 distance A initial/final/min/max/mean: 3.74522/3.79117/3.05314/4.54728/3.65769
- Key-residue min distance A initial/final/min/max/mean: 3.74522/3.7453/3.05314/4.19166/3.58712
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 1.09496/1.55146/2.44408
- Ligand RMSD vs original R2 minimized input A initial/final/max: 2.75423/2.01015/4.0735

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
