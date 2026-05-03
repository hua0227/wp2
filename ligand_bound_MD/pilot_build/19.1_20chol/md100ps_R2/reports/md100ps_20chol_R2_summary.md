# 20chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 19.1 / L003 / 20chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md50ps_R2\outputs\TRKB_20chol_L003_19.1_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md50ps_R2\outputs\TRKB_20chol_L003_19.1_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md50ps_R2\outputs\TRKB_20chol_L003_19.1_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 310.968/312.407/307.396/313.033/310.223
- Potential energy kJ/mol initial/final/mean: -454568.395539/-456382.600128/-456074.898864

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.70425/3.7754/3.33008/4.68348/3.87721
- VAL17 distance A initial/final/min/max/mean: 3.6503/3.19104/3.10259/4.80157/3.81703
- SER20 distance A initial/final/min/max/mean: 3.28552/3.65355/3.0656/4.36151/3.64137
- Key-residue min distance A initial/final/min/max/mean: 3.28552/3.19104/3.0656/4.10447/3.51387
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.986619/1.88172/2.48709
- Ligand RMSD vs original R2 minimized input A initial/final/max: 1.68814/2.14695/2.76888

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
