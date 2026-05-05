# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 9.2 / L004 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md50ps_R2\outputs\TRKB_40chol_L004_9.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md50ps_R2\outputs\TRKB_40chol_L004_9.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md50ps_R2\outputs\TRKB_40chol_L004_9.2_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 310.617/311.039/306.928/314.071/310.452
- Potential energy kJ/mol initial/final/mean: -446610.935543/-446521.780758/-447385.344357

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.0016/3.3584/3.0016/4.07209/3.49724
- VAL17 distance A initial/final/min/max/mean: 4.14722/6.19057/3.43751/6.94213/4.83366
- SER20 distance A initial/final/min/max/mean: 4.05344/6.19966/3.28306/6.38796/4.83775
- Key-residue min distance A initial/final/min/max/mean: 3.0016/3.3584/3.0016/4.07209/3.4888
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.600848/1.11044/1.57828
- Ligand RMSD vs original R2 minimized input A initial/final/max: 1.53797/2.15594/2.48472

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
