# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 2.3 / L006 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md50ps_R2\outputs\TRKB_40chol_L006_2.3_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md50ps_R2\outputs\TRKB_40chol_L006_2.3_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md50ps_R2\outputs\TRKB_40chol_L006_2.3_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 308.556/313.032/307.393/313.232/310.55
- Potential energy kJ/mol initial/final/mean: -446038.149544/-449638.141732/-447528.98222

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.72347/3.60753/2.94422/4.12284/3.57987
- VAL17 distance A initial/final/min/max/mean: 4.33477/4.19607/3.34962/6.83733/4.72281
- SER20 distance A initial/final/min/max/mean: 5.95041/5.63985/4.86857/8.03557/6.48268
- Key-residue min distance A initial/final/min/max/mean: 3.72347/3.60753/2.94422/4.12284/3.56648
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 1.01087/1.59628/2.27603
- Ligand RMSD vs original R2 minimized input A initial/final/max: 3.56399/3.55301/4.46815

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
