# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 6.2 / L008 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- 50 ps checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_40chol\md50ps_R2\outputs\TRKB_40chol_L008_6.2_50ps_R2.chk`
- 50 ps state XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_40chol\md50ps_R2\outputs\TRKB_40chol_L008_6.2_50ps_R2_state.xml`
- 50 ps final PDB input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\6.2_40chol\md50ps_R2\outputs\TRKB_40chol_L008_6.2_50ps_R2_final.pdb`
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
- Temperature K initial/final/min/max/mean: 310.025/308.83/307.357/314.443/310.43
- Potential energy kJ/mol initial/final/mean: -445010.143871/-447549.047191/-447061.571273

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.74405/3.49956/3.09964/4.06799/3.56118
- VAL17 distance A initial/final/min/max/mean: 5.5132/6.31934/5.30606/7.09401/6.25162
- SER20 distance A initial/final/min/max/mean: 6.31445/6.5774/5.51066/7.74437/6.47023
- Key-residue min distance A initial/final/min/max/mean: 3.74405/3.49956/3.09964/4.06799/3.56118
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.540049/0.982462/1.92747
- Ligand RMSD vs original R2 minimized input A initial/final/max: 2.92386/3.43107/4.03875

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
