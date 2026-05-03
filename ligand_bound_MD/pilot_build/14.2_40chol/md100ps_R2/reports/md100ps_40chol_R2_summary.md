# 40chol 100 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 100 ps restrained MD pilot continuation for 14.2 / L002 / 40chol after the successful 50 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md50ps_R2\outputs\TRKB_40chol_L002_14.2_50ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md50ps_R2\outputs\TRKB_40chol_L002_14.2_50ps_R2_state.xml`
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

- Completed 100 ps: True
- Last step: 200000/200000
- Stop reason: none
- NaN detected: False
- Potential energy finite across monitor rows: True
- Temperature K initial/final/min/max/mean: 311.039/311.505/307.013/314.195/310.457
- Potential energy kJ/mol initial/final/mean: -446508/-447460/-446411

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.82823/4.47124/3.06787/5.5268/3.80623
- VAL17 distance A initial/final/min/max/mean: 7.75557/8.74193/7.27015/9.06403/8.25231
- SER20 distance A initial/final/min/max/mean: 6.50974/7.98056/6.50974/9.16933/8.12394
- Key-residue min distance A initial/final/min/max/mean: 3.82823/4.47124/3.06787/5.5268/3.80623
- Ligand near chains always 2: True
- Ligand RMSD vs 100 ps input A initial/final/max: 0.000487925/1.55302/3.00028

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend 14.2 20chol vs 40chol comparison analysis: True
- This stage is still restrained MD pilot work, not production MD.
