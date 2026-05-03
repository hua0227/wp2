# 40chol 50 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 50 ps restrained MD pilot for 14.2 / L002 / 40chol after the successful 10 ps R2 restrained MD pilot.

## Continuation

- Continued from checkpoint: True
- Checkpoint input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\short_md_R2\outputs\TRKB_40chol_L002_14.2_10ps_R2.chk`
- State XML input: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\short_md_R2\outputs\TRKB_40chol_L002_14.2_10ps_R2_state.xml`
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
- Potential energy finite across monitor rows: True
- Temperature K initial/final/min/max/mean: 309.543/311.039/305.997/313.026/310.151
- Potential energy kJ/mol initial/final/mean: -444541/-445546/-444705

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.78204/3.82823/3.12865/4.10192/3.57109
- VAL17 distance A initial/final/min/max/mean: 4.80769/7.75557/4.57741/8.75018/7.41172
- SER20 distance A initial/final/min/max/mean: 6.65544/6.50974/6.05695/9.27324/7.63254
- Key-residue min distance A initial/final/min/max/mean: 3.78204/3.82823/3.12865/4.10192/3.57109
- Ligand near chains always 2: True
- Ligand RMSD vs 50 ps input A initial/final/max: 0.000440667/2.46876/2.94562

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 100 ps 40chol restrained MD: True
- This stage is still restrained MD pilot work, not production MD.
