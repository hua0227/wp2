# 40chol 10 ps R2 Restrained MD Pilot Summary

## Purpose

Run a 10 ps restrained MD pilot for 14.2 / L002 / 40chol after topology/read test, targeted cleanup, R2 minimization, and very-short probe.

## Inputs

- PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound_cleaned.psf`
- Input PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_R2_minimized.pdb`
- GRO box: `C:\TRKB_WP2\TRKB_40CHOL\gromacs\step5_input.gro`
- CHARMM toppar: `C:\TRKB_WP2\TRKB_40CHOL\toppar`
- Ligand str: `C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str`

## Method

- OpenCL used: True
- Platform: OpenCL
- Timestep: 0.0005 ps (0.5 fs)
- Integrator: LangevinMiddleIntegrator
- Nonbonded method: PME
- Constraints: HBonds
- rigidWater: True
- Temperature schedule: 50 K 2 ps, 100 K 2 ps, 200 K 2 ps, 310 K 4 ps
- R2 restraints: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.
- PBC-aware restraint expression: `0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2`

## Outcome

- Completed 10 ps: True
- Last step: 20000/20000
- Stop reason: none
- NaN detected: False
- Potential energy finite across monitor rows: True
- Temperature K initial/final/min/max/mean: 50.2156/309.543/34.2285/313.042/191.015
- Potential energy kJ/mol initial/final/min/max/mean: -615276/-444189/-615276/-443747/-519056

## Ligand Retention

- TYR13 distance A initial/final/min/max/mean: 3.60697/3.78204/3.03179/4.2273/3.63808
- VAL17 distance A initial/final/min/max/mean: 3.69188/4.80769/3.3749/4.86692/3.8653
- SER20 distance A initial/final/min/max/mean: 5.12268/6.65544/4.99965/7.15012/5.92529
- Key-residue min distance A initial/final/min/max/mean: 3.60697/3.78204/3.03179/4.21873/3.61544
- Ligand near chains always 2: True
- Ligand heavy-atom RMSD A initial/final/max: 5.24844e-06/1.21958/1.42974

## Clash Monitoring

- Severe clashes appeared: False
- Close contacts appeared: False
- Final severe/close contacts: 0/0

## Recommendation

- Recommend next 50 ps 40chol restrained MD pilot: True
- This stage is still restrained MD pilot work, not production MD.
