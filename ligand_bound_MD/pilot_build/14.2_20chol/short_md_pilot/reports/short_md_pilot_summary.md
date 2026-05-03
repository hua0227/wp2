# Short Restrained MD Pilot Summary

## Purpose

This pilot tests whether the minimized ligand-bound 14.2 / L002 / 20chol system can tolerate a very short restrained NVT run while retaining the docked ligand near the expected TRKB pocket.

## Inputs

- Cleaned PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_bound_cleaned.psf`
- Minimized PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_minimized.pdb`
- GRO box: `C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro`
- CHARMM toppar: `C:\TRKB_WP2\TRKB_20CHOL\toppar`
- Ligand stream: `C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str`

## Settings

- OpenMM platform requested: OpenCL
- OpenMM platform used: OpenCL
- Timestep: 0.500 fs
- Nonbonded method: PME
- Nonbonded cutoff: 1.2 nm
- Constraints: HBonds
- rigidWater: True
- Protein backbone restraint: 1000.0 kJ/mol/nm^2
- Lipid/cholesterol heavy atom restraint: 1000.0 kJ/mol/nm^2
- Ligand heavy atom restraint: 150.0 kJ/mol/nm^2
- Water/ion restraints: none

## Outcome

- Completed 10 ps: NO
- Last completed step: 1 / 20000
- Last finite snapshot step: 1
- Safety stop reason: temperature_above_500K
- NaN/nonfinite detected: NO
- Maximum sampled temperature: 1797.49 K
- Temperature stable below 500 K: NO
- Potential energy remained finite: YES
- Last potential energy: 4.16529e+08 kJ/mol
- Last total energy: 4.17371e+08 kJ/mol
- Last finite potential energy: 4.16529e+08 kJ/mol
- Ligand min distance to TYR13/VAL17/SER20: 2.78641 A
- Ligand near two chains: YES
- Ligand heavy atom RMSD relative to minimized pose: 0.137304 A
- Recommend next longer restrained MD: NO

## Instability Diagnostics

- The official OpenCL run entered MD and stopped at step 1 because the sampled temperature reached 1797.49 K.
- The constrained initial system had potential energy `4.05834e+08 kJ/mol`; force diagnostics found the largest initial forces on LIG H9/H11 and nearby MEMB:94 POPC atoms, with max force about `2.10e+06 kJ/mol/nm`.
- Before HBond constraint projection, the minimized PDB had 4019 constrained bonds deviating by more than 0.01 nm, with maximum deviation about 0.2598 nm. After `applyConstraints`, the maximum constraint deviation dropped to about `9.44e-07 nm`, so the constraint solver itself is not the limiting issue.
- A diagnostic constrained restrained relaxation reduced max force to about `2.30e+04 kJ/mol/nm` and retained the ligand near the key residues and both chains, but a 0.5 fs test still exceeded 500 K by step 5. This supports further local relaxation/stabilization before any longer restrained MD.

## Final Key-Residue Distances

| Key residue | Segid | Resid | Min distance A | Within 6 A |
|---|---|---|---:|---:|
| TYR13 | PROA | 434 | 7.66618 | 0 |
| VAL17 | PROA | 438 | 6.97547 | 0 |
| SER20 | PROA | 441 | 2.78641 | 1 |
| TYR13 | PROB | 434 | 3.36728 | 1 |
| VAL17 | PROB | 438 | 3.57028 | 1 |
| SER20 | PROB | 441 | 6.73379 | 0 |

## Final Chain Distances

| Chain | Min distance A | Within 6 A |
|---|---:|---:|
| PROA | 2.78641 | 1 |
| PROB | 3.36728 | 1 |

## Notes

- This stage is a short restrained MD pilot only, not production MD.
- No production MD was run.
- No original apo files were modified.
- The ligand `.str` file was not modified.

## Next Step

- Review the warnings in the retention CSV and log before deciding whether to repeat this pilot with adjusted restraints or shorter heating increments.
