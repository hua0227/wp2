# Clash Resolution and Restrained Minimization Summary

## Scope

This stage handled only ligand 14.2 / L002 / 20chol. It used the pilot ligand-bound PSF/PDB from the previous topology build and performed clash classification, pilot-local cleaned-system staging, and restrained minimization. No production MD, short MD, long MD, or trajectory generation was run.

## Original Bound-System Clashes

- Severe clashes, ligand heavy to non-ligand heavy < 1.0 A: 9
- Close contacts, ligand heavy to non-ligand heavy < 1.5 A: 17
- Severe clash categories: lipid 9; protein 0; cholesterol 0; water 0; ion 0; other/unknown 0
- Close contact categories: lipid 16; cholesterol 1; protein 0; water 0; ion 0; other/unknown 0

Detailed classification:
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\reports\bound_clash_classification.csv`
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\reports\bound_clash_classification.md`

## Water/Ion Cleaning

- Severe water/ion clash residues found: 0
- Water/ion residues deleted: 0
- Protein/lipid/cholesterol/ligand deleted: 0

No water or ion deletion was allowed or needed. Because all severe clashes were lipid contacts, no atoms were removed. The cleaned PSF/PDB are pilot-local copies of the bound PSF/PDB.

Cleaned outputs:
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_bound_cleaned.psf`
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_bound_cleaned.pdb`

Cleaned-system VMD log:
`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\logs\cleaned_system_build.log`

## Restrained Minimization

- OpenMM Python: `C:\Users\14566\miniconda3\envs\trkb_openmm\python.exe`
- OpenMM platform used: OpenCL
- CUDA used: NO
- Strong restraints: protein backbone and lipid/cholesterol heavy atoms
- Weak restraints: ligand heavy atoms
- Water/ions restrained: NO
- Minimization max iterations: 1000
- Restrained minimization success: YES
- Initial potential energy: `19539656635280.184 kJ/mol`
- Final potential energy: `396752606.00629056 kJ/mol`
- Initial/final energy finite: YES
- Energy decreased: YES
- NaN energy observed: NO

Minimized output:
`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_minimized.pdb`

Minimization log:
`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\logs\restrained_minimization.log`

## Post-Minimization Contact Check

- Post-minimized severe clashes < 1.0 A: 0
- Post-minimized close contacts < 1.5 A: 0

## Ligand Retention

- Ligand heavy atom RMSD before vs after minimization: 1.985 A
- Ligand within 6 A of at least one TYR13/VAL17/SER20 key residue: YES
- Ligand within 6 A of two protein chains: YES

Key-residue distances after minimization:

| Key residue | Segid | Resid | Min distance A | Within 6 A |
|---|---|---:|---:|---:|
| TYR13 | PROA | 434 | 7.690 | 0 |
| VAL17 | PROA | 438 | 7.055 | 0 |
| SER20 | PROA | 441 | 2.573 | 1 |
| TYR13 | PROB | 434 | 3.365 | 1 |
| VAL17 | PROB | 438 | 3.570 | 1 |
| SER20 | PROB | 441 | 6.576 | 0 |
| nearest protein chain | PROA |  | 2.573 | 1 |
| nearest protein chain | PROB |  | 3.365 | 1 |

Detailed retention reports:
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\reports\minimization_retention_check.csv`
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\reports\minimization_retention_check.md`

## Short MD Pilot Readiness

Recommendation: YES, cautiously proceed to a short restrained MD pilot only if the next stage remains explicitly non-production and monitored.

Rationale: severe/close ligand contacts were removed by restrained minimization, the energy is finite and decreased substantially, and the ligand remains near key residues and both protein chains. The final energy remains large in absolute terms, so this should not be treated as production-ready.

## Explicit Non-Actions

- No MD was run in this stage.
- No production MD was run.
- No trajectory was generated.
- Original apo files were not modified.
- The original `step5_assembly.psf` was not modified.
- The ligand `.str` file was not modified.
- No protein, lipid, cholesterol, ligand, water, or ion atoms were deleted from the original files.
- No Top10 batch processing was performed.
