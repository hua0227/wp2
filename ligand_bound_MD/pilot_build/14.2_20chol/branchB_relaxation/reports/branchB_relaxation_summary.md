# Branch B Relaxation Summary

## Scope

- ligand_id: 14.2
- resname: L002
- system: 20chol
- branch: Branch B remove_POPC94
- This stage used only Branch B pilot-local files.
- No production MD was run.
- No longer MD was run.
- Original apo files, receptor files, docking, and ligand `.str` were not modified.

## Why Round-Trip And Staged Relaxation Were Needed

The targeted cleanup removed POPC MEMB:94 and eliminated severe/close ligand clashes, but the ultra-conservative 1 ps 50 K check stopped after 10 steps because temperature exceeded 200 K. The next diagnostic layer therefore checked whether the Branch B minimized PDB is handed to OpenMM exactly, then inspected the force distribution and relaxed the system in staged restraint schedules without starting dynamics.

## Coordinate Handoff

- Coordinate handoff passed: YES.
- PSF atom count: 52271.
- Cleaned PDB atom count: 52271.
- Minimized PDB atom count: 52271.
- Roundtrip heavy-atom RMSD: 0 A.
- Roundtrip ligand heavy-atom RMSD: 0 A.
- Periodic box max delta: 0 A.

## Initial Force Diagnostics

- Initial Branch B potential energy: -332119.090573 kJ/mol.
- Force NaN/nonfinite present: NO.
- Highest-force atom: C36 POPC MEMB:49 index 7084.
- Highest-force magnitude: 11003.2118392 kJ/mol/nm.
- Highest forces concentrated in ligand H atoms: NO.

## Highest Force Atoms

| Rank | Atom | Residue | Category | Force kJ/mol/nm |
|---:|---|---|---|---:|
| 1 | C36 index 7084 | POPC MEMB:49 | lipid | 11003.212 |
| 2 | C38 index 5750 | POPC MEMB:39 | lipid | 10017.316 |
| 3 | C22 index 17908 | POPC MEMB:138 | lipid | 8975.6267 |
| 4 | C35 index 7081 | POPC MEMB:49 | lipid | 8723.1313 |
| 5 | C22 index 13084 | POPC MEMB:102 | lipid | 8510.5433 |
| 6 | C36 index 12616 | POPC MEMB:98 | lipid | 8390.4539 |
| 7 | C26 index 8109 | POPC MEMB:57 | lipid | 7955.7669 |
| 8 | C27 index 15118 | POPC MEMB:117 | lipid | 7702.2679 |
| 9 | C24 index 15109 | POPC MEMB:117 | lipid | 7655.5484 |
| 10 | C37 index 7087 | POPC MEMB:49 | lipid | 7615.3443 |

## Staged Minimization

| Stage | Physical energy kJ/mol | Physical max force kJ/mol/nm | Max force atom | Severe/close | Key min A | Near chains | Ligand RMSD A |
|---:|---:|---:|---|---:|---:|---:|---:|
| 1 | -218465.80151 | 8497.4351 | C25 POPC MEMB:36 index=5292 | 0/0 | 3.42301 | 2 | 1.57033 |
| 2 | -143535.160885 | 9832.8889 | C25 POPC MEMB:91 index=11762 | 0/0 | 2.93018 | 2 | 1.31672 |
| 3 | -143534.455807 | 9832.8593 | C25 POPC MEMB:91 index=11762 | 0/0 | 2.93018 | 2 | 1.31672 |
| 4 | -143534.610104 | 9832.8554 | C25 POPC MEMB:91 index=11762 | 0/0 | 2.93018 | 2 | 1.31672 |

## Ligand Retention And Clashes

- Stage 4 severe/close clashes: 0 / 0.
- Stage 4 ligand min distance to TYR13/VAL17/SER20: 2.93018 A.
- Stage 4 ligand near chains count: 2.
- Stage 4 ligand RMSD vs Branch B starting ligand: 1.31672 A.
- Stage 4 max force reduction vs initial diagnostic: 10.6%.

## Very-Short Probe

- Probe executed: NO.
- Probe passed: NO.
- Probe was skipped because:
  - Max force reduction was only 10.6%, below the 50% significance threshold used for this gate

## Recommendation

- Recommend next 10 ps restrained MD attempt: NO.
- Reason: staged minimization preserved ligand retention and kept clashes at zero, but Stage 4 did not significantly reduce the maximum force, so the requested gate for the 0.2 ps probe was not met.
- Continue diagnosing the residual high-force POPC region locally before another short restrained MD attempt.
- No production MD was run in this stage.
