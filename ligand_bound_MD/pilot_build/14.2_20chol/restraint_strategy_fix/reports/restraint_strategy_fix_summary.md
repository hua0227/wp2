# Restraint Strategy Fix Summary

## Scope

- ligand_id: 14.2
- resname: L002
- system: 20chol
- branch: Branch B remove_POPC94
- No production MD was run.
- No longer MD was run.
- No files outside pilot-local Branch B/fix outputs were modified.

## Why MEMB:49 Was Not Further Processed

Residual force diagnostics showed MEMB:49 was far from ligand L002, had no abnormal bonds, and dropped out of the no-restraint top force atoms. Its previous high force was dominated by absolute-coordinate CustomExternalForce restraints, so deleting or processing MEMB:49 would not address the root cause.

## PROA437 Force Source

- PROA:437 O dominant force group: torsion (CustomTorsionForce, group 4).
- PROA:437 O dominant group force: 31707.7442227 kJ/mol/nm.
- Interpretation: this is mainly internal torsion/custom torsion strain, not the old MEMB49 restraint artifact.

## Old Restraint Strategy Problem

The old strategy strongly restrained all lipid/cholesterol heavy atoms with absolute Cartesian displacement. In this periodic system, that produced restraint artifacts when coordinates/reference positions landed in different periodic images. The corrected strategy uses `periodicdistance(...)` and avoids strong lipid/cholesterol restraints.

## R1 vs R2

- R1: protein backbone medium PBC-aware restraint + ligand heavy weak PBC-aware restraint; no lipid/cholesterol restraint.
- R2: same as R1 plus extremely weak lipid/cholesterol heavy atom PBC-aware restraint.

| Branch | Final energy kJ/mol | Max force kJ/mol/nm | Top force atom | MEMB49 max | PROA437/438 max | Severe/close | Key min A | Near chains | Ligand RMSD A |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| R1 | -604351.86694 | 2438.54048764 | OH2 TIP3 TIP3:7718 index=42769 | 279.790370882 | 363.473054469 | 0/0 | 3.29666205834 | 1 | 3.24019156608 |
| R2 | -601401.197018 | 2438.89312781 | OH2 TIP3 TIP3:4391 index=32788 | 368.749201497 | 237.200192171 | 0/0 | 3.62569172881 | 2 | 2.01382186488 |

## Ligand Retention

- R1 ligand near chains: 1; R1 is less preferred.
- R2 ligand near chains: 2; R2 retained the ligand near both chains.
- R2 severe/close ligand clashes: 0 / 0.

## Very-Short Probe

- Probe branch: R2.
- Conditions: 50 K, 0.05 fs timestep, 0.2 ps total, OpenCL, PBC-aware restraints, not production MD.
- Probe executed: YES.
- Probe passed: YES.
- Final probe step/status: 4000 / ok.
- Final probe temperature: 46.5284 K.
- Final probe ligand near chains: 2.

## Recommendation

- Recommend reattempting 10 ps restrained MD: YES.
- Use the R2 PBC-aware restraint strategy, not the old absolute lipid/cholesterol restraints.
- This recommendation is for the next short restrained MD attempt only, not production MD.
- PyMOL view: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\restraint_strategy_fix\outputs\view_R1_R2_PROA437_L002.pml`
