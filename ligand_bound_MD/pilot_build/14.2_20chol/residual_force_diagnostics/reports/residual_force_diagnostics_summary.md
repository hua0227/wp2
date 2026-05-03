# Residual Force Diagnostics Summary

## Scope

- ligand_id: 14.2
- resname: L002
- system: 20chol
- branch: Branch B remove_POPC94
- This stage is diagnostic only; no structure cleanup or deletion was performed.
- No MD of any kind was run.
- No original apo files, receptor files, docking outputs, or ligand `.str` files were modified.

## MEMB:49 High Force Status

- MEMB:49 high force with Stage 4 restraints: YES, max 43018.2617188 kJ/mol/nm.
- Dominant with-restraint force group: CustomExternalForce/restraints (CustomExternalForce, group 9).
- Dominant with-restraint MEMB:49 atom: N POPC MEMB:49 index=6984.
- Dominant non-restraint MEMB:49 force group: HarmonicAngleForce (HarmonicAngleForce), max 4797.37415152 kJ/mol/nm.
- No-restraint MEMB:49 max force: 2751.0359106 kJ/mol/nm.
- No-restraint MEMB:49 in top 50: NO.

## Force Source Interpretation

- Restraints: YES, the MEMB:49 high force is dominated by the CustomExternalForce restraint group when Stage 4-style restraints are present.
- Bond/angle/internal geometry: no abnormal MEMB:49 bonds were found; the largest non-restraint MEMB:49 group is HarmonicAngleForce.
- Nonbonded contacts: nonbonded contribution exists but is not the dominant with-restraint source; no-restraint top force atom is O VAL PROA:437 index 252.
- Ligand-lipid interaction: not supported by the distance/contact diagnostics.

## MEMB:49 Geometry And Contacts

- Abnormal MEMB:49 bond count: 0.
- MEMB:49 to ligand L002 minimum distance: 25.9557 A.
- MEMB:49 still clashes with ligand: NO.
- MEMB:49 has NaN coordinates: NO.
- MEMB:49 near box boundary: NO; min boundary margin 5.74105 A.
- Original apo MEMB:49 heavy-atom RMSD: nan A; mapped atoms 0.

## Top With-Restraint Force Atoms

| Rank | Atom | Residue | Category | Force kJ/mol/nm | MEMB49 |
|---:|---|---|---|---:|---:|
| 1 | C214 index 12323 | POPC MEMB:96 | lipid | 91855.009 | 0 |
| 2 | C212 index 12317 | POPC MEMB:96 | lipid | 91160.378 | 0 |
| 3 | C19 index 11603 | CHL1 MEMB:89 | cholesterol | 90916.238 | 0 |
| 4 | C13 index 17612 | POPC MEMB:136 | lipid | 90777.178 | 0 |
| 5 | C14 index 13328 | POPC MEMB:104 | lipid | 90729.676 | 0 |
| 6 | C211 index 12314 | POPC MEMB:96 | lipid | 90654.373 | 0 |
| 7 | C314 index 16258 | POPC MEMB:125 | lipid | 90631.812 | 0 |
| 8 | C14 index 12658 | POPC MEMB:99 | lipid | 90623.396 | 0 |
| 9 | O13 index 15886 | POPC MEMB:123 | lipid | 90595.326 | 0 |
| 10 | O13 index 16154 | POPC MEMB:125 | lipid | 90575.387 | 0 |

## Top No-Restraint Force Atoms

| Rank | Atom | Residue | Category | Force kJ/mol/nm | MEMB49 |
|---:|---|---|---|---:|---:|
| 1 | O index 252 | VAL PROA:437 | protein | 35006.486 | 0 |
| 2 | N index 253 | VAL PROA:438 | protein | 25020.037 | 0 |
| 3 | C25 index 11762 | POPC MEMB:91 | lipid | 9864.8587 | 0 |
| 4 | C212 index 4909 | POPC MEMB:33 | lipid | 9547.0683 | 0 |
| 5 | C35 index 6143 | POPC MEMB:42 | lipid | 9346.0167 | 0 |
| 6 | C216 index 19163 | POPC MEMB:147 | lipid | 9307.3043 | 0 |
| 7 | C211 index 4906 | POPC MEMB:33 | lipid | 8992.8406 | 0 |
| 8 | C313 index 8177 | POPC MEMB:57 | lipid | 8139.4421 | 0 |
| 9 | C index 251 | VAL PROA:437 | protein | 7948.7195 | 0 |
| 10 | C26 index 6501 | POPC MEMB:45 | lipid | 7595.2636 | 0 |

## Recommendations

- A. Continue minimization: NO - MEMB:49 has no abnormal bonds and is not the no-restraint top force source; blind minimization is unlikely to target the root cause.
- B. Adjust restraints: YES - the with-restraint MEMB:49 force is dominated by CustomExternalForce/restraints, consistent with absolute position restraints interacting with periodic/wrapped coordinates.
- C. Locally delete/process MEMB:49: NO - MEMB:49 is far from ligand, has no abnormal bonds, and drops out of the no-restraint top 50.
- D. Abandon current 14.2 pilot: NO - ligand pose remains clean; diagnose restraint/reference imaging and the no-restraint top protein force before abandoning.

Recommended next step: fix the restraint/reference-coordinate handling before any new short restrained MD attempt. In particular, avoid absolute lipid position restraints on coordinates that may be imaged across the periodic box, or regenerate restraint reference coordinates in the same periodic image as the OpenMM Context. Also inspect the no-restraint top protein force around PROA:437 before restarting a stability probe.

PyMOL view script: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\residual_force_diagnostics\outputs\view_L002_MEMB49.pml`

No MD was run in this diagnostic stage.
