# Bound-System Clash Classification

## Scope

Pilot: ligand 14.2 / L002 / 20chol bound PSF/PDB. This report classifies ligand heavy atom contacts against non-ligand heavy atoms.

## Summary

- Bound PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_L002_14.2_bound.pdb`
- Ligand heavy atoms checked: 21
- Non-ligand heavy atoms checked: 18685
- Severe clashes < 1.0 A: 9
- Close contacts < 1.5 A: 17

Severe clashes are included in the close-contact count.

## Severe Clash Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 9 |
| cholesterol | 0 |
| water | 0 |
| ion | 0 |
| other/unknown | 0 |

## Close Contact Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 16 |
| cholesterol | 1 |
| water | 0 |
| ion | 0 |
| other/unknown | 0 |

## Water/Ion Deletion Decision

No severe-clashing water or ion residues were found. No atoms will be deleted; cleaned outputs will be pilot-local copies of the bound PSF/PDB.

## Closest Contacts

| Type | Distance A | Category | Ligand atom | Other atom | Other residue |
|---|---:|---|---|---|---|
| severe_clash | 0.271 | lipid | C7 LIG:1 | C39 12224 | POPC MEMB:94 |
| severe_clash | 0.625 | lipid | C2 LIG:1 | C34 12209 | POPC MEMB:94 |
| severe_clash | 0.736 | lipid | C8 LIG:1 | C311 12230 | POPC MEMB:94 |
| severe_clash | 0.791 | lipid | C6 LIG:1 | C37 12218 | POPC MEMB:94 |
| severe_clash | 0.873 | lipid | C6 LIG:1 | C38 12221 | POPC MEMB:94 |
| severe_clash | 0.915 | lipid | C3 LIG:1 | C36 12215 | POPC MEMB:94 |
| severe_clash | 0.955 | lipid | C19 LIG:1 | C36 12215 | POPC MEMB:94 |
| severe_clash | 0.972 | lipid | C5 LIG:1 | C37 12218 | POPC MEMB:94 |
| severe_clash | 0.991 | lipid | N LIG:1 | C310 12227 | POPC MEMB:94 |
| close_contact | 1.004 | lipid | C3 LIG:1 | C35 12212 | POPC MEMB:94 |
| close_contact | 1.200 | lipid | C18 LIG:1 | C38 12221 | POPC MEMB:94 |
| close_contact | 1.309 | lipid | N LIG:1 | C39 12224 | POPC MEMB:94 |
| close_contact | 1.332 | lipid | C4 LIG:1 | C35 12212 | POPC MEMB:94 |
| close_contact | 1.398 | lipid | C7 LIG:1 | C38 12221 | POPC MEMB:94 |
| close_contact | 1.439 | lipid | C18 LIG:1 | C37 12218 | POPC MEMB:94 |
| close_contact | 1.466 | cholesterol | C14 LIG:1 | C21 11468 | CHL1 MEMB:87 |
| close_contact | 1.472 | lipid | C11 LIG:1 | C311 12230 | POPC MEMB:94 |
