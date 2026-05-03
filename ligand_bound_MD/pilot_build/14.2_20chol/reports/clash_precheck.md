# Ligand-Apo Clash Precheck

## Scope

Pilot: ligand 14.2 / L002 / 20chol. This check reports heavy-atom contacts only and does not delete or alter any apo atoms.

## Inputs

- Apo PDB: `C:\TRKB_WP2\TRKB_20CHOL\openmm_short_md_output\short_md_final.pdb`
- Ligand all-atom PDB: `C:\TRKB_WP2\ligand_bound_MD\preflight\ligand_pose_allatom_named\14.2\20chol\L002_14.2_20chol_allatom_named.pdb`

## Summary

- Apo atoms parsed: 52357
- Apo heavy atoms checked: 18685
- Ligand atoms parsed: 48
- Ligand heavy atoms checked: 21
- Severe clashes, ligand heavy to apo heavy < 1.0 A: 9
- Close contacts, ligand heavy to apo heavy < 1.5 A: 17

Severe clashes are a subset of the close-contact count.

## Severe Clash Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 9 |
| water | 0 |
| ion | 0 |
| cholesterol | 0 |
| unknown | 0 |

## Close Contact Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 16 |
| water | 0 |
| ion | 0 |
| cholesterol | 1 |
| unknown | 0 |

## Closest Contacts

| Type | Distance A | Apo category | Ligand atom | Apo atom | Apo residue |
|---|---:|---|---|---|---|
| severe_clash | 0.271 | lipid | C 7 | C39 12226 | POP C:94 |
| severe_clash | 0.625 | lipid | C 2 | C34 12211 | POP C:94 |
| severe_clash | 0.736 | lipid | C 9 | C311 12232 | POP C:94 |
| severe_clash | 0.791 | lipid | C 6 | C37 12220 | POP C:94 |
| severe_clash | 0.873 | lipid | C 6 | C38 12223 | POP C:94 |
| severe_clash | 0.915 | lipid | C 3 | C36 12217 | POP C:94 |
| severe_clash | 0.955 | lipid | C 21 | C36 12217 | POP C:94 |
| severe_clash | 0.972 | lipid | C 5 | C37 12220 | POP C:94 |
| severe_clash | 0.991 | lipid | N 8 | C310 12229 | POP C:94 |
| close_contact | 1.004 | lipid | C 3 | C35 12214 | POP C:94 |
| close_contact | 1.200 | lipid | C 20 | C38 12223 | POP C:94 |
| close_contact | 1.309 | lipid | N 8 | C39 12226 | POP C:94 |
| close_contact | 1.332 | lipid | C 4 | C35 12214 | POP C:94 |
| close_contact | 1.398 | lipid | C 7 | C38 12223 | POP C:94 |
| close_contact | 1.439 | lipid | C 20 | C37 12220 | POP C:94 |
| close_contact | 1.466 | cholesterol | C 15 | C21 11470 | CHL C:87 |
| close_contact | 1.472 | lipid | C 12 | C311 12232 | POP C:94 |
