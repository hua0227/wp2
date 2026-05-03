# 40chol Ligand Clash Diagnosis

## Scope

Pilot: 14.2 / L002 / 40chol ligand-bound system. Ligand heavy atoms were checked against non-ligand heavy atoms.

## Summary

- Bound PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound.pdb`
- Ligand heavy atoms checked: 21
- Non-ligand heavy atoms checked: 18163
- Severe clashes < 1.0 A: 8
- Close contacts < 1.5 A: 15
- Severe clashes are included in the close-contact count.

## Severe Clash Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 8 |
| cholesterol | 0 |
| water | 0 |
| ion | 0 |
| unknown | 0 |

## Close Contact Categories

| Category | Count |
|---|---:|
| protein | 0 |
| lipid | 15 |
| cholesterol | 0 |
| water | 0 |
| ion | 0 |
| unknown | 0 |

## Targeted Cleanup Decision

- Decision: deleted dominant severe-clashing lipid molecule POPC MEMB:113
- Selected whole molecules/residues:
  - lipid: POPC MEMB:113

## Severe Clash Molecule Counts

| Category | Molecule | Severe count |
|---|---|---:|
| lipid | POPC MEMB:113 | 8 |

## Closest Contacts

| Type | Distance A | Category | Ligand atom | Other atom | Other residue |
|---|---:|---|---|---|---|
| severe_clash | 0.585 | lipid | C7 LIG:1 | C27 12677 | POPC MEMB:113 |
| severe_clash | 0.721 | lipid | C18 LIG:1 | C210 12685 | POPC MEMB:113 |
| severe_clash | 0.806 | lipid | C3 LIG:1 | C211 12687 | POPC MEMB:113 |
| severe_clash | 0.872 | lipid | C19 LIG:1 | C210 12685 | POPC MEMB:113 |
| severe_clash | 0.893 | lipid | C19 LIG:1 | C211 12687 | POPC MEMB:113 |
| severe_clash | 0.926 | lipid | O LIG:1 | C25 12671 | POPC MEMB:113 |
| severe_clash | 0.981 | lipid | C17 LIG:1 | C26 12674 | POPC MEMB:113 |
| severe_clash | 0.990 | lipid | C18 LIG:1 | C29 12683 | POPC MEMB:113 |
| close_contact | 1.158 | lipid | C17 LIG:1 | C27 12677 | POPC MEMB:113 |
| close_contact | 1.285 | lipid | C2 LIG:1 | C213 12693 | POPC MEMB:113 |
| close_contact | 1.299 | lipid | C7 LIG:1 | C28 12680 | POPC MEMB:113 |
| close_contact | 1.309 | lipid | O LIG:1 | C26 12674 | POPC MEMB:113 |
| close_contact | 1.325 | lipid | C6 LIG:1 | C28 12680 | POPC MEMB:113 |
| close_contact | 1.343 | lipid | C3 LIG:1 | C212 12690 | POPC MEMB:113 |
| close_contact | 1.431 | lipid | C2 LIG:1 | C212 12690 | POPC MEMB:113 |
