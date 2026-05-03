# Detailed Ligand-Lipid Clash Diagnosis

## Scope

- Pilot: ligand_id 14.2 / resname L002 / system 20chol.
- Bound PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_L002_14.2_bound.psf`
- Bound PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_L002_14.2_bound.pdb`
- Previous minimized PDB checked: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\clash_resolution\outputs\TRKB_20chol_L002_14.2_minimized.pdb`
- Contacts include all ligand atoms, including hydrogens, against all non-ligand atoms.
- Severe cutoff: < 1.0 A.
- Close-contact cutoff: < 1.5 A, with severe contacts included in the total close count.

## Summary

| Source | Severe <1.0 A | Close <1.5 A | POPC94 severe | POPC94 close | Ligand atoms severe with POPC94 | H9/H11 severe with POPC94 | MEMB:94 POPC atoms | POPC94 dominant severe source |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| bound_initial | 27 | 81 | 22 | 63 | C12, C19, C2, C3, C5, C6, C7, C8, H14, H19, H26, H27, H7, H9, N | 1 | 134 | 1 |
| previous_minimized | 0 | 1 | 0 | 0 | none | 0 | 134 | 0 |

## Atom Pair Classes

| Source | Cutoff | heavy-heavy | heavy-H | H-heavy | H-H |
|---|---|---:|---:|---:|---:|
| bound_initial | severe_lt_1.0A | 9 | 7 | 5 | 6 |
| bound_initial | close_1.0_to_1.5A | 8 | 17 | 9 | 20 |
| previous_minimized | severe_lt_1.0A | 0 | 0 | 0 | 0 |
| previous_minimized | close_1.0_to_1.5A | 0 | 0 | 0 | 1 |

## Top Severe-Clashing Residues

### bound_initial

| Rank | Residue | Severe contact count | POPC94 |
|---:|---|---:|---:|
| 1 | POPC MEMB:94 | 22 | 1 |
| 2 | CHL1 MEMB:87 | 5 | 0 |

### previous_minimized

| Rank | Residue | Severe contact count | POPC94 |
|---:|---|---:|---:|
|  | none | 0 | 0 |

## LIG H9/H11 Contacts

### bound_initial

| Cutoff | Distance A | Ligand | Other atom | Other residue | Pair class | POPC94 |
|---|---:|---|---|---|---|---:|
| severe_lt_1.0A | 0.8605 | H9 | C312 12233 | POPC MEMB:94 | H-heavy | 1 |
| close_1.0_to_1.5A | 1.1326 | H9 | H12Y 12235 | POPC MEMB:94 | H-H | 1 |
| close_1.0_to_1.5A | 1.1338 | H9 | C311 12230 | POPC MEMB:94 | H-heavy | 1 |

### previous_minimized

| Cutoff | Distance A | Ligand | Other atom | Other residue | Pair class | POPC94 |
|---|---:|---|---|---|---|---:|
| none |  |  |  |  |  | 0 |

## Closest Contacts

| Source | Cutoff | Distance A | Pair class | Ligand atom | Other atom | Other residue | Category | POPC94 |
|---|---|---:|---|---|---|---|---|---:|
| bound_initial | severe_lt_1.0A | 0.2710 | heavy-heavy | C7 LIG:1 | C39 12224 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.4632 | H-H | H24 LIG:1 | H21B 11470 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | severe_lt_1.0A | 0.5338 | H-H | H19 LIG:1 | H12Y 12235 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.5808 | H-heavy | H25 LIG:1 | C21 11468 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | severe_lt_1.0A | 0.6101 | heavy-H | C8 LIG:1 | H11Y 12232 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.6255 | heavy-heavy | C2 LIG:1 | C34 12209 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.6396 | H-heavy | H7 LIG:1 | C38 12221 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.6678 | heavy-H | C19 LIG:1 | H6X 12216 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.6832 | heavy-H | C2 LIG:1 | H4Y 12211 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.7363 | heavy-heavy | C8 LIG:1 | C311 12230 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.7541 | heavy-H | C12 LIG:1 | H9Y 12226 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.7912 | heavy-heavy | C6 LIG:1 | C37 12218 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.8022 | heavy-H | C5 LIG:1 | H7Y 12220 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.8123 | H-heavy | H26 LIG:1 | C36 12215 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.8447 | heavy-H | C15 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | severe_lt_1.0A | 0.8483 | H-heavy | H26 LIG:1 | C35 12212 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.8605 | H-heavy | H9 LIG:1 | C312 12233 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.8673 | H-H | H16 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | severe_lt_1.0A | 0.8735 | heavy-heavy | C6 LIG:1 | C38 12221 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9148 | heavy-heavy | C3 LIG:1 | C36 12215 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9309 | H-H | H14 LIG:1 | H11X 12231 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9551 | heavy-heavy | C19 LIG:1 | C36 12215 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9622 | H-H | H27 LIG:1 | H5Y 12214 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9719 | heavy-heavy | C5 LIG:1 | C37 12218 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9830 | heavy-H | C5 LIG:1 | H7X 12219 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9910 | heavy-heavy | N LIG:1 | C310 12227 | POPC MEMB:94 | lipid | 1 |
| bound_initial | severe_lt_1.0A | 0.9918 | H-H | H18 LIG:1 | H22A 11473 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.0041 | heavy-heavy | C3 LIG:1 | C35 12212 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.0117 | H-H | H2 LIG:1 | HB2 295 | ALA PROA:440 | protein | 0 |
| bound_initial | close_1.0_to_1.5A | 1.0312 | heavy-H | C7 LIG:1 | H9Y 12226 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.0611 | H-heavy | H27 LIG:1 | C35 12212 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.0620 | H-heavy | H7 LIG:1 | C39 12224 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.0756 | H-H | H23 LIG:1 | H21C 11471 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.0792 | heavy-H | N LIG:1 | H10Y 12229 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.0937 | H-H | H5 LIG:1 | H5Y 12214 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1039 | H-H | H21 LIG:1 | H9Y 12226 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1326 | H-H | H9 LIG:1 | H12Y 12235 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1338 | H-heavy | H9 LIG:1 | C311 12230 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1342 | heavy-H | C9 LIG:1 | H11Y 12232 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1378 | heavy-H | C3 LIG:1 | H6X 12216 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1438 | heavy-H | C4 LIG:1 | H5Y 12214 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1440 | H-H | H6 LIG:1 | H7Y 12220 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1496 | H-heavy | H6 LIG:1 | C37 12218 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1499 | H-heavy | H25 LIG:1 | C20 11466 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.1656 | heavy-H | C14 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.1983 | H-H | H7 LIG:1 | H8Y 12223 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.1985 | H-heavy | H24 LIG:1 | C21 11468 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.2004 | heavy-heavy | C18 LIG:1 | C38 12221 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2094 | H-H | H8 LIG:1 | H9X 12225 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2152 | H-heavy | H6 LIG:1 | C38 12221 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2239 | H-H | H25 LIG:1 | H21C 11471 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.2324 | heavy-H | C11 LIG:1 | H11X 12231 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2387 | H-heavy | H27 LIG:1 | C34 12209 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2455 | heavy-H | C11 LIG:1 | H11Y 12232 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.2565 | H-H | H13 LIG:1 | H26A 2482 | CHL1 MEMB:14 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.2609 | H-H | H20 LIG:1 | H22A 11473 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.2977 | heavy-H | C18 LIG:1 | H8Y 12223 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3069 | heavy-H | C19 LIG:1 | H6Y 12217 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3094 | heavy-heavy | N LIG:1 | C39 12224 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3211 | H-H | H25 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.3322 | heavy-heavy | C4 LIG:1 | C35 12212 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3389 | H-H | H26 LIG:1 | H6Y 12217 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3530 | heavy-H | C17 LIG:1 | H9X 12225 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3675 | heavy-H | C7 LIG:1 | H9X 12225 | POPC MEMB:94 | lipid | 1 |
| previous_minimized | close_1.0_to_1.5A | 1.3742 | H-H | H2 LIG:1 | HG13 231 | VAL PROA:436 | protein | 0 |
| bound_initial | close_1.0_to_1.5A | 1.3905 | H-H | H6 LIG:1 | H8X 12222 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3950 | H-H | H26 LIG:1 | H5X 12213 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3967 | heavy-H | N LIG:1 | H9X 12225 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.3983 | heavy-heavy | C7 LIG:1 | C38 12221 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4079 | heavy-H | C2 LIG:1 | H4X 12210 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4272 | heavy-H | C6 LIG:1 | H7X 12219 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4391 | heavy-heavy | C18 LIG:1 | C37 12218 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4428 | H-H | H7 LIG:1 | H8X 12222 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4464 | heavy-H | C6 LIG:1 | H7Y 12220 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4509 | H-H | H17 LIG:1 | H9X 12225 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4645 | H-H | H23 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.4657 | heavy-heavy | C14 LIG:1 | C21 11468 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.4722 | heavy-heavy | C11 LIG:1 | C311 12230 | POPC MEMB:94 | lipid | 1 |
| bound_initial | close_1.0_to_1.5A | 1.4724 | heavy-H | C16 LIG:1 | H21A 11469 | CHL1 MEMB:87 | cholesterol | 0 |
| bound_initial | close_1.0_to_1.5A | 1.4731 | H-H | H22 LIG:1 | H9Y 12226 | POPC MEMB:94 | lipid | 1 |

## Deletion Decision Input

Branch B is permitted to remove only the whole POPC molecule MEMB:94 in a pilot-local cleaned system. The deletion decision uses the bound Branch B starting coordinates plus the prior failed-MD force diagnostic that implicated LIG H9/H11 and MEMB:94 POPC. The previous minimized coordinates are also listed because they show whether simple geometric contacts remained after the earlier minimization.
- MEMB:94 POPC dominant severe source in bound_initial: YES
- MEMB:94 POPC dominant severe source in previous_minimized: NO
