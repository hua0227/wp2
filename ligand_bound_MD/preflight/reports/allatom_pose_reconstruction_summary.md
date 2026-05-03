# All-Atom Docked Ligand Pose Reconstruction Summary

## Why reconstruction was needed

The docked poses are not CGenFF all-atom ligand coordinate sets. They are docking outputs with heavy atoms plus limited hydrogens and may use a docking/root-branch atom order. CGenFF MOL2/STR files describe all atoms with CGenFF atom names, so an atom-count mismatch blocks the next ligand-bound topology/read-test step.

This reconstruction maps docked heavy atom coordinates back onto the CGenFF all-atom molecule and writes preliminary all-atom, CGenFF-named ligand PDB files. Heavy atom coordinates come from the docked pose. Hydrogen coordinates are preliminary local placements from the CGenFF/reference geometry and still require later minimization/relaxation before any production simulation workflow.

## Tools used

- RDKit available: YES (available)
- Open Babel available: YES (Open Babel 3.1.1 -- May 16 2020 -- 11:57:55)
- Reconstruction methods used: rdkit: 20

## Results

- Ligand-system combinations processed: 20
- All-atom named poses generated: 20
- Visual-only complex previews generated: 20

## Failures

- NONE

## Generated complex previews

- 8.1 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\8.1\20chol\complex_preview_allatom_8.1_20chol.pdb
- 8.1 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\8.1\40chol\complex_preview_allatom_8.1_40chol.pdb
- 14.2 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\14.2\20chol\complex_preview_allatom_14.2_20chol.pdb
- 14.2 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\14.2\40chol\complex_preview_allatom_14.2_40chol.pdb
- 19.1 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\19.1\20chol\complex_preview_allatom_19.1_20chol.pdb
- 19.1 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\19.1\40chol\complex_preview_allatom_19.1_40chol.pdb
- 9.2 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\9.2\20chol\complex_preview_allatom_9.2_20chol.pdb
- 9.2 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\9.2\40chol\complex_preview_allatom_9.2_40chol.pdb
- 12.2 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\12.2\20chol\complex_preview_allatom_12.2_20chol.pdb
- 12.2 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\12.2\40chol\complex_preview_allatom_12.2_40chol.pdb
- 2.3 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\2.3\20chol\complex_preview_allatom_2.3_20chol.pdb
- 2.3 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\2.3\40chol\complex_preview_allatom_2.3_40chol.pdb
- 4.3 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\4.3\20chol\complex_preview_allatom_4.3_20chol.pdb
- 4.3 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\4.3\40chol\complex_preview_allatom_4.3_40chol.pdb
- 6.2 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\6.2\20chol\complex_preview_allatom_6.2_20chol.pdb
- 6.2 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\6.2\40chol\complex_preview_allatom_6.2_40chol.pdb
- 8.3 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\8.3\20chol\complex_preview_allatom_8.3_20chol.pdb
- 8.3 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\8.3\40chol\complex_preview_allatom_8.3_40chol.pdb
- 17.2 20chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\17.2\20chol\complex_preview_allatom_17.2_20chol.pdb
- 17.2 40chol: C:\TRKB_WP2\ligand_bound_MD\preflight\complex_previews_allatom\17.2\40chol\complex_preview_allatom_17.2_40chol.pdb

## Important limitation

The complex preview PDB files are only for visualization. They are apo ATOM/HETATM records plus reconstructed ligand coordinates, not PSF/PDB/topology builds, not validated OpenMM systems, and not simulation-ready topologies.

## Next step

The next step is to attempt ligand-bound topology construction and an OpenMM read test using the CGenFF parameters and these preliminary all-atom named ligand coordinates.
