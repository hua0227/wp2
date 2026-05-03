# Ligand-Bound MD Workflow Template

This is a workflow note only. It does not run MD automatically.

## 1. Select Ligands

Use `selected_candidates_top10.csv`, `selected_candidates_top6.csv`, or `selected_candidates_top3.csv` as the selection source. The current staging folders are under `selected_inputs/<ligand_id>/<system>`.

## 2. Prepare Ligand Parameters

For each ligand, prepare CHARMM-compatible ligand parameters before any MD run:

- Confirm the source SDF/PDBQT corresponds to the intended ligand.
- Generate or curate a CGenFF-compatible `.mol2`.
- Generate the matching CGenFF/CHARMM `.str`.
- Keep parameterization status in `parameterization_needed/ligand_parameterization_checklist.csv`.

## 3. Build Ligand-Bound Systems

For each ligand and system:

- Start from a copied working system, not the original apo OpenMM directory.
- Merge the docked ligand pose into the apo 20% or 40% CHOL membrane coordinates.
- Preserve receptor, membrane, solvent, ions, and ligand coordinate consistency.
- Check residue names and atom names against the ligand `.str` file.

## 4. Pre-MD Validation

- Confirm OpenMM can read the combined topology and coordinates.
- Check total charge, ligand atom count, and missing parameters.
- Visually inspect the merged system in PyMOL/VMD.

## 5. Short Sanity Run

After parameters and topology pass validation:

- Energy minimization.
- 5000-step short MD.
- Confirm ligand remains near TYR13/VAL17/SER20 and does not create severe clashes.

## 6. Longer Production MD

Only consider longer ligand-bound MD after the short sanity run is stable and the ligand parameters have been reviewed.
