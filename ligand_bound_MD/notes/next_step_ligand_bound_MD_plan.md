# Ligand-Bound MD Preparation Plan

## Current Stage Goal

Prepare selected ligand-bound MD input folders for the top-ranked TRKB-TMD docking candidates without running ligand-bound MD and without modifying the original apo OpenMM systems.

## Completed Upstream Work

- Apo TRKB-TMD short OpenMM MD has been completed for 20% CHOL and 40% CHOL systems.
- AutoDock Vina docking has been completed for candidate ligands and reference ligands.
- All-mode contact analysis has been completed.
- Integrated scoring has been generated from docking affinity, TYR13/VAL17/SER20 contacts, two-chain proximity, pose mode reliability, 20%/40% CHOL robustness, and reference anchors.
- PyMOL packages are available for candidate and reference pose inspection.

## Scope of This Stage

This stage only prepares ligand-bound MD input folders. It does not run ligand-bound MD, does not regenerate receptors, and does not modify the original apo system directories.

## Ligand Parameterization Requirement

Each ligand needs CGenFF/CHARMM-compatible parameters before OpenMM MD. At minimum, each selected ligand should have:

- A ligand `.mol2` file with checked atom names and chemistry.
- A CGenFF/CHARMM `.str` stream file containing ligand topology and parameters.

Do not run OpenMM ligand-bound MD directly without ligand force-field parameters. Missing parameter files should remain marked as `pending` in the checklist.

## Later Formal Ligand-Bound MD Workflow

1. Ligand parameterization.
2. Merge the docked ligand pose back into the 20%/40% CHOL apo membrane system.
3. Check topology, residue name, atom naming, protonation, charge, and coordinate consistency.
4. Run an OpenMM read/load test.
5. Run energy minimization.
6. Run a 5000-step short MD sanity check.
7. Analyze whether the ligand remains near TYR13/VAL17/SER20.
