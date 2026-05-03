# 14.2 / L002 / 40chol Non-Production Pilot Summary

## Purpose

This stage tests whether the ligand-bound 40% CHOL system can be built, read by OpenMM, diagnosed for clashes, and minimized with R2 PBC-aware restraints before any longer restrained MD.

## Relationship To 20chol Pilot

The 20chol pilot is used only as a workflow reference. No 20chol coordinates, cleaned systems, or extended 20chol MD products were copied into this 40chol build.

## Build And Read Test

- Ligand-bound PSF generated: True
- Ligand-bound PDB generated: True
- psfgen log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\logs\psfgen_build_40chol.log`
- OpenMM read test log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\logs\openmm_readtest_40chol.log`
- OpenMM initial energy: 387732754963599.44 kJ/mol

## Clash Diagnosis And Cleanup

- Severe ligand/non-ligand heavy-atom clashes < 1.0 A: 8
- Close ligand/non-ligand heavy-atom contacts < 1.5 A: 15
- Targeted cleanup: deleted dominant severe-clashing lipid molecule POPC MEMB:113; deleted lipid: POPC MEMB:113
- Cleaned PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound_cleaned.psf`
- Cleaned PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound_cleaned.pdb`

## R2 PBC-Aware Minimization

- R2 minimization success: True
- Minimized PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_R2_minimized.pdb`
- Final energy: -615048.277368 kJ/mol
- Max force: 2343.79254074 kJ/mol/nm
- Top force atom: OH2 TIP3 TIP3:7831 index=42408
- Ligand min distance to TYR13: 3.6061 A
- Ligand min distance to VAL17: 3.69268 A
- Ligand min distance to SER20: 5.12341 A
- Ligand near chains count: 2
- Severe/close after minimization: 0/0

## Very-Short Probe

- Very-short 0.2 ps 50 K probe executed: True
- Very-short probe passed: True
- Very-short probe status: last step 4000 status ok
- Very-short probe log: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\logs\very_short_probe_40chol.log`

## Recommendation

- Recommend entering 10 ps 40chol restrained MD: True
- This is not production MD.
