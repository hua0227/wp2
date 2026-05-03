# 40chol R2 Minimization Summary

## Scope

14.2 / L002 / 40chol non-production pilot. OpenCL was used. No MD was run in this minimization step.

## Results

- Cleaned PSF: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound_cleaned.psf`
- Cleaned PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_bound_cleaned.pdb`
- Minimized PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\outputs\TRKB_40chol_L002_14.2_R2_minimized.pdb`
- Initial restrained energy: 3.87725518402e+14 kJ/mol
- Final restrained energy: -615048.277368 kJ/mol
- Energy finite: True
- Force NaN: False
- Coordinate NaN: False
- Max force: 2343.79254074 kJ/mol/nm
- Top force atom: OH2 TIP3 TIP3:7831 index=42408
- Ligand min distance to TYR13: 3.6061 A
- Ligand min distance to VAL17: 3.69268 A
- Ligand min distance to SER20: 5.12341 A
- Ligand near chains count: 2
- Severe clashes < 1.0 A after minimization: 0
- Close contacts < 1.5 A after minimization: 0
- Very-short probe gate passed: True

R2 restraints used periodicdistance anchors: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.
