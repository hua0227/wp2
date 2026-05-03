# 14.2 / L002 20chol vs 40chol R2 Pilot Comparison

## Purpose

Compare restrained MD pilot stability for 14.2 / L002 in 20% CHOL and 40% CHOL TRKB-TMD systems.

## Input Data

- 20chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\reports\md100ps_R2_monitor.csv`
- 20chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2.dcd`
- 20chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2_final.pdb`
- 20chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\md100ps_R2\reports\md100ps_R2_summary.md`
- 40chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv`
- 40chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md100ps_R2\outputs\TRKB_40chol_L002_14.2_100ps_R2.dcd`
- 40chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md100ps_R2\outputs\TRKB_40chol_L002_14.2_100ps_R2_final.pdb`
- 40chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md`

## Simulation Status

Both systems are 100 ps R2 PBC-aware restrained MD pilots. These are not production MD simulations.

## Data Handling

- 40chol monitor includes individual TYR13, VAL17, and SER20 distances.
- 20chol monitor includes only the overall ligand minimum distance to TYR13/VAL17/SER20.
- 20chol per-residue distances were recomputed from DCD plus final-PDB topology where trajectory frames exist. recomputed from DCD: 50 frames at 2 ps spacing, mapped to 2-100 ps
- 20chol per-residue values are available only at DCD frame times; missing 0 ps and odd-ps values were left blank and were not fabricated.

## Stability Results

- 20chol completed ps: 100; 40chol completed ps: 100.
- Temperature range 20chol: 306.937-313.268 K; 40chol: 307.013-314.195 K.
- NaN detected 20chol: False; 40chol: False.
- Severe clashes appeared 20chol: False; 40chol: False.
- Close contacts appeared 20chol: False; 40chol: False.
- Near chains remained 2 in both systems: True.

## Ligand Retention

- Ligand RMSD final/max 20chol: 1.19815/1.52269 A.
- Ligand RMSD final/max 40chol: 1.55302/3.00028 A.
- Key-residue minimum distance final/mean 20chol: 3.6668/3.74823 A.
- Key-residue minimum distance final/mean 40chol: 4.47124/3.80623 A.
- Mean TYR13/VAL17/SER20 distances 20chol: 6.87328/3.86602/4.10325 A.
- Mean TYR13/VAL17/SER20 distances 40chol: 3.80623/8.25231/8.12394 A.

## Key Comparison

- 20chol keeps L002 close to the key-residue region, with a lower overall key-min distance and lower ligand RMSD during this restrained pilot.
- 40chol still retains L002 at the interface because near chains remain 2 and TYR13 remains close, but VAL17 and SER20 contacts are weaker than in 20chol based on available per-residue distances.
- 40chol VAL17/SER20 clearly farther than 20chol: YES. 40chol VAL17/SER20 mean distances are 8.252/8.124 A versus 20chol 3.866/4.103 A.
- This pattern suggests 40chol may alter the ligand contact pose within the interface rather than fully ejecting the ligand.

## Limitations

- 100 ps restrained MD is not production MD.
- Both systems used restraints, so these data are pilot stability checks rather than equilibrium conclusions.
- CGenFF parameter penalty for L002 must remain part of interpretation.
- Only one ligand was analyzed.
- More ligands, longer simulations, reduced restraints, and replicate trajectories are needed before drawing production-level conclusions.

## Next Steps

1. Treat 14.2 as the first focused candidate for cross-cholesterol comparison.
2. Extend this workflow to another moderate-penalty candidate such as 19.1 or 2.3.
3. Perform a more systematic 14.2 contact occupancy / distance occupancy analysis.
4. Do not describe these restrained pilot results as production MD conclusions.

## Outputs

- Summary statistics table: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\tables\comparison_summary_statistics.csv`
- Combined time series table: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\tables\combined_time_series_14.2.csv`
- Figures directory: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\figures`
- PyMOL visualization script: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\pymol\view_14.2_20_vs_40_final.pml`
