# Contact Occupancy / Distance Occupancy Analysis

## Scope

This analysis uses only the existing 14.2 / L002 `20chol` and `40chol` 100 ps restrained MD pilot outputs from `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\tables\combined_time_series_14.2.csv` and `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\tables\comparison_summary_statistics.csv`.

## Important Interpretation

- This is restrained MD pilot analysis, not production MD.
- 40chol still shows interface retention because key-min distance remains low overall and near chains stay at 2 for the full analyzed trajectory.
- 20chol per-residue occupancies use only frames with real per-residue distances available from the prior DCD-based reconstruction. No missing per-residue values were fabricated.

## Occupancy Summary

- TYR13 occupancy is stronger in 40chol than in 20chol:
  40chol TYR13 `<5 A` occupancy = 96.0%, versus 20chol = 0.0%.
- VAL17 occupancy drops in 40chol:
  20chol VAL17 `<5 A` occupancy = 100.0%, versus 40chol = 0.0%.
- SER20 occupancy drops in 40chol:
  20chol SER20 `<5 A` occupancy = 98.0%, versus 40chol = 0.0%.
- Key-min-distance occupancy remains high in both systems:
  20chol `<5 A` = 100.0%, 40chol `<5 A` = 96.0%.
- Near chains = 2 fraction:
  20chol = 100.0%, 40chol = 100.0%.

## RMSD Summary

- 20chol RMSD mean / median / max / 95th percentile = 0.833 / 0.813 / 1.523 / 1.198 A.
- 40chol RMSD mean / median / max / 95th percentile = 1.683 / 1.584 / 3.000 / 2.640 A.

## Interpretation

- 40chol retains the ligand at the interface, but the contact pattern shifts.
- VAL17 and SER20 contact occupancies are clearly lower in 40chol than in 20chol, with `<5 A` occupancy drops of 100.0 and 98.0 percentage points.
- TYR13 becomes more prominent in 40chol, with a `<5 A` occupancy gain of 96.0 percentage points.
- Together, these results suggest that 40chol may alter the interfacial contact pose of 14.2 rather than removing interface retention.

## Outputs

- Occupancy summary table: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\tables\contact_occupancy_summary.csv`
- Contact occupancy bar plot: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\figures\contact_occupancy_barplot.png`
- Per-residue distance box plot: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\figures\per_residue_distance_boxplot.png`
- RMSD distribution plot: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40\figures\rmsd_distribution_20_vs_40.png`
