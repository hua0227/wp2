# 2.3 / L006 20chol vs 40chol R2 Pilot Comparison

## 分析目的

比较 2.3 / L006 在 20% CHOL 与 40% CHOL TRKB-TMD 体系中的 restrained MD pilot 稳定性。

## 输入数据

- 20chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md100ps_R2\reports\md100ps_20chol_R2_monitor.csv`
- 20chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md100ps_R2\outputs\TRKB_20chol_L006_2.3_100ps_R2.dcd`
- 20chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md100ps_R2\outputs\TRKB_20chol_L006_2.3_100ps_R2_final.pdb`
- 20chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_20chol\md100ps_R2\reports\md100ps_20chol_R2_summary.md`
- 40chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv`
- 40chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md100ps_R2\outputs\TRKB_40chol_L006_2.3_100ps_R2.dcd`
- 40chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md100ps_R2\outputs\TRKB_40chol_L006_2.3_100ps_R2_final.pdb`
- 40chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\2.3_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md`

## 模拟状态说明

两个体系都属于 100 ps R2 PBC-aware restrained MD pilot。它们不是 production MD 模拟，本报告也不是 production MD 结论。

## 稳定性结果

- 20chol completed ps: 100.000; 40chol completed ps: 100.000
- 20chol temperature range: 307.158-314.423 K
- 40chol temperature range: 307.393-313.232 K
- NaN detected, 20chol/40chol: NO / NO
- severe clash ever, 20chol/40chol: NO / NO
- close clash ever, 20chol/40chol: NO / NO
- near chains remained 2 throughout, 20chol/40chol: YES / YES

## Ligand 保留情况

- 20chol ligand RMSD final/max: 1.291/2.206 A
- 40chol ligand RMSD final/max: 1.596/2.276 A
- 20chol TYR13 / VAL17 / SER20 final: 3.631 / 3.978 / 4.694 A
- 40chol TYR13 / VAL17 / SER20 final: 3.608 / 4.196 / 5.640 A
- 20chol TYR13 / VAL17 / SER20 mean: 3.724 / 4.088 / 5.802 A
- 40chol TYR13 / VAL17 / SER20 mean: 3.580 / 4.723 / 6.483 A
- 20chol key-min distance final/mean: 3.631/3.696 A
- 40chol key-min distance final/mean: 3.608/3.566 A

## Contact Occupancy 结果

- TYR13 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 85.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 95.0% / 100.0% / 100.0% / 100.0%
- VAL17 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 41.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 19.0% / 69.0% / 92.0% / 100.0%
- SER20 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 2.0% / 16.0% / 55.0% / 100.0%
  - 40chol: 0.0% / 1.0% / 22.0% / 99.0%
- key-min-distance `<4 / <5 / <6 A` occupancy
  - 20chol: 89.0% / 100.0% / 100.0%
  - 40chol: 98.0% / 100.0% / 100.0%
- near chains = 2 fraction
  - 20chol: 100.0%
  - 40chol: 100.0%
- ligand RMSD mean / median / max / 95th percentile
  - 20chol: 1.439 / 1.477 / 2.206 / 2.001 A
  - 40chol: 1.710 / 1.700 / 2.276 / 2.136 A

## 20chol vs 40chol 接触模式差异

- 20chol contact mode label: `TYR13_shifted` (contact shift)
- 40chol contact mode label: `TYR13_shifted` (contact shift)
- Cholesterol response label: `interface_retained_but_weaker_contact`
- Candidate interpretation: 2.3 / L006 keeps interface retention in both systems, but one or more key-residue occupancies soften in 40chol.

## 局限性

- 100 ps restrained MD 不是 production MD。
- 两个体系都使用了 restraints，因此这更适合作为 pilot stability / retention 检查，而不是平衡态结论。
- CGenFF 参数 penalty 仍然需要和本分析一起解读。
- 目前仍只覆盖单个 ligand、单条 trajectory。
- 后续仍需要更多 ligand、更长时间、更少 restraints，以及 replicate trajectories 支持。

## 下一步建议

1. 如果 2.3 双体系 occupancy 保持良好，可作为后续跨胆固醇横向比较的重点候选。
2. 后续可进入 14.2 vs 19.1 vs 2.3 的横向比较。
3. 暂时不要称为 production MD 结论。

## 输出

- Summary statistics table: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40\tables\comparison_summary_statistics.csv`
- Combined time series table: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40\tables\combined_time_series_2.3.csv`
- Contact occupancy summary table: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40\tables\contact_occupancy_summary.csv`
- Figures directory: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40\figures`
- PyMOL visualization script: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40\pymol\view_2.3_20_vs_40_final.pml`
