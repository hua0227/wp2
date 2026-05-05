# 9.2 / L004 20chol vs 40chol R2 Pilot Comparison

## 分析目的

比较 9.2 / L004 在 20% CHOL 与 40% CHOL TRKB-TMD 体系中的 restrained MD pilot 稳定性。

## 输入数据

- 20chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md100ps_R2\reports\md100ps_20chol_R2_monitor.csv`
- 20chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md100ps_R2\outputs\TRKB_20chol_L004_9.2_100ps_R2.dcd`
- 20chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md100ps_R2\outputs\TRKB_20chol_L004_9.2_100ps_R2_final.pdb`
- 20chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_20chol\md100ps_R2\reports\md100ps_20chol_R2_summary.md`
- 40chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv`
- 40chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md100ps_R2\outputs\TRKB_40chol_L004_9.2_100ps_R2.dcd`
- 40chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md100ps_R2\outputs\TRKB_40chol_L004_9.2_100ps_R2_final.pdb`
- 40chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\9.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md`

## 模拟状态说明

两个体系都属于 100 ps R2 PBC-aware restrained MD pilot。它们不是 production MD 模拟，本报告也不是 production MD 结论。

## 稳定性结果

- 20chol completed ps: 100.000; 40chol completed ps: 100.000
- 20chol temperature range: 306.879-313.012 K
- 40chol temperature range: 306.928-314.071 K
- NaN detected, 20chol/40chol: NO / NO
- severe clash ever, 20chol/40chol: NO / NO
- close clash ever, 20chol/40chol: NO / NO
- near chains remained 2 throughout, 20chol/40chol: YES / YES

## Ligand 保留情况

- 20chol ligand RMSD final/max: 1.005/2.631 A
- 40chol ligand RMSD final/max: 1.110/1.578 A
- 20chol TYR13 / VAL17 / SER20 final: 4.088 / 3.700 / 4.179 A
- 40chol TYR13 / VAL17 / SER20 final: 3.358 / 6.191 / 6.200 A
- 20chol TYR13 / VAL17 / SER20 mean: 3.900 / 3.956 / 5.128 A
- 40chol TYR13 / VAL17 / SER20 mean: 3.497 / 4.834 / 4.838 A
- 20chol key-min distance final/mean: 3.700/3.694 A
- 40chol key-min distance final/mean: 3.358/3.489 A

## Contact Occupancy 结果

- TYR13 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 60.0% / 98.0% / 100.0% / 100.0%
  - 40chol: 98.0% / 100.0% / 100.0% / 100.0%
- VAL17 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 62.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 11.0% / 62.0% / 94.0% / 100.0%
- SER20 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 13.0% / 42.0% / 85.0% / 100.0%
  - 40chol: 8.0% / 65.0% / 96.0% / 100.0%
- key-min-distance `<4 / <5 / <6 A` occupancy
  - 20chol: 87.0% / 100.0% / 100.0%
  - 40chol: 98.0% / 100.0% / 100.0%
- near chains = 2 fraction
  - 20chol: 100.0%
  - 40chol: 100.0%
- ligand RMSD mean / median / max / 95th percentile
  - 20chol: 1.660 / 1.643 / 2.631 / 2.444 A
  - 40chol: 0.922 / 0.919 / 1.578 / 1.315 A

## 20chol vs 40chol 接触模式差异

- 20chol contact mode label: `weak_or_mixed_contact` (interface retention only)
- 40chol contact mode label: `weak_or_mixed_contact` (interface retention only)
- Cholesterol response label: `interface_retained_but_weaker_contact`
- Candidate interpretation: 9.2 / L004 keeps interface retention in both systems, but one or more key-residue occupancies soften in 40chol.

## 局限性

- 100 ps restrained MD 不是 production MD。
- 两个体系都使用了 restraints，因此这更适合作为 pilot stability / retention 检查，而不是平衡态结论。
- CGenFF 参数 penalty 仍然需要和本分析一起解读。
- 目前仍只覆盖单个 ligand、单条 trajectory。
- 后续仍需要更多 ligand、更长时间、更少 restraints，以及 replicate trajectories 支持。

## 下一步建议

1. 如果 9.2 双体系 occupancy 保持良好，可作为后续跨胆固醇横向比较的重点候选。
2. 后续可进入 14.2 vs 19.1 vs 9.2 的横向比较。
3. 暂时不要称为 production MD 结论。

## 输出

- Summary statistics table: `C:\TRKB_WP2\ligand_bound_MD\comparison_9.2_20_vs_40\tables\comparison_summary_statistics.csv`
- Combined time series table: `C:\TRKB_WP2\ligand_bound_MD\comparison_9.2_20_vs_40\tables\combined_time_series_9.2.csv`
- Contact occupancy summary table: `C:\TRKB_WP2\ligand_bound_MD\comparison_9.2_20_vs_40\tables\contact_occupancy_summary.csv`
- Figures directory: `C:\TRKB_WP2\ligand_bound_MD\comparison_9.2_20_vs_40\figures`
- PyMOL visualization script: `C:\TRKB_WP2\ligand_bound_MD\comparison_9.2_20_vs_40\pymol\view_9.2_20_vs_40_final.pml`
