# 19.1 / L003 20chol vs 40chol R2 Pilot Comparison

## 分析目的

比较 19.1 / L003 在 20% CHOL 与 40% CHOL TRKB-TMD 体系中的 restrained MD pilot 稳定性。

## 输入数据

- 20chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\reports\md100ps_20chol_R2_monitor.csv`
- 20chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\outputs\TRKB_20chol_L003_19.1_100ps_R2.dcd`
- 20chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\outputs\TRKB_20chol_L003_19.1_100ps_R2_final.pdb`
- 20chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\reports\md100ps_20chol_R2_summary.md`
- 40chol monitor CSV: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv`
- 40chol DCD: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md100ps_R2\outputs\TRKB_40chol_L003_19.1_100ps_R2.dcd`
- 40chol final PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md100ps_R2\outputs\TRKB_40chol_L003_19.1_100ps_R2_final.pdb`
- 40chol summary: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md`

## 模拟状态

两个体系都属于 100 ps R2 PBC-aware restrained MD pilot。它们不是 production MD 模拟，本报告也不是 production MD 结论。

## 稳定性结果

- 20chol completed ps: 100.000; 40chol completed ps: 100.000
- 20chol temperature range: 307.396-313.033 K
- 40chol temperature range: 307.428-313.780 K
- NaN detected, 20chol/40chol: False / False
- severe clash ever, 20chol/40chol: False / False
- close clash ever, 20chol/40chol: False / False
- near chains remained 2 throughout, 20chol/40chol: True / True

## Ligand 保留情况

- 20chol ligand RMSD final/max: 1.882/2.487 A
- 40chol ligand RMSD final/max: 1.551/2.444 A
- 20chol TYR13 / VAL17 / SER20 final: 3.775 / 3.191 / 3.654 A
- 40chol TYR13 / VAL17 / SER20 final: 3.745 / 4.085 / 3.791 A
- 20chol TYR13 / VAL17 / SER20 mean: 3.877 / 3.817 / 3.641 A
- 40chol TYR13 / VAL17 / SER20 mean: 4.026 / 3.949 / 3.658 A
- 20chol key-min distance final/mean: 3.191/3.514 A
- 40chol key-min distance final/mean: 3.745/3.587 A

## Contact Occupancy 结果

- TYR13 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 73.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 55.0% / 97.0% / 100.0% / 100.0%
- VAL17 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 70.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 58.0% / 99.0% / 100.0% / 100.0%
- SER20 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: 94.0% / 100.0% / 100.0% / 100.0%
  - 40chol: 91.0% / 100.0% / 100.0% / 100.0%
- key-min-distance `<4 / <5 / <6 A` occupancy
  - 20chol: 99.0% / 100.0% / 100.0%
  - 40chol: 97.0% / 100.0% / 100.0%
- near chains = 2 fraction
  - 20chol: 100.0%
  - 40chol: 100.0%
- ligand RMSD mean / median / max / 95th percentile
  - 20chol: 1.692 / 1.744 / 2.487 / 2.204 A
  - 40chol: 1.417 / 1.350 / 2.444 / 2.181 A

## 关键比较

- 19.1 在 20chol 和 40chol 中都保持了 near chains = 2，并且 key-min distance occupancy 保持很高，说明 ligand 在两种胆固醇条件下都维持 interface retention。
- 与 14.2 的 40chol 模式不同，19.1 在 40chol 中并没有明显退化成偏 TYR13、同时失去 VAL17/SER20 的接触构象。
- 19.1 是否在 40chol 中仍保持 VAL17/SER20 contact: YES
- 结合最终距离、平均距离和 `<5 A` occupancy，19.1 看起来在两个胆固醇条件下都更均衡地保持了 TYR13 / VAL17 / SER20 三关键残基接触。

## 局限性

- 100 ps restrained MD 不是 production MD。
- 两个体系都使用了 restraints，因此这更适合作为 pilot stability / retention 检查，而不是平衡态结论。
- CGenFF 参数 penalty 对 L003 / 19.1 的解释仍然重要，需要与本分析一起看。
- 目前仍只覆盖单个 ligand、单条 trajectory。
- 后续仍需要更多 ligand、更长时间、更少 restraints，以及 replicate trajectories 支持。

## 下一步建议

1. 如果 19.1 双体系 occupancy 保持良好，可作为比 14.2 更均衡的重点候选。
2. 后续可进入 14.2 vs 19.1 横向比较。
3. 暂时不要称为 production MD 结论。

## 输出

- Summary statistics table: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40\tables\comparison_summary_statistics.csv`
- Combined time series table: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40\tables\combined_time_series_19.1.csv`
- Contact occupancy summary table: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40\tables\contact_occupancy_summary.csv`
- Figures directory: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40\figures`
- PyMOL visualization script: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40\pymol\view_19.1_20_vs_40_final.pml`
