# WP2 Gate 1 Summary Report

## A. 分析目的

总结五个 moderate-penalty 候选物在 20% CHOL 和 40% CHOL TRKB-TMD 体系中的 100 ps R2 PBC-aware restrained MD pilot 表现。

## B. 输入数据

- 14.2: `C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40`
- 19.1: `C:\TRKB_WP2\ligand_bound_MD\comparison_19.1_20_vs_40`
- 2.3: `C:\TRKB_WP2\ligand_bound_MD\comparison_2.3_20_vs_40`
- 6.2: `C:\TRKB_WP2\ligand_bound_MD\comparison_6.2_20_vs_40`
- 17.2: `C:\TRKB_WP2\ligand_bound_MD\comparison_17.2_20_vs_40`

主汇总使用了每个 comparison 目录中的 `comparison_summary_statistics.csv` 与 `contact_occupancy_summary.csv`。

## C. 方法说明

所有候选物均经过：
- ligand-bound topology build
- OpenMM read test
- targeted cleanup if needed
- R2 PBC-aware minimization
- very-short probe
- 10 ps / 50 ps / 100 ps restrained MD pilot
- contact occupancy analysis

本汇总只做离线整合，不新增任何 MD。本报告中的所有结论都只适用于 100 ps non-production R2 PBC-aware restrained MD pilot，不是 production MD 结论。

## D. 稳定性总结

- 五个候选物的 20chol / 40chol 体系都完成了 100 ps。
- 本轮汇总中没有检测到 NaN。
- temperature range 整体维持在约 307-314 K 的 pilot 区间。
- near chains 在 10 个 candidate-system 中都保持为 2。
- severe / close clash 在五个候选物的 comparison 表中都没有再次出现。

## E. 接触模式总结

- 14.2 / L002: 20chol `VAL17_SER20_dominant` (0.0% / 100.0% / 98.0%); 40chol `TYR13_shifted` (96.0% / 0.0% / 0.0%); response `cholesterol_sensitive_contact_shift`.
- 19.1 / L003: 20chol `balanced_three_residue_contact` (100.0% / 100.0% / 100.0%); 40chol `balanced_three_residue_contact` (97.0% / 99.0% / 100.0%); response `cholesterol_robust_balanced`.
- 2.3 / L006: 20chol `TYR13_shifted` (100.0% / 100.0% / 16.0%); 40chol `TYR13_shifted` (100.0% / 69.0% / 1.0%); response `TYR13_dominant_interface_retention`.
- 6.2 / L008: 20chol `weak_or_mixed_contact` (99.0% / 65.0% / 100.0%); 40chol `TYR13_shifted` (100.0% / 0.0% / 0.0%); response `cholesterol_sensitive_contact_shift`.
- 17.2 / L010: 20chol `balanced_three_residue_contact` (99.0% / 100.0% / 100.0%); 40chol `weak_or_mixed_contact` (100.0% / 98.0% / 63.0%); response `interface_retained_but_weaker_contact`.

## F. Gate 1 分类

- Priority A: 19.1 / L003, 17.2 / L010
- Priority B: 14.2 / L002
- Priority C: 2.3 / L006, 6.2 / L008
- Deferred / Reserve: none

## G. 推荐候选物

- Priority A 候选物适合作为下一阶段最优先稳定候选。
- Priority B 候选物更适合作为 cholesterol-sensitive contact-shift 机制候选。
- Priority C 候选物仍可保留，但不建议在 Gate 2 前就继续大规模扩展。

## H. 局限性

- 这些都是 100 ps non-production restrained MD pilot。
- 不是 production MD。
- 体系使用了 R2 PBC-aware restraints。
- 每个体系目前只有单次轨迹，没有 replicate。
- CGenFF penalty 仍需与本表一起解读。
- 仍需要更长时间、更弱 restraints 或重复轨迹验证。

## I. 下一步建议

1. 将 Priority A 候选物作为下一阶段优先稳定候选。
2. 将 Priority B 候选物作为 cholesterol-sensitive contact-shift 机制候选。
3. 将 Priority C 候选物保留为 reserve candidates。
4. Gate 1 之后建议先暂停继续大规模 MD，先给老师评估。
5. 若进入 Gate 2，可对 Priority A/B 做更长时间、更弱 restraints 或 replicate trajectories。

## 输出

- Candidate summary table: `C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\tables\WP2_Gate1_candidate_summary.csv`
- Ligand-level summary table: `C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\tables\WP2_Gate1_ligand_level_summary.csv`
- Report-ready summary table: `C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\tables\WP2_Gate1_report_ready_summary.csv`
- Figures directory: `C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\figures`
- PyMOL script: `C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\pymol\view_gate1_final_structures.pml`
