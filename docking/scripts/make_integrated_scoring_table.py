from pathlib import Path
import csv
import statistics
from collections import defaultdict

ROOT = Path(r"C:\TRKB_WP2\docking")

BEST_CSV = ROOT / "best_pose_per_ligand.csv"
ALLMODES_CSV = ROOT / "all_modes_contact_analysis.csv"

OUT_SYSTEM_CSV = ROOT / "integrated_scoring_by_system.csv"
OUT_CANDIDATE_CSV = ROOT / "integrated_candidate_summary.csv"
OUT_REFERENCE_CSV = ROOT / "integrated_reference_summary.csv"

REFERENCE_LIGANDS = {
    "FLX_R",
    "FLX_S",
    "IMI",
    "HNK_RR",
    "HNK_SS",
    "KET_R",
    "KET_S",
    "VEN",
    "MOC",
    "CPZ",
    "DPH",
    "ISO",
}

ANCHOR_REFS = {
    "FLX_R",
    "FLX_S",
    "IMI",
    "HNK_RR",
}

COMPARATOR_REFS = REFERENCE_LIGANDS - ANCHOR_REFS


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def to_float(x, default=None):
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


def to_int(x, default=0):
    try:
        if x is None or x == "":
            return default
        return int(float(x))
    except Exception:
        return default


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def median_or_none(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return statistics.median(values)


def calc_mode_reliability(mode, affinity, mode1_affinity):
    """
    10分：最佳pose越靠前、与mode1能量差越小，越可靠。
    """
    mode = to_int(mode, 10)
    affinity = to_float(affinity, None)
    mode1_affinity = to_float(mode1_affinity, None)

    if mode <= 1:
        mode_rank_points = 5
    elif mode <= 3:
        mode_rank_points = 4
    elif mode <= 5:
        mode_rank_points = 3
    elif mode <= 7:
        mode_rank_points = 2
    else:
        mode_rank_points = 1

    if affinity is None or mode1_affinity is None:
        energy_gap_points = 2
    else:
        # affinity 越负越好。若最佳接触pose比mode1差太多，需要降权
        gap = affinity - mode1_affinity
        if gap <= 0.25:
            energy_gap_points = 5
        elif gap <= 0.50:
            energy_gap_points = 4
        elif gap <= 0.80:
            energy_gap_points = 3
        elif gap <= 1.20:
            energy_gap_points = 2
        else:
            energy_gap_points = 1

    return mode_rank_points + energy_gap_points


def classify(total, interface_flag, key_contact_flag):
    if total >= 75 and interface_flag and key_contact_flag:
        return "A_high_priority"
    if total >= 65 and interface_flag:
        return "B_good"
    if total >= 55:
        return "C_check_manually"
    return "D_low_priority"


def main():
    best_rows = read_csv(BEST_CSV)
    allmode_rows = read_csv(ALLMODES_CSV)

    # 获取每个 ligand-system 的 mode1 affinity，作为 mode 可靠性参考
    mode1_affinity = {}
    for r in allmode_rows:
        key = (r["group"], r["ligand_id"], r["system"])
        if to_int(r.get("mode")) == 1:
            mode1_affinity[key] = to_float(r.get("affinity_kcal_mol"))

    # 计算每个系统的 reference anchor
    ref_by_system = defaultdict(list)
    anchor_ref_by_system = defaultdict(list)
    comparator_ref_by_system = defaultdict(list)

    for r in best_rows:
        if r["group"] != "reference":
            continue
        ligand = r["ligand_id"]
        system = r["system"]
        aff = to_float(r.get("affinity_kcal_mol"))

        if aff is None:
            continue

        ref_by_system[system].append(aff)

        if ligand in ANCHOR_REFS:
            anchor_ref_by_system[system].append(aff)
        if ligand in COMPARATOR_REFS:
            comparator_ref_by_system[system].append(aff)

    anchors = {}
    for system in sorted(set(r["system"] for r in best_rows)):
        anchor_med = median_or_none(anchor_ref_by_system[system])
        comparator_med = median_or_none(comparator_ref_by_system[system])
        all_med = median_or_none(ref_by_system[system])

        anchors[system] = {
            "anchor_refs": ";".join(sorted(ANCHOR_REFS)),
            "comparator_refs": ";".join(sorted(COMPARATOR_REFS)),
            "anchor_median": anchor_med,
            "comparator_median": comparator_med,
            "all_reference_median": all_med,
        }

    scored_rows = []

    for r in best_rows:
        group = r["group"]
        ligand = r["ligand_id"]
        system = r["system"]

        aff = to_float(r.get("affinity_kcal_mol"))
        mode = to_int(r.get("mode"))
        minD = to_float(r.get("min_distance_to_key_residues_A"), 999.0)

        n_contact_res_4A = to_int(r.get("n_contacted_key_residues_4A"))
        n_near_res_6A = to_int(r.get("n_near_key_residues_6A"))
        n_contact_chains_4A = to_int(r.get("n_contacted_chains_4A"))
        n_near_chains_6A = to_int(r.get("n_near_chains_6A"))

        # 1. affinity score，25分
        # 使用用户指定 reference anchors：FLX_R / FLX_S / IMI / HNK_RR 的 median 作为较好锚点；
        # 其他 reference ligands 的 median 作为较弱/对照锚点。
        anchor_med = anchors[system]["anchor_median"]
        comparator_med = anchors[system]["comparator_median"]
        all_med = anchors[system]["all_reference_median"]

        # 如果 comparator median 缺失，使用 all reference median 作为较弱锚点。
        weak_anchor = comparator_med
        if weak_anchor is None:
            weak_anchor = all_med

        if aff is None or anchor_med is None or weak_anchor is None or abs(weak_anchor - anchor_med) < 1e-6:
            affinity_points = 12.5
        else:
            # affinity 越负越好
            # 若 aff 优于 anchor median，则接近满分
            # 若 aff 弱于 weak anchor，则接近 0
            raw = (weak_anchor - aff) / (weak_anchor - anchor_med)
            affinity_points = clamp(raw, 0.0, 1.2) / 1.2 * 25.0

        # 2. key residue contact score，30分
        # contact4A 更重要；near6A 作为补充；minD 作为几何合理性
        contact_points = clamp(n_contact_res_4A / 4.0, 0, 1) * 14.0
        near_points = clamp(n_near_res_6A / 5.0, 0, 1) * 10.0

        if minD <= 2.5:
            distance_points = 6.0
        elif minD <= 3.0:
            distance_points = 5.0
        elif minD <= 4.0:
            distance_points = 3.5
        elif minD <= 6.0:
            distance_points = 2.0
        else:
            distance_points = 0.0

        keyres_points = contact_points + near_points + distance_points

        # 3. interface / two-chain score，20分
        if n_near_chains_6A >= 2:
            near_chain_points = 8.0
        elif n_near_chains_6A == 1:
            near_chain_points = 3.0
        else:
            near_chain_points = 0.0

        if n_contact_chains_4A >= 2:
            contact_chain_points = 8.0
        elif n_contact_chains_4A == 1:
            contact_chain_points = 3.0
        else:
            contact_chain_points = 0.0

        # 如果同时满足接近两条链且至少有关键残基4A接触，给额外界面奖励
        interface_bonus = 4.0 if (n_near_chains_6A >= 2 and n_contact_res_4A >= 2) else 0.0
        interface_points = near_chain_points + contact_chain_points + interface_bonus

        # 4. mode reliability，10分
        key = (group, ligand, system)
        mode_reliability_points = calc_mode_reliability(
            mode=mode,
            affinity=aff,
            mode1_affinity=mode1_affinity.get(key),
        )

        # 5. reference role / sanity flag，不直接给候选加分，只标记
        if group == "reference" and ligand in ANCHOR_REFS:
            ref_role = "anchor_reference"
        elif group == "reference" and ligand in COMPARATOR_REFS:
            ref_role = "reference_comparator"
        elif group == "reference":
            ref_role = "reference_other"
        else:
            ref_role = "candidate"

        # system score 最高 85，剩余15分给跨环境一致性，在 candidate summary 中计算
        system_score_85 = (
            affinity_points
            + keyres_points
            + interface_points
            + mode_reliability_points
        )

        interface_flag = n_near_chains_6A >= 2
        key_contact_flag = n_contact_res_4A >= 2 or n_near_res_6A >= 4

        scored = dict(r)
        scored.update({
            "ref_role": ref_role,
            "anchor_ref_panel": anchors[system]["anchor_refs"],
            "comparator_ref_panel": anchors[system]["comparator_refs"],
            "anchor_ref_median_affinity": f"{anchor_med:.3f}" if anchor_med is not None else "",
            "comparator_ref_median_affinity": f"{comparator_med:.3f}" if comparator_med is not None else "",
            "mode1_affinity_kcal_mol": mode1_affinity.get(key, ""),
            "affinity": r.get("affinity_kcal_mol", ""),
            "min_distance_A": f"{minD:.3f}" if minD != 999.0 else "",
            "contacts_4A": n_contact_res_4A,
            "contacts_4A_detail": r.get("contacted_key_residues_4A", ""),
            "contacts_6A": n_near_res_6A,
            "contacts_6A_detail": r.get("near_key_residues_6A", ""),
            "chain_contacts": n_contact_chains_4A,
            "chain_near_6A": n_near_chains_6A,
            "affinity_points_25": f"{affinity_points:.2f}",
            "keyres_points_30": f"{keyres_points:.2f}",
            "interface_points_20": f"{interface_points:.2f}",
            "mode_reliability_points_10": f"{mode_reliability_points:.2f}",
            "system_score": f"{system_score_85:.2f}",
            "system_score_85": f"{system_score_85:.2f}",
            "interface_flag": int(interface_flag),
            "key_contact_flag": int(key_contact_flag),
        })

        scored_rows.append(scored)

    # 写 system-level 综合评分
    system_fieldnames = [
        "group",
        "ref_role",
        "ligand_id",
        "system",
        "mode",
        "affinity",
        "affinity_kcal_mol",
        "mode1_affinity_kcal_mol",
        "rmsd_lb",
        "rmsd_ub",
        "min_distance_A",
        "min_distance_to_key_residues_A",
        "contacts_4A",
        "contacts_4A_detail",
        "n_contacted_key_residues_4A",
        "contacted_key_residues_4A",
        "contacts_6A",
        "contacts_6A_detail",
        "n_near_key_residues_6A",
        "near_key_residues_6A",
        "chain_contacts",
        "chain_near_6A",
        "n_contacted_chains_4A",
        "n_near_chains_6A",
        "anchor_ref_panel",
        "comparator_ref_panel",
        "anchor_ref_median_affinity",
        "comparator_ref_median_affinity",
        "affinity_points_25",
        "keyres_points_30",
        "interface_points_20",
        "mode_reliability_points_10",
        "system_score",
        "system_score_85",
        "interface_flag",
        "key_contact_flag",
        "pose_file",
    ]

    with open(OUT_SYSTEM_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=system_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(scored_rows)

    # candidate 跨 20/40 汇总
    cand_rows = [r for r in scored_rows if r["group"] == "candidate"]
    by_ligand = defaultdict(dict)

    for r in cand_rows:
        by_ligand[r["ligand_id"]][r["system"]] = r

    candidate_summary = []

    for ligand, sysmap in by_ligand.items():
        r20 = sysmap.get("20chol")
        r40 = sysmap.get("40chol")

        s20 = to_float(r20.get("system_score_85")) if r20 else None
        s40 = to_float(r40.get("system_score_85")) if r40 else None

        # 跨环境一致性 15分
        # 两个环境都有好表现时加分；只在一个环境好则中等
        robust_points = 0.0
        if s20 is not None and s40 is not None:
            min_s = min(s20, s40)
            avg_s = (s20 + s40) / 2.0

            if min_s >= 60:
                robust_points = 15.0
            elif min_s >= 52:
                robust_points = 12.0
            elif avg_s >= 60:
                robust_points = 9.0
            elif avg_s >= 52:
                robust_points = 6.0
            else:
                robust_points = 3.0
        elif s20 is not None or s40 is not None:
            robust_points = 4.0

        avg_system = statistics.mean([x for x in [s20, s40] if x is not None])
        min_system = min([x for x in [s20, s40] if x is not None])

        final_score_100 = avg_system + robust_points

        best_system = ""
        best_mode = ""
        best_aff = ""
        if r20 and r40:
            if s20 >= s40:
                best = r20
            else:
                best = r40
        else:
            best = r20 or r40

        if best:
            best_system = best["system"]
            best_mode = best["mode"]
            best_aff = best["affinity_kcal_mol"]

        # 分类
        interface_both = (
            r20 is not None and r40 is not None
            and to_int(r20.get("interface_flag")) == 1
            and to_int(r40.get("interface_flag")) == 1
        )
        keycontact_both = (
            r20 is not None and r40 is not None
            and to_int(r20.get("key_contact_flag")) == 1
            and to_int(r40.get("key_contact_flag")) == 1
        )

        priority = classify(final_score_100, interface_both, keycontact_both)

        candidate_summary.append({
            "ligand_id": ligand,
            "final_score_100": f"{final_score_100:.2f}",
            "priority_class": priority,
            "robustness_points_15": f"{robust_points:.2f}",
            "avg_system_score_85": f"{avg_system:.2f}",
            "min_system_score_85": f"{min_system:.2f}",
            "score_20chol": f"{s20:.2f}" if s20 is not None else "",
            "score_20chol_85": f"{s20:.2f}" if s20 is not None else "",
            "mode_20chol": r20["mode"] if r20 else "",
            "affinity_20chol": r20["affinity_kcal_mol"] if r20 else "",
            "contact4A_20chol": r20["n_contacted_key_residues_4A"] if r20 else "",
            "near6A_20chol": r20["n_near_key_residues_6A"] if r20 else "",
            "nearChains_20chol": r20["n_near_chains_6A"] if r20 else "",
            "score_40chol": f"{s40:.2f}" if s40 is not None else "",
            "score_40chol_85": f"{s40:.2f}" if s40 is not None else "",
            "mode_40chol": r40["mode"] if r40 else "",
            "affinity_40chol": r40["affinity_kcal_mol"] if r40 else "",
            "contact4A_40chol": r40["n_contacted_key_residues_4A"] if r40 else "",
            "near6A_40chol": r40["n_near_key_residues_6A"] if r40 else "",
            "nearChains_40chol": r40["n_near_chains_6A"] if r40 else "",
            "best_system": best_system,
            "best_mode": best_mode,
            "best_affinity": best_aff,
        })

    candidate_summary.sort(key=lambda r: to_float(r["final_score_100"], 0), reverse=True)

    candidate_fieldnames = [
        "ligand_id",
        "final_score_100",
        "priority_class",
        "robustness_points_15",
        "avg_system_score_85",
        "min_system_score_85",
        "score_20chol",
        "score_20chol_85",
        "mode_20chol",
        "affinity_20chol",
        "contact4A_20chol",
        "near6A_20chol",
        "nearChains_20chol",
        "score_40chol",
        "score_40chol_85",
        "mode_40chol",
        "affinity_40chol",
        "contact4A_40chol",
        "near6A_40chol",
        "nearChains_40chol",
        "best_system",
        "best_mode",
        "best_affinity",
    ]

    with open(OUT_CANDIDATE_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=candidate_fieldnames)
        writer.writeheader()
        writer.writerows(candidate_summary)

    # reference 汇总，方便验证
    ref_summary = [r for r in scored_rows if r["group"] == "reference"]
    ref_summary.sort(
        key=lambda r: (
            r["system"],
            -to_float(r["system_score_85"], 0)
        )
    )

    with open(OUT_REFERENCE_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=system_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ref_summary)

    print("Wrote:", OUT_SYSTEM_CSV)
    print("Wrote:", OUT_CANDIDATE_CSV)
    print("Wrote:", OUT_REFERENCE_CSV)

    print("\nReference anchors used for scoring:")
    print("anchor_refs:", ", ".join(sorted(ANCHOR_REFS)))
    print("comparator_refs:", ", ".join(sorted(COMPARATOR_REFS)))
    for system, a in anchors.items():
        anchor_med = a["anchor_median"]
        comparator_med = a["comparator_median"]
        all_med = a["all_reference_median"]
        print(
            system,
            "anchor_median=", f"{anchor_med:.3f}" if anchor_med is not None else "",
            "comparator_median=", f"{comparator_med:.3f}" if comparator_med is not None else "",
            "all_reference_median=", f"{all_med:.3f}" if all_med is not None else "",
        )

    print("\nTop 20 candidates by final integrated score:")
    for r in candidate_summary[:20]:
        print(
            r["ligand_id"],
            "final=", r["final_score_100"],
            "class=", r["priority_class"],
            "20=", r["score_20chol"],
            "mode20=", r["mode_20chol"],
            "40=", r["score_40chol"],
            "mode40=", r["mode_40chol"],
            "best=", r["best_system"], "mode", r["best_mode"],
        )

    print("\nSuggested Top 6 for manual PyMOL review:")
    for r in candidate_summary[:6]:
        print(r["ligand_id"], r["final_score_100"], r["priority_class"])


if __name__ == "__main__":
    main()
