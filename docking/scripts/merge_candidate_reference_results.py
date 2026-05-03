from pathlib import Path
import csv

root = Path(r"C:\TRKB_WP2\docking")

candidate_csv = root / "vina_contact_analysis.csv"
reference_csv = root / "reference" / "reference_contact_analysis.csv"
out_csv = root / "merged_candidate_reference_ranking.csv"

rows = []

with open(candidate_csv, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for r in reader:
        r["group"] = "candidate"
        rows.append(r)

with open(reference_csv, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for r in reader:
        r["group"] = "reference"
        rows.append(r)

def to_float(x, default=999.0):
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
        return int(x)
    except Exception:
        return default

# 简单综合排序分：
# 越高越好
# 主要奖励：score 强、靠近关键残基、靠近更多关键残基、接触两条链
for r in rows:
    score = to_float(r.get("best_affinity_kcal_mol"))
    minD = to_float(r.get("min_distance_to_key_residues_A"))
    nearRes = to_int(r.get("n_near_key_residues_6A"))
    contactRes = to_int(r.get("n_contacted_key_residues_4A"))
    nearChains = to_int(r.get("n_near_chains_6A"))
    contactChains = to_int(r.get("n_contacted_chains_4A"))

    composite = 0.0
    composite += -score * 1.0
    composite += nearRes * 0.6
    composite += contactRes * 0.8
    composite += nearChains * 1.0
    composite += contactChains * 1.2
    composite += max(0.0, 6.0 - minD) * 0.3

    r["composite_score"] = f"{composite:.3f}"

fieldnames = [
    "group",
    "ligand_id",
    "system",
    "best_affinity_kcal_mol",
    "min_distance_to_key_residues_A",
    "n_contacted_key_residues_4A",
    "contacted_key_residues_4A",
    "n_near_key_residues_6A",
    "near_key_residues_6A",
    "n_contacted_chains_4A",
    "n_near_chains_6A",
    "composite_score",
    "pose_file",
]

with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print("Wrote:", out_csv)
print("Rows:", len(rows))

for system in ["20chol", "40chol"]:
    sub = [r for r in rows if r["system"] == system]
    sub.sort(key=lambda r: to_float(r["composite_score"], default=-999), reverse=True)

    print(f"\nTop 20 merged ranking for {system}:")
    for r in sub[:20]:
        print(
            r["group"],
            r["ligand_id"],
            "score=", r["best_affinity_kcal_mol"],
            "minD=", r["min_distance_to_key_residues_A"],
            "contact4A=", r["n_contacted_key_residues_4A"],
            "near6A=", r["n_near_key_residues_6A"],
            "nearChains=", r["n_near_chains_6A"],
            "composite=", r["composite_score"],
        )