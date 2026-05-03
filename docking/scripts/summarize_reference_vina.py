from pathlib import Path
import csv
import re

root = Path(r"C:\TRKB_WP2\docking")
log_dir = root / "reference" / "logs"
out_csv = root / "reference" / "reference_vina_summary.csv"

score_pattern = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")

def read_log_text(path: Path) -> str:
    encodings = ["utf-8", "utf-16", "utf-16-le", "gbk", "cp936", "latin1"]
    for enc in encodings:
        try:
            return path.read_text(encoding=enc)
        except Exception:
            pass
    return path.read_bytes().decode(errors="ignore")

rows = []

for system in ["20chol", "40chol"]:
    for log_file in sorted(log_dir.glob(f"*_{system}.log")):
        ligand_id = log_file.name.replace(f"_{system}.log", "")
        best_affinity = None

        text = read_log_text(log_file)
        for line in text.splitlines():
            m = score_pattern.match(line)
            if m:
                best_affinity = float(m.group(1))
                break

        rows.append({
            "ligand_id": ligand_id,
            "system": system,
            "best_affinity_kcal_mol": best_affinity,
            "log_file": str(log_file),
        })

with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "ligand_id",
            "system",
            "best_affinity_kcal_mol",
            "log_file",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

print("Wrote:", out_csv)
print("Rows:", len(rows))

missing = [r for r in rows if r["best_affinity_kcal_mol"] is None]
print("Missing parsed scores:", len(missing))

if missing:
    print("\nFirst missing rows:")
    for r in missing[:10]:
        print(r["ligand_id"], r["system"], r["log_file"])

for system in ["20chol", "40chol"]:
    sub = [
        r for r in rows
        if r["system"] == system and r["best_affinity_kcal_mol"] is not None
    ]
    sub.sort(key=lambda r: r["best_affinity_kcal_mol"])

    print(f"\nReference ranking for {system}:")
    for r in sub:
        print(r["ligand_id"], r["best_affinity_kcal_mol"])
