from pathlib import Path
import csv
import math

root = Path(r"C:\TRKB_WP2\docking")

summary_csv = root / "reference" / "reference_vina_summary.csv"
out_csv = root / "reference" / "reference_contact_analysis.csv"

systems = {
    "20chol": {
        "receptor_pdb": root / "receptors" / "20chol" / "TRKB_20chol_protein.pdb",
        "pose_dir": root / "reference" / "results" / "20chol",
        "pose_suffix": "_20chol_out.pdbqt",
    },
    "40chol": {
        "receptor_pdb": root / "receptors" / "40chol" / "TRKB_40chol_protein.pdb",
        "pose_dir": root / "reference" / "results" / "40chol",
        "pose_suffix": "_40chol_out.pdbqt",
    },
}

target_resseq = {13, 17, 20}
contact_cutoff = 4.0
near_cutoff = 6.0


def dist(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def parse_receptor_target_atoms(pdb_path):
    atoms = []

    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip() or "NA"
                res_seq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue

            if res_seq in target_resseq:
                atoms.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain_id": chain_id,
                    "res_seq": res_seq,
                    "xyz": (x, y, z),
                })

    if not atoms:
        raise RuntimeError(f"No target residues 13/17/20 found in {pdb_path}")

    return atoms


def parse_first_model_ligand_atoms(pdbqt_path):
    atoms = []
    saw_model = False
    in_model = False

    with open(pdbqt_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                saw_model = True
                in_model = True
                continue

            if line.startswith("ENDMDL") and saw_model:
                break

            if saw_model and not in_model:
                continue

            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    continue
                atoms.append((x, y, z))

    if not atoms:
        raise RuntimeError(f"No ligand atoms found in {pdbqt_path}")

    return atoms


with open(summary_csv, newline="", encoding="utf-8-sig") as f:
    summary_rows = list(csv.DictReader(f))

target_atoms_by_system = {
    system: parse_receptor_target_atoms(info["receptor_pdb"])
    for system, info in systems.items()
}

out_rows = []

for row in summary_rows:
    ligand_id = row["ligand_id"]
    system = row["system"]
    affinity = row["best_affinity_kcal_mol"]

    info = systems[system]
    pose_file = info["pose_dir"] / f"{ligand_id}{info['pose_suffix']}"

    if not pose_file.exists():
        print("Missing pose:", pose_file)
        continue

    ligand_atoms = parse_first_model_ligand_atoms(pose_file)
    target_atoms = target_atoms_by_system[system]

    min_dist_all = 999.0
    contacted_residues = set()
    near_residues = set()
    contacted_chains = set()
    near_chains = set()

    for latom in ligand_atoms:
        for ratom in target_atoms:
            d = dist(latom, ratom["xyz"])
            key = f"{ratom['res_name']}{ratom['res_seq']}_chain{ratom['chain_id']}"

            min_dist_all = min(min_dist_all, d)

            if d <= contact_cutoff:
                contacted_residues.add(key)
                contacted_chains.add(ratom["chain_id"])

            if d <= near_cutoff:
                near_residues.add(key)
                near_chains.add(ratom["chain_id"])

    out_rows.append({
        "ligand_id": ligand_id,
        "system": system,
        "best_affinity_kcal_mol": affinity,
        "min_distance_to_key_residues_A": f"{min_dist_all:.3f}",
        "n_contacted_key_residues_4A": len(contacted_residues),
        "contacted_key_residues_4A": ";".join(sorted(contacted_residues)),
        "n_near_key_residues_6A": len(near_residues),
        "near_key_residues_6A": ";".join(sorted(near_residues)),
        "n_contacted_chains_4A": len(contacted_chains),
        "n_near_chains_6A": len(near_chains),
        "pose_file": str(pose_file),
    })

with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
    fieldnames = [
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
        "pose_file",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(out_rows)

print("Wrote:", out_csv)
print("Rows:", len(out_rows))

for system in ["20chol", "40chol"]:
    sub = [r for r in out_rows if r["system"] == system]
    sub.sort(
        key=lambda r: (
            -int(r["n_near_key_residues_6A"]),
            -int(r["n_near_chains_6A"]),
            float(r["best_affinity_kcal_mol"]) if r["best_affinity_kcal_mol"] else 999,
            float(r["min_distance_to_key_residues_A"]),
        )
    )

    print(f"\nReference contact ranking for {system}:")
    for r in sub:
        print(
            r["ligand_id"],
            "score=", r["best_affinity_kcal_mol"],
            "minD=", r["min_distance_to_key_residues_A"],
            "nearResidues6A=", r["n_near_key_residues_6A"],
            "nearChains6A=", r["n_near_chains_6A"],
            "near=", r["near_key_residues_6A"],
        )