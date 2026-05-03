from pathlib import Path
import csv
import math
import re

root = Path(r"C:\TRKB_WP2\docking")

out_all_csv = root / "all_modes_contact_analysis.csv"
out_best_csv = root / "best_pose_per_ligand.csv"
best_pdbqt_dir = root / "best_pose_review_pdbqt"
pdb_dir = root / "best_pose_review_pdb"
pml_dir = root / "figure_best" / "pml"
view_dir = root / "figure_best" / "views"
render_jobs_csv = root / "figure_best" / "render_jobs_best.csv"

best_pdbqt_dir.mkdir(parents=True, exist_ok=True)
pdb_dir.mkdir(parents=True, exist_ok=True)
pml_dir.mkdir(parents=True, exist_ok=True)
view_dir.mkdir(parents=True, exist_ok=True)

target_resseq = {13, 17, 20}
contact_cutoff = 4.0
near_cutoff = 6.0

datasets = [
    {
        "group": "candidate",
        "system": "20chol",
        "receptor_pdb": root / "receptors" / "20chol" / "TRKB_20chol_protein.pdb",
        "pose_dir": root / "results_batch" / "20chol",
        "pose_suffix": "_20chol_out.pdbqt",
    },
    {
        "group": "candidate",
        "system": "40chol",
        "receptor_pdb": root / "receptors" / "40chol" / "TRKB_40chol_protein.pdb",
        "pose_dir": root / "results_batch" / "40chol",
        "pose_suffix": "_40chol_out.pdbqt",
    },
    {
        "group": "reference",
        "system": "20chol",
        "receptor_pdb": root / "receptors" / "20chol" / "TRKB_20chol_protein.pdb",
        "pose_dir": root / "reference" / "results" / "20chol",
        "pose_suffix": "_20chol_out.pdbqt",
    },
    {
        "group": "reference",
        "system": "40chol",
        "receptor_pdb": root / "receptors" / "40chol" / "TRKB_40chol_protein.pdb",
        "pose_dir": root / "reference" / "results" / "40chol",
        "pose_suffix": "_40chol_out.pdbqt",
    },
]

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

def parse_all_models_from_pdbqt(pdbqt_path):
    models = []
    current = None

    vina_remark = re.compile(
        r"REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)"
    )

    with open(pdbqt_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                parts = line.split()
                try:
                    mode = int(parts[1])
                except Exception:
                    mode = len(models) + 1

                current = {
                    "mode": mode,
                    "affinity": None,
                    "rmsd_lb": None,
                    "rmsd_ub": None,
                    "atoms": [],
                    "lines": [line],
                }
                continue

            if current is None:
                current = {
                    "mode": 1,
                    "affinity": None,
                    "rmsd_lb": None,
                    "rmsd_ub": None,
                    "atoms": [],
                    "lines": [],
                }

            current["lines"].append(line)

            m = vina_remark.search(line)
            if m:
                current["affinity"] = float(m.group(1))
                current["rmsd_lb"] = float(m.group(2))
                current["rmsd_ub"] = float(m.group(3))
                continue

            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    current["atoms"].append((x, y, z))
                except Exception:
                    continue

            if line.startswith("ENDMDL"):
                if current and current["atoms"]:
                    models.append(current)
                current = None

    if current and current["atoms"]:
        models.append(current)

    if not models:
        raise RuntimeError(f"No models found in {pdbqt_path}")

    return models

def analyze_model(ligand_atoms, target_atoms):
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

    return {
        "min_distance_to_key_residues_A": min_dist_all,
        "n_contacted_key_residues_4A": len(contacted_residues),
        "contacted_key_residues_4A": ";".join(sorted(contacted_residues)),
        "n_near_key_residues_6A": len(near_residues),
        "near_key_residues_6A": ";".join(sorted(near_residues)),
        "n_contacted_chains_4A": len(contacted_chains),
        "n_near_chains_6A": len(near_chains),
    }

def composite(row):
    affinity = row["affinity_kcal_mol"]
    minD = row["min_distance_to_key_residues_A"]
    nearRes = row["n_near_key_residues_6A"]
    contactRes = row["n_contacted_key_residues_4A"]
    nearChains = row["n_near_chains_6A"]
    contactChains = row["n_contacted_chains_4A"]

    val = 0.0
    if affinity is not None:
        val += -float(affinity) * 1.0
    val += int(nearRes) * 0.6
    val += int(contactRes) * 0.8
    val += int(nearChains) * 1.0
    val += int(contactChains) * 1.2
    val += max(0.0, 6.0 - float(minD)) * 0.3
    return val

def safe_label(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)

rows = []
model_cache = {}

target_atoms_cache = {}
for ds in datasets:
    system = ds["system"]
    if system not in target_atoms_cache:
        target_atoms_cache[system] = parse_receptor_target_atoms(ds["receptor_pdb"])

for ds in datasets:
    group = ds["group"]
    system = ds["system"]
    pose_dir = ds["pose_dir"]
    suffix = ds["pose_suffix"]
    target_atoms = target_atoms_cache[system]

    pose_files = sorted(pose_dir.glob(f"*{suffix}"))
    print(f"{group} {system}: pose files = {len(pose_files)}")

    for pose_file in pose_files:
        ligand_id = pose_file.name.replace(suffix, "")
        models = parse_all_models_from_pdbqt(pose_file)
        model_cache[str(pose_file)] = models

        for model in models:
            result = analyze_model(model["atoms"], target_atoms)

            row = {
                "group": group,
                "ligand_id": ligand_id,
                "system": system,
                "mode": model["mode"],
                "affinity_kcal_mol": model["affinity"],
                "rmsd_lb": model["rmsd_lb"],
                "rmsd_ub": model["rmsd_ub"],
                "min_distance_to_key_residues_A": round(result["min_distance_to_key_residues_A"], 3),
                "n_contacted_key_residues_4A": result["n_contacted_key_residues_4A"],
                "contacted_key_residues_4A": result["contacted_key_residues_4A"],
                "n_near_key_residues_6A": result["n_near_key_residues_6A"],
                "near_key_residues_6A": result["near_key_residues_6A"],
                "n_contacted_chains_4A": result["n_contacted_chains_4A"],
                "n_near_chains_6A": result["n_near_chains_6A"],
                "pose_file": str(pose_file),
            }
            row["composite_score"] = round(composite(row), 3)
            rows.append(row)

fieldnames = [
    "group",
    "ligand_id",
    "system",
    "mode",
    "affinity_kcal_mol",
    "rmsd_lb",
    "rmsd_ub",
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

with open(out_all_csv, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Wrote:", out_all_csv)
print("All mode rows:", len(rows))

best_rows = []

for group in sorted(set(r["group"] for r in rows)):
    for system in sorted(set(r["system"] for r in rows)):
        ligand_ids = sorted(set(r["ligand_id"] for r in rows if r["group"] == group and r["system"] == system))
        for ligand_id in ligand_ids:
            sub = [
                r for r in rows
                if r["group"] == group and r["system"] == system and r["ligand_id"] == ligand_id
            ]
            sub.sort(
                key=lambda r: (
                    float(r["composite_score"]),
                    -float(r["affinity_kcal_mol"]) if r["affinity_kcal_mol"] is not None else -999,
                    -float(r["min_distance_to_key_residues_A"]),
                ),
                reverse=True,
            )
            best_rows.append(sub[0])

with open(out_best_csv, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(best_rows)

print("Wrote:", out_best_csv)
print("Best pose rows:", len(best_rows))

# Select visualization jobs:
# 1. all reference best poses
# 2. top 12 candidate best poses per system
selected = []

for r in best_rows:
    if r["group"] == "reference":
        selected.append(r)

for system in ["20chol", "40chol"]:
    cand = [r for r in best_rows if r["group"] == "candidate" and r["system"] == system]
    cand.sort(key=lambda r: float(r["composite_score"]), reverse=True)
    selected.extend(cand[:12])

# Extract selected best mode as individual PDBQT
render_rows = []

for r in selected:
    pose_file = Path(r["pose_file"])
    models = model_cache.get(str(pose_file))
    if models is None:
        models = parse_all_models_from_pdbqt(pose_file)

    mode = int(r["mode"])
    chosen = None
    for m in models:
        if int(m["mode"]) == mode:
            chosen = m
            break

    if chosen is None:
        print("Could not find chosen mode:", r)
        continue

    label = safe_label(f"{r['group']}_{r['ligand_id']}_{r['system']}_mode{mode}")
    out_pdbqt = best_pdbqt_dir / f"{label}.pdbqt"
    out_pdb = pdb_dir / f"{label}.pdb"

    with open(out_pdbqt, "w", encoding="utf-8") as f:
        for line in chosen["lines"]:
            f.write(line)

    receptor_pdb = root / "receptors" / r["system"] / f"TRKB_{r['system'].replace('chol','chol')}_protein.pdb"
    if r["system"] == "20chol":
        receptor_pdb = root / "receptors" / "20chol" / "TRKB_20chol_protein.pdb"
    else:
        receptor_pdb = root / "receptors" / "40chol" / "TRKB_40chol_protein.pdb"

    render_rows.append({
        "label": label,
        "system": r["system"],
        "receptor_pdb": receptor_pdb.as_posix(),
        "ligand_pdbqt": out_pdbqt.as_posix(),
        "ligand_pdb": out_pdb.as_posix(),
        "group": r["group"],
        "ligand_id": r["ligand_id"],
        "mode": mode,
        "affinity_kcal_mol": r["affinity_kcal_mol"],
        "composite_score": r["composite_score"],
    })

with open(render_jobs_csv, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "label",
            "system",
            "receptor_pdb",
            "ligand_pdbqt",
            "ligand_pdb",
            "group",
            "ligand_id",
            "mode",
            "affinity_kcal_mol",
            "composite_score",
        ],
    )
    writer.writeheader()
    writer.writerows(render_rows)

print("Wrote:", render_jobs_csv)
print("Render jobs:", len(render_rows))

# Generate PML files
def make_pml(row):
    label = row["label"]
    receptor = row["receptor_pdb"]
    ligand = row["ligand_pdb"]

    front = (view_dir / f"{label}_front.png").as_posix()
    front45 = (view_dir / f"{label}_front45.png").as_posix()
    side = (view_dir / f"{label}_side.png").as_posix()
    back = (view_dir / f"{label}_back.png").as_posix()

    pml = f"""
reinitialize
bg_color white
set ray_opaque_background, off
set opaque_background, off
set orthoscopic, on
set antialias, 2
set depth_cue, off
set ray_trace_mode, 1
set ray_shadows, off
set specular, 0.20
set ambient, 0.45
set direct, 0.55
set two_sided_lighting, on
set cartoon_fancy_helices, on
set cartoon_transparency, 0.12
set stick_radius, 0.18
set sphere_scale, 0.22
set valence, 0
set label_size, 18
set label_color, black
set label_outline_color, white
set label_shadow_mode, 0

load {receptor}, rec
load {ligand}, lig

hide everything, all
show cartoon, rec
color cyan, rec
color cyan, rec and chain A
color orange, rec and chain B

select keyres, rec and resi 13+17+20
show sticks, keyres
color magenta, keyres
label keyres and name CA, "%s%s" % (resn, resi)

show sticks, lig
color yellow, lig

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8

ray 1600, 1400
png {front}, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 45
ray 1600, 1400
png {front45}, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 90
ray 1600, 1400
png {side}, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 180
ray 1600, 1400
png {back}, dpi=300

quit
"""
    pml_path = pml_dir / f"{label}.pml"
    pml_path.write_text(pml.strip() + "\n", encoding="utf-8")
    return pml_path

for row in render_rows:
    make_pml(row)

print("PML files written to:", pml_dir)

for system in ["20chol", "40chol"]:
    sub = [r for r in best_rows if r["system"] == system]
    sub.sort(key=lambda r: float(r["composite_score"]), reverse=True)
    print(f"\nTop 20 best-pose-per-ligand for {system}:")
    for r in sub[:20]:
        print(
            r["group"],
            r["ligand_id"],
            "mode=", r["mode"],
            "affinity=", r["affinity_kcal_mol"],
            "minD=", r["min_distance_to_key_residues_A"],
            "contact4A=", r["n_contacted_key_residues_4A"],
            "near6A=", r["n_near_key_residues_6A"],
            "nearChains=", r["n_near_chains_6A"],
            "composite=", r["composite_score"],
        )
