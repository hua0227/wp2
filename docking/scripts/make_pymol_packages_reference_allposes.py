from pathlib import Path
import shutil
import subprocess
import csv

ROOT = Path(r"C:\TRKB_WP2\docking")
OUT_ROOT = ROOT / "pymol_packages_reference"
INDEX_CSV = OUT_ROOT / "pymol_reference_packages_index.csv"

DATASETS = [
    {
        "system": "20chol",
        "receptor_pdb": ROOT / "receptors" / "20chol" / "TRKB_20chol_protein.pdb",
        "pose_dir": ROOT / "reference" / "results" / "20chol",
        "suffix": "_20chol_out.pdbqt",
    },
    {
        "system": "40chol",
        "receptor_pdb": ROOT / "receptors" / "40chol" / "TRKB_40chol_protein.pdb",
        "pose_dir": ROOT / "reference" / "results" / "40chol",
        "suffix": "_40chol_out.pdbqt",
    },
]

POSE_COLORS = [
    "yellow",
    "cyan",
    "magenta",
    "salmon",
    "lime",
    "orange",
    "marine",
    "violet",
    "wheat",
    "tv_blue",
]

def find_obabel():
    candidates = [
        shutil.which("obabel"),
        r"D:\software\OpenBabel-3.1.1\obabel.exe",
        r"C:\Users\14566\miniconda3\envs\trkb_openmm\Library\bin\obabel.exe",
        r"C:\Users\14566\miniconda3\envs\trkb_openmm\Scripts\obabel.exe",
    ]
    for c in candidates:
        if c and Path(c).exists():
            return str(c)
    raise FileNotFoundError("Cannot find obabel.exe.")

def parse_vina_models(pdbqt_path: Path):
    lines = pdbqt_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)

    models = []
    current = []
    in_model = False

    for line in lines:
        if line.startswith("MODEL"):
            if current:
                if any(x.startswith(("ATOM", "HETATM")) for x in current):
                    models.append(current)
            current = [line]
            in_model = True
            continue

        if line.startswith("ENDMDL"):
            if in_model:
                current.append(line)
                if any(x.startswith(("ATOM", "HETATM")) for x in current):
                    models.append(current)
                current = []
                in_model = False
            continue

        if in_model:
            current.append(line)

    if not models:
        if any(x.startswith(("ATOM", "HETATM")) for x in lines):
            models = [lines]

    return models

def write_pml(pkg_dir: Path, label: str, nposes: int):
    pml_path = pkg_dir / f"open_{label}.pml"

    lines = []
    lines.append("reinitialize")
    lines.append("bg_color white")
    lines.append("set orthoscopic, on")
    lines.append("set depth_cue, off")
    lines.append("set ray_opaque_background, off")
    lines.append("set cartoon_fancy_helices, on")
    lines.append("set cartoon_transparency, 0.10")
    lines.append("set stick_radius, 0.18")
    lines.append("set label_size, 18")
    lines.append("set label_color, black")
    lines.append("set valence, 0")
    lines.append("")
    lines.append("load receptor.pdb, receptor")
    lines.append("hide everything, all")
    lines.append("show cartoon, receptor")
    lines.append("color cyan, receptor")
    lines.append("color cyan, receptor and chain A")
    lines.append("color lightorange, receptor and chain B")
    lines.append("")
    lines.append("select keyres, receptor and resi 13+17+20")
    lines.append("show sticks, keyres")
    lines.append("color hotpink, keyres")
    lines.append('label keyres and name CA, "%s%s" % (resn, resi)')
    lines.append("")

    pose_names = []
    for i in range(1, nposes + 1):
        obj = f"pose{i:02d}"
        pose_names.append(obj)
        color = POSE_COLORS[(i - 1) % len(POSE_COLORS)]
        lines.append(f"load pose{i:02d}.pdb, {obj}")
        lines.append(f"show sticks, {obj}")
        lines.append(f"color {color}, {obj}")

    lines.append("")
    if pose_names:
        lines.append("group poses, " + " ".join(pose_names))
        for obj in pose_names[1:]:
            lines.append(f"disable {obj}")
        lines.append("enable pose01")
        lines.append("zoom keyres or pose01, 10")
    else:
        lines.append("zoom keyres, 10")

    lines.append("")
    lines.append("# Commands:")
    lines.append("# disable pose01")
    lines.append("# enable pose02")
    lines.append("# enable pose*")
    lines.append("# disable pose*")
    lines.append("# enable pose07")
    lines.append("# zoom keyres or pose07, 10")
    lines.append("# show surface, receptor")
    lines.append("# set transparency, 0.65, receptor")

    pml_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return pml_path

def convert_with_obabel(obabel_exe: str, in_pdbqt: Path, out_pdb: Path):
    cmd = [obabel_exe, str(in_pdbqt), "-O", str(out_pdb)]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel conversion failed:\n{result.stdout}")

def main():
    obabel_exe = find_obabel()
    print("Using obabel:", obabel_exe)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    index_rows = []

    for ds in DATASETS:
        system = ds["system"]
        receptor_pdb = ds["receptor_pdb"]
        pose_dir = ds["pose_dir"]
        suffix = ds["suffix"]

        system_out = OUT_ROOT / system
        system_out.mkdir(parents=True, exist_ok=True)

        pose_files = sorted(pose_dir.glob(f"*{suffix}"))
        print(f"\n=== {system} ===")
        print("Reference pose files:", len(pose_files))

        for pose_file in pose_files:
            ligand_id = pose_file.name.replace(suffix, "")
            label = f"{ligand_id}-{system}"
            pkg_dir = system_out / label
            pkg_dir.mkdir(parents=True, exist_ok=True)

            receptor_out = pkg_dir / "receptor.pdb"
            shutil.copy2(receptor_pdb, receptor_out)

            models = parse_vina_models(pose_file)
            if not models:
                print("WARNING: no models found in", pose_file)
                continue

            for i, model_lines in enumerate(models, start=1):
                pose_pdbqt = pkg_dir / f"pose{i:02d}.pdbqt"
                pose_pdb = pkg_dir / f"pose{i:02d}.pdb"

                pose_pdbqt.write_text("".join(model_lines), encoding="utf-8")
                convert_with_obabel(obabel_exe, pose_pdbqt, pose_pdb)

            pml_path = write_pml(pkg_dir, label, len(models))

            index_rows.append({
                "system": system,
                "ligand_id": ligand_id,
                "label": label,
                "package_dir": str(pkg_dir),
                "pml_file": str(pml_path),
                "n_poses": len(models),
            })

            print("Done:", label, f"poses={len(models)}")

    with open(INDEX_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["system", "ligand_id", "label", "package_dir", "pml_file", "n_poses"]
        )
        writer.writeheader()
        writer.writerows(index_rows)

    print("\nAll reference packages done.")
    print("Index written to:", INDEX_CSV)

if __name__ == "__main__":
    main()