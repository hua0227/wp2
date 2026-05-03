from pathlib import Path
import csv
import mdtraj as md
import numpy as np

systems = {
    "20chol": {
        "full": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_full.pdb"),
        "protein": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_protein.pdb"),
        "box": Path(r"C:\TRKB_WP2\docking\configs\TRKB_20chol_box.txt"),
        "report": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_residue_report.csv"),
    },
    "40chol": {
        "full": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_full.pdb"),
        "protein": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_protein.pdb"),
        "box": Path(r"C:\TRKB_WP2\docking\configs\TRKB_40chol_box.txt"),
        "report": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_residue_report.csv"),
    },
}

target_resseq = {433, 437, 440}
target_names = {"TYR", "VAL", "SER"}

for name, paths in systems.items():
    print(f"\n=== Processing {name} ===")
    traj = md.load(str(paths["full"]))

    protein_atoms = traj.topology.select("protein")
    if len(protein_atoms) == 0:
        raise RuntimeError(f"No protein atoms found in {paths['full']}")

    protein = traj.atom_slice(protein_atoms)
    protein.save_pdb(str(paths["protein"]))
    print("Saved protein receptor:", paths["protein"])
    print("Protein atoms:", len(protein_atoms))

    # Save residue report for manual checking
    with open(paths["report"], "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chain_index", "resSeq", "residue_name", "residue_index", "n_atoms"])
        for res in protein.topology.residues:
            if res.is_protein:
                writer.writerow([
                    res.chain.index,
                    res.resSeq,
                    res.name,
                    res.index,
                    len(list(res.atoms)),
                ])
    print("Saved residue report:", paths["report"])

    target_atoms = []
    target_residues_found = []

    for res in protein.topology.residues:
        if res.is_protein and res.resSeq in target_resseq:
            for atom in res.atoms:
                target_atoms.append(atom.index)
            target_residues_found.append((res.chain.index, res.resSeq, res.name))

    if len(target_atoms) == 0:
        print("WARNING: Could not find residues by resSeq 433/437/440.")
        print("Please inspect residue report manually:", paths["report"])
        continue

    xyz_nm = protein.xyz[0, target_atoms, :]
    center_nm = xyz_nm.mean(axis=0)
    center_angstrom = center_nm * 10.0

    print("Target residues found:")
    for item in target_residues_found:
        print("  chain_index, resSeq, name =", item)

    # Docking box size: enough for TMD dimer crevice; can later refine after visualization
    size_x = 22.0
    size_y = 22.0
    size_z = 24.0

    with open(paths["box"], "w") as f:
        f.write(f"center_x = {center_angstrom[0]:.3f}\n")
        f.write(f"center_y = {center_angstrom[1]:.3f}\n")
        f.write(f"center_z = {center_angstrom[2]:.3f}\n")
        f.write(f"size_x = {size_x:.1f}\n")
        f.write(f"size_y = {size_y:.1f}\n")
        f.write(f"size_z = {size_z:.1f}\n")

    print("Saved docking box:", paths["box"])
    print(f"Center Angstrom: {center_angstrom[0]:.3f}, {center_angstrom[1]:.3f}, {center_angstrom[2]:.3f}")