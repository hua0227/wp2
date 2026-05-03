from pathlib import Path
import mdtraj as md

systems = {
    "20chol": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_protein.pdb"),
    "40chol": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_protein.pdb"),
}

for system_name, pdb_path in systems.items():
    print(f"\n=== {system_name} ===")
    traj = md.load(str(pdb_path))
    top = traj.topology

    chains = {}
    for res in top.residues:
        if not res.is_protein:
            continue
        chains.setdefault(res.chain.index, []).append(res)

    for chain_id, residues in chains.items():
        print(f"\nChain {chain_id}, residue count = {len(residues)}")
        print("Index-in-chain | resSeq | name | mdtraj_residue_index")
        for i, res in enumerate(residues):
            print(f"{i:3d} | {res.resSeq:5d} | {res.name:3s} | {res.index:5d}")

        print("\nCandidate Y433/V437/S440 pattern in this chain:")
        found = False

        # Y433 -> V437 is +4 residues
        # Y433 -> S440 is +7 residues
        for i, res in enumerate(residues):
            if i + 7 >= len(residues):
                continue

            r0 = residues[i]
            r4 = residues[i + 4]
            r7 = residues[i + 7]

            if r0.name == "TYR" and r4.name == "VAL" and r7.name == "SER":
                found = True
                print("FOUND candidate:")
                print(f"  Y433-like: chain={chain_id}, index_in_chain={i},   resSeq={r0.resSeq}, name={r0.name}, mdtraj_index={r0.index}")
                print(f"  V437-like: chain={chain_id}, index_in_chain={i+4}, resSeq={r4.resSeq}, name={r4.name}, mdtraj_index={r4.index}")
                print(f"  S440-like: chain={chain_id}, index_in_chain={i+7}, resSeq={r7.resSeq}, name={r7.name}, mdtraj_index={r7.index}")

        if not found:
            print("No Y +4 V +7 S pattern found in this chain.")