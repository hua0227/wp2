from pathlib import Path

root = Path(r"C:\TRKB_WP2\docking\pymol_packages")

pml_files = list(root.rglob("*.pml"))

print("PML files found:", len(pml_files))

for pml in pml_files:
    text = pml.read_text(encoding="utf-8", errors="ignore")

    fixed_lines = []
    for line in text.splitlines():
        # Remove non-ASCII characters to avoid PyMOL GBK decoding errors
        clean = line.encode("ascii", errors="ignore").decode("ascii")
        fixed_lines.append(clean)

    pml.write_text("\n".join(fixed_lines) + "\n", encoding="ascii")

print("Done. All PML files were rewritten as ASCII.")