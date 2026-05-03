import csv
from pathlib import Path
import pymol2

ROOT = Path(r"C:\TRKB_WP2\docking")
JOB_CSV = ROOT / "figure_batch" / "configs" / "render_jobs.csv"
OUT_DIR = ROOT / "figure_batch" / "views"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# 你需要把这里替换成你在 PyMOL 中 get_view 得到的值
VIEW_20 = (
 -0.890625060,    0.332770050,    0.309920400,\
     0.296211958,   -0.092574105,    0.950625658,\
     0.345032126,    0.938448191,   -0.016118638,\
     0.000000000,    0.000000000,  -82.703971863,\
    -3.484630585,   -0.545082092,   -9.712956429,\
    54.398319244,  111.009651184,  -20.000000000
)

VIEW_40 = (
   -0.828766704,    0.157212898,    0.537048161,\
     0.534601450,   -0.061135426,    0.842882991,\
     0.165347666,    0.985670209,   -0.033383202,\
     0.000000000,    0.000000000, -106.664848328,\
    -6.618614197,   -3.383709908,  -10.366297722,\
    77.076637268,  136.252883911,  -20.000000000
)

def apply_style(cmd, receptor_obj, ligand_obj):
    cmd.bg_color("white")
    cmd.hide("everything", "all")

    cmd.show("cartoon", receptor_obj)
    cmd.color("cyan", receptor_obj)

    # 两条链分色（如果链信息正确）
    cmd.color("cyan", f"{receptor_obj} and chain A")
    cmd.color("lightorange", f"{receptor_obj} and chain B")

    # 关键残基
    cmd.select("keyres", f"{receptor_obj} and resi 13+17+20")
    cmd.show("sticks", "keyres")
    cmd.color("deeppink", "keyres")

    # ligand
    cmd.show("sticks", ligand_obj)
    cmd.color("yellow", ligand_obj)

    # 外观优化
    cmd.set("cartoon_transparency", 0.15)
    cmd.set("stick_radius", 0.18)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("specular", 0.25)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)

def render_views(cmd, label, receptor_obj, system):
    if system == "20chol":
        cmd.set_view(VIEW_20)
    else:
        cmd.set_view(VIEW_40)

    cmd.zoom("keyres or ligand", 10)

    views = [
        ("front", []),
        ("front45", [("turn", "y", 45)]),
        ("side", [("turn", "y", 90)]),
        ("back", [("turn", "y", 180)]),
    ]

    for view_name, ops in views:
        if system == "20chol":
            cmd.set_view(VIEW_20)
        else:
            cmd.set_view(VIEW_40)

        cmd.zoom("keyres or ligand", 10)

        for op in ops:
            if op[0] == "turn":
                _, axis, ang = op
                cmd.turn(axis, ang)

        out_png = OUT_DIR / f"{label}_{view_name}.png"
        cmd.ray(1600, 1400)
        cmd.png(str(out_png), dpi=300)

def main():
    with open(JOB_CSV, newline="", encoding="utf-8-sig") as f:
        jobs = list(csv.DictReader(f))

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.feedback("disable", "all", "everything")

        for row in jobs:
            label = row["label"]
            system = row["system"]
            receptor_pdb = row["receptor_pdb"]
            ligand_pdb = row["ligand_pdb"]

            print(f"Rendering {label} ...")

            cmd.reinitialize()

            cmd.load(receptor_pdb, "rec")
            cmd.load(ligand_pdb, "ligand")

            apply_style(cmd, "rec", "ligand", )
            render_views(cmd, label, "rec", system)

    print("Done.")

if __name__ == "__main__":
    main()