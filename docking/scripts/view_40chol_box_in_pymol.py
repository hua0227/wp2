from pymol import cmd
from pymol.cgo import LINEWIDTH, BEGIN, LINES, COLOR, VERTEX, END

# ===== user settings =====
pdb_path = r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_protein.pdb"
obj_name = "rec40"

center_x = -2.374
center_y = 0.341
center_z = -8.937
size_x = 22.0
size_y = 22.0
size_z = 24.0

# 已按 mapping 文件修正
key_resi = "13+17+20"

def draw_box(name, cx, cy, cz, sx, sy, sz, color=(1.0, 1.0, 0.0), width=3.0):
    hx = sx / 2.0
    hy = sy / 2.0
    hz = sz / 2.0

    corners = [
        [cx-hx, cy-hy, cz-hz],
        [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz],
        [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz],
        [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz],
        [cx-hx, cy+hy, cz+hz],
    ]

    edges = [
        (0,1),(1,2),(2,3),(3,0),
        (4,5),(5,6),(6,7),(7,4),
        (0,4),(1,5),(2,6),(3,7)
    ]

    obj = [LINEWIDTH, width, BEGIN, LINES, COLOR, color[0], color[1], color[2]]
    for i, j in edges:
        obj.extend([VERTEX, *corners[i], VERTEX, *corners[j]])
    obj.append(END)

    cmd.load_cgo(obj, name)

cmd.load(pdb_path, obj_name)
cmd.hide("everything", obj_name)
cmd.show("cartoon", obj_name)
cmd.color("cyan", f"{obj_name} and chain A")
cmd.color("orange", f"{obj_name} and chain B")

cmd.select("keyres40", f"{obj_name} and resi {key_resi}")
cmd.show("sticks", "keyres40")
cmd.color("hotpink", "keyres40")
cmd.label("keyres40 and name CA", '"%s%s" % (resn, resi)')

cmd.show("spheres", "keyres40 and name CA")
cmd.set("sphere_scale", 0.35, "keyres40 and name CA")

cmd.pseudoatom("box_center40", pos=[center_x, center_y, center_z], label="BOX_CENTER")
cmd.show("spheres", "box_center40")
cmd.color("yellow", "box_center40")
cmd.set("sphere_scale", 0.5, "box_center40")

draw_box("dockbox40", center_x, center_y, center_z, size_x, size_y, size_z)

cmd.zoom("keyres40 or dockbox40", 12)