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

load C:/TRKB_WP2/docking/receptors/20chol/TRKB_20chol_protein.pdb, rec
load C:/TRKB_WP2/docking/best_pose_review_pdb/reference_ISO_20chol_mode2.pdb, lig

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
png C:/TRKB_WP2/docking/figure_best/views/reference_ISO_20chol_mode2_front.png, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 45
ray 1600, 1400
png C:/TRKB_WP2/docking/figure_best/views/reference_ISO_20chol_mode2_front45.png, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 90
ray 1600, 1400
png C:/TRKB_WP2/docking/figure_best/views/reference_ISO_20chol_mode2_side.png, dpi=300

orient rec
zoom keyres or lig, 9
center keyres or lig
turn z, 10
turn x, -8
turn y, 180
ray 1600, 1400
png C:/TRKB_WP2/docking/figure_best/views/reference_ISO_20chol_mode2_back.png, dpi=300

quit
