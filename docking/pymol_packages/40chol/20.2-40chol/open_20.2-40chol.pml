reinitialize
bg_color white
set orthoscopic, on
set depth_cue, off
set ray_opaque_background, off
set cartoon_fancy_helices, on
set cartoon_transparency, 0.10
set stick_radius, 0.18
set label_size, 18
set label_color, black
set valence, 0

load receptor.pdb, receptor
hide everything, all
show cartoon, receptor
color cyan, receptor
color cyan, receptor and chain A
color lightorange, receptor and chain B

select keyres, receptor and resi 13+17+20
show sticks, keyres
color hotpink, keyres
label keyres and name CA, "%s%s" % (resn, resi)

load pose01.pdb, pose01
show sticks, pose01
color yellow, pose01
load pose02.pdb, pose02
show sticks, pose02
color cyan, pose02
load pose03.pdb, pose03
show sticks, pose03
color magenta, pose03
load pose04.pdb, pose04
show sticks, pose04
color salmon, pose04
load pose05.pdb, pose05
show sticks, pose05
color lime, pose05
load pose06.pdb, pose06
show sticks, pose06
color orange, pose06
load pose07.pdb, pose07
show sticks, pose07
color marine, pose07
load pose08.pdb, pose08
show sticks, pose08
color violet, pose08
load pose09.pdb, pose09
show sticks, pose09
color wheat, pose09
load pose10.pdb, pose10
show sticks, pose10
color tv_blue, pose10

group poses, pose01 pose02 pose03 pose04 pose05 pose06 pose07 pose08 pose09 pose10
disable pose02
disable pose03
disable pose04
disable pose05
disable pose06
disable pose07
disable pose08
disable pose09
disable pose10
enable pose01
zoom keyres or pose01, 10

# 
# enable pose02
# disable pose01
# enable pose*
# show surface, receptor
# set transparency, 0.65, receptor
