reinitialize
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/md100ps_R2/outputs/TRKB_20chol_L002_14.2_100ps_R2_final.pdb, 14_2_20chol_L002
color marine, 14_2_20chol_L002
show sticks, 14_2_20chol_L002 and resn L002
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/md100ps_R2/outputs/TRKB_40chol_L002_14.2_100ps_R2_final.pdb, 14_2_40chol_L002
color bluewhite, 14_2_40chol_L002
show sticks, 14_2_40chol_L002 and resn L002
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_20chol/md100ps_R2/outputs/TRKB_20chol_L003_19.1_100ps_R2_final.pdb, 19_1_20chol_L003
color forest, 19_1_20chol_L003
show sticks, 19_1_20chol_L003 and resn L003
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_40chol/md100ps_R2/outputs/TRKB_40chol_L003_19.1_100ps_R2_final.pdb, 19_1_40chol_L003
color limegreen, 19_1_40chol_L003
show sticks, 19_1_40chol_L003 and resn L003
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_20chol/md100ps_R2/outputs/TRKB_20chol_L006_2.3_100ps_R2_final.pdb, 2_3_20chol_L006
color violet, 2_3_20chol_L006
show sticks, 2_3_20chol_L006 and resn L006
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_40chol/md100ps_R2/outputs/TRKB_40chol_L006_2.3_100ps_R2_final.pdb, 2_3_40chol_L006
color magenta, 2_3_40chol_L006
show sticks, 2_3_40chol_L006 and resn L006
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/md100ps_R2/outputs/TRKB_20chol_L008_6.2_100ps_R2_final.pdb, 6_2_20chol_L008
color orange, 6_2_20chol_L008
show sticks, 6_2_20chol_L008 and resn L008
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_40chol/md100ps_R2/outputs/TRKB_40chol_L008_6.2_100ps_R2_final.pdb, 6_2_40chol_L008
color tv_orange, 6_2_40chol_L008
show sticks, 6_2_40chol_L008 and resn L008
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/17.2_20chol/md100ps_R2/outputs/TRKB_20chol_L010_17.2_100ps_R2_final.pdb, 17_2_20chol_L010
color teal, 17_2_20chol_L010
show sticks, 17_2_20chol_L010 and resn L010
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/17.2_40chol/md100ps_R2/outputs/TRKB_40chol_L010_17.2_100ps_R2_final.pdb, 17_2_40chol_L010
color cyan, 17_2_40chol_L010
show sticks, 17_2_40chol_L010 and resn L010
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/9.2_20chol/md100ps_R2/outputs/TRKB_20chol_L004_9.2_100ps_R2_final.pdb, 9_2_20chol_L004
color salmon, 9_2_20chol_L004
show sticks, 9_2_20chol_L004 and resn L004
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/9.2_40chol/md100ps_R2/outputs/TRKB_40chol_L004_9.2_100ps_R2_final.pdb, 9_2_40chol_L004
color red, 9_2_40chol_L004
show sticks, 9_2_40chol_L004 and resn L004
hide everything
show cartoon, polymer.protein
show sticks, resn L002+L003+L004+L006+L008+L010
show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)
color yellow, resn TYR+VAL+SER and resi 434+438+441
set cartoon_transparency, 0.45
set stick_radius, 0.18
zoom (resn L002+L003+L004+L006+L008+L010) or (resi 434+438+441), 10
orient (resn L002+L003+L004+L006+L008+L010) or (resi 434+438+441)
# Visualization only; not used for simulation.
