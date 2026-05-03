reinitialize
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/md100ps_R2/outputs/TRKB_20chol_L008_6.2_100ps_R2_final.pdb, TRKB_20chol_L008_6_2
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_40chol/md100ps_R2/outputs/TRKB_40chol_L008_6.2_100ps_R2_final.pdb, TRKB_40chol_L008_6_2
hide everything
show cartoon, polymer.protein
show sticks, resn L008
show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)
color marine, TRKB_20chol_L008_6_2
color firebrick, TRKB_40chol_L008_6_2
color cyan, TRKB_20chol_L008_6_2 and resn L008
color orange, TRKB_40chol_L008_6_2 and resn L008
color yellow, resn TYR+VAL+SER and resi 434+438+441
set stick_radius, 0.18
set cartoon_transparency, 0.35
zoom (resn L008) or (resi 434+438+441), 8
orient (resn L008) or (resi 434+438+441)
# Visualization only; not used for simulation.
