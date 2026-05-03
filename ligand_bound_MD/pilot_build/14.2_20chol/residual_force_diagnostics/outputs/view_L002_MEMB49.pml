reinitialize
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/branchB_relaxation/outputs/branchB_stage4_minimized.pdb, branchB_stage4

hide everything, all
show cartoon, polymer.protein
color gray70, polymer.protein

select ligand_L002, (segi LIG or resn L002)
show sticks, ligand_L002
color yellow, ligand_L002

select MEMB49_POPC, (segi MEMB and resi 49 and resn POPC)
show sticks, MEMB49_POPC
color red, MEMB49_POPC

select near_MEMB49_5A, byres ((all within 5.0 of MEMB49_POPC) and not MEMB49_POPC and not ligand_L002 and not solvent)
show sticks, near_MEMB49_5A
color marine, near_MEMB49_5A

select key_residues, ((segi PROA or segi PROB) and ((resi 434 and resn TYR) or (resi 438 and resn VAL) or (resi 441 and resn SER)))
show sticks, key_residues
color hotpink, key_residues

show spheres, inorganic
set stick_radius, 0.18
set sphere_scale, 0.28
set cartoon_transparency, 0.35, polymer.protein
set valence, on

distance L002_to_MEMB49_min, ligand_L002, MEMB49_POPC, mode=2
hide labels, L002_to_MEMB49_min

zoom (ligand_L002 or MEMB49_POPC or near_MEMB49_5A or key_residues), 8
orient (ligand_L002 or MEMB49_POPC)
