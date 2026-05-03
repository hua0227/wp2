reinitialize
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/restraint_strategy_fix/outputs/R1_pbcaware_minimized.pdb, R1_pbcaware
load C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/restraint_strategy_fix/outputs/R2_pbcaware_minimized.pdb, R2_pbcaware

hide everything, all
show cartoon, R1_pbcaware and polymer.protein
show cartoon, R2_pbcaware and polymer.protein
color gray80, R1_pbcaware and polymer.protein
color gray45, R2_pbcaware and polymer.protein

select R1_L002, R1_pbcaware and (segi LIG or resn L002)
select R2_L002, R2_pbcaware and (segi LIG or resn L002)
show sticks, R1_L002 or R2_L002
color yellow, R1_L002 or R2_L002

select R1_PROA437_438, R1_pbcaware and segi PROA and resi 437+438
select R2_PROA437_438, R2_pbcaware and segi PROA and resi 437+438
show sticks, R1_PROA437_438 or R2_PROA437_438
color red, R1_PROA437_438 or R2_PROA437_438

select key_residues, (segi PROA or segi PROB) and ((resi 434 and resn TYR) or (resi 438 and resn VAL) or (resi 441 and resn SER))
show sticks, key_residues
color hotpink, key_residues

select near_PROA437_5A, byres ((all within 5.0 of (R2_pbcaware and segi PROA and resi 437)) and not solvent)
show sticks, near_PROA437_5A
color marine, near_PROA437_5A and not (R2_L002 or R2_PROA437_438 or key_residues)

set stick_radius, 0.18
set cartoon_transparency, 0.45
zoom (R2_L002 or R2_PROA437_438 or near_PROA437_5A), 8
orient (R2_L002 or R2_PROA437_438)
