# Pilot-local Branch B cleanup for ligand 14.2 / L002 / 20chol.
# Deletes only the whole POPC molecule MEMB:94 from the branch-local copy.
# Original apo files, bound inputs, and ligand stream are not modified.

set input_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.psf"
set input_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.pdb"
set output_dir "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/targeted_clash_cleanup/outputs/branch_B_remove_POPC94"
set output_psf "$output_dir/branch_B_cleaned.psf"
set output_pdb "$output_dir/branch_B_cleaned.pdb"

file mkdir $output_dir
mol new $input_psf type psf waitfor all
mol addfile $input_pdb type pdb waitfor all

set del [atomselect top "segid MEMB and resid 94 and resname POPC"]
set del_count [$del num]
puts "Selected MEMB:94 POPC atoms for deletion: $del_count"
if {$del_count <= 0} {
    puts "ERROR: no atoms selected for MEMB:94 POPC; refusing to write branch_B cleaned system."
    quit 2
}

set forbidden [atomselect top "(protein or segid LIG or resname CHL or resname CHOL or resname CLR) and (segid MEMB and resid 94 and resname POPC)"]
if {[$forbidden num] > 0} {
    puts "ERROR: deletion selection overlaps forbidden protein/ligand/cholesterol atoms."
    quit 3
}

puts "Deleted residue atom listing:"
puts [$del get {index serial name resname segid resid}]

set keep [atomselect top "not (segid MEMB and resid 94 and resname POPC)"]
puts "Atoms retained: [$keep num]"
$keep writepsf $output_psf
$keep writepdb $output_pdb

$forbidden delete
$del delete
$keep delete
puts "Wrote Branch B cleaned PSF: $output_psf"
puts "Wrote Branch B cleaned PDB: $output_pdb"
quit
