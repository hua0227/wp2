# Pilot-local cleaned system builder for 14.2 / L002 / 40chol.
# Deletes only dominant severe-clashing whole water/ion/lipid/cholesterol molecules selected by classify_40chol_ligand_clashes.py.
# Does not modify original apo files or ligand .str.
set input_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound.psf"
set input_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound.pdb"
set output_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound_cleaned.psf"
set output_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound_cleaned.pdb"

mol new $input_psf type psf waitfor all
mol addfile $input_pdb type pdb waitfor all
set del [atomselect top {(segid "MEMB" and resid "113" and resname "POPC")}]
puts "Deleting selected whole molecules/residues:"
puts [$del get {segid resid resname name index}]
set keep [atomselect top {not ((segid "MEMB" and resid "113" and resname "POPC"))}]
$keep writepsf $output_psf
$keep writepdb $output_pdb
$del delete
$keep delete
puts "Wrote cleaned PSF: $output_psf"
puts "Wrote cleaned PDB: $output_pdb"
quit
