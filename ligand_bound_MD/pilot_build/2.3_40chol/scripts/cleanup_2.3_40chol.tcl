# Cleanup Tcl for 2.3 / L006 / 40chol
set input_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_40chol/outputs/TRKB_40chol_L006_2.3_bound.psf"
set input_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_40chol/outputs/TRKB_40chol_L006_2.3_bound.pdb"
set output_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_40chol/outputs/TRKB_40chol_L006_2.3_bound_cleaned.psf"
set output_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/2.3_40chol/outputs/TRKB_40chol_L006_2.3_bound_cleaned.pdb"

mol new $input_psf type psf waitfor all
mol addfile $input_pdb type pdb waitfor all
set keep [atomselect top {not ((segid "MEMB" and resid "113" and resname "POPC"))}]
$keep writepsf $output_psf
$keep writepdb $output_pdb
$keep delete
quit
