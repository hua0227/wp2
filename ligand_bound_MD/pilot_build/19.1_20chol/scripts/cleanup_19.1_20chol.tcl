# Cleanup Tcl for 19.1 / L003 / 20chol
set input_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_20chol/outputs/TRKB_20chol_L003_19.1_bound.psf"
set input_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_20chol/outputs/TRKB_20chol_L003_19.1_bound.pdb"
set output_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_20chol/outputs/TRKB_20chol_L003_19.1_bound_cleaned.psf"
set output_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_20chol/outputs/TRKB_20chol_L003_19.1_bound_cleaned.pdb"

mol new $input_psf type psf waitfor all
mol addfile $input_pdb type pdb waitfor all
set keep [atomselect top {not ((segid "MEMB" and resid "94" and resname "POPC"))}]
$keep writepsf $output_psf
$keep writepdb $output_pdb
$keep delete
quit
