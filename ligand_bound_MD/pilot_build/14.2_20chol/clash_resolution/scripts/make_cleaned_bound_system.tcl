# Pilot-local cleaned system builder for 14.2 / L002 / 20chol.
# Deletes only severe-clashing whole water/ion residues when present.
# Does not modify original apo files or ligand .str.
set input_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.psf"
set input_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.pdb"
set output_psf "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/clash_resolution/outputs/TRKB_20chol_L002_14.2_bound_cleaned.psf"
set output_pdb "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/clash_resolution/outputs/TRKB_20chol_L002_14.2_bound_cleaned.pdb"

puts "No severe water/ion clashes found. Copying bound PSF/PDB to cleaned outputs."
file copy -force $input_psf $output_psf
file copy -force $input_pdb $output_pdb
puts "Wrote cleaned PSF: $output_psf"
puts "Wrote cleaned PDB: $output_pdb"
quit
