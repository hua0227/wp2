# Pilot VMD/psfgen build script for ligand 14.2 / L002 / 20chol.
# This script does not run minimization or MD. It only attempts to append a
# standalone ligand segment to the existing apo PSF/PDB and write a new PSF/PDB.
#
# Notes:
# - Some psfgen/VMD builds do not fully support CHARMM stream files that contain
#   both topology and parameter sections. If topology commands for .str files
#   fail, run this manually in a VMD/psfgen environment that supports stream
#   topology sections, or split topology handling outside this pilot.
# - The ligand coordinate PDB below is a pilot-local derivative with atom names
#   copied from L002_14.2.str. The original reconstructed ligand PDB is unchanged.

package require psfgen
resetpsf

set toppar_dir [file normalize "C:/TRKB_WP2/TRKB_20CHOL/toppar"]
set ligand_str [file normalize "C:/TRKB_WP2/ligand_bound_MD/cgenff_parameterization/returned_str/L002_14.2.str"]
set apo_psf [file normalize "C:/TRKB_WP2/TRKB_20CHOL/step5_assembly.psf"]
set apo_shortmd_pdb [file normalize "C:/TRKB_WP2/TRKB_20CHOL/openmm_short_md_output/short_md_final.pdb"]
set apo_shortmd_psfgen_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_shortmd_psfgen_ready.pdb"]
set apo_charmmgui_pdb [file normalize "C:/TRKB_WP2/TRKB_20CHOL/step5_assembly.pdb"]
set ligand_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/L002_14.2_20chol_psfgen_ready.pdb"]
set out_psf [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.psf"]
set out_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_20chol/outputs/TRKB_20chol_L002_14.2_bound.pdb"]

topology [file join $toppar_dir "top_all36_prot.rtf"]
topology [file join $toppar_dir "top_all36_lipid.rtf"]
topology [file join $toppar_dir "top_all36_cgenff.rtf"]
topology [file join $toppar_dir "toppar_water_ions.str"]
topology [file join $toppar_dir "toppar_all36_lipid_cholesterol.str"]
topology $ligand_str

# Prefer the post-short-MD coordinates. The readpsf ... pdb form is used so the
# existing apo PSF can be paired with a coordinate PDB without regenerating apo.
# short_md_final.pdb itself uses PDB chain/resid naming that does not match the
# CHARMM-GUI PSF, so a pilot-local psfgen-ready PDB carries the short-MD xyz
# coordinates with step5_assembly.pdb atom identity/segid/resid columns.
if {[catch {readpsf $apo_psf pdb $apo_shortmd_psfgen_pdb} read_msg]} {
    puts "WARNING: readpsf with psfgen-ready short-MD coordinate PDB failed: $read_msg"
    puts "Trying CHARMM-GUI step5_assembly.pdb as coordinate fallback."
    resetpsf
    topology [file join $toppar_dir "top_all36_prot.rtf"]
    topology [file join $toppar_dir "top_all36_lipid.rtf"]
    topology [file join $toppar_dir "top_all36_cgenff.rtf"]
    topology [file join $toppar_dir "toppar_water_ions.str"]
    topology [file join $toppar_dir "toppar_all36_lipid_cholesterol.str"]
    topology $ligand_str
    readpsf $apo_psf pdb $apo_charmmgui_pdb
}

segment LIG {
    first none
    last none
    residue 1 L002
}
coordpdb $ligand_pdb LIG

guesscoord
writepsf $out_psf
writepdb $out_pdb

puts "Wrote ligand-bound PSF: $out_psf"
puts "Wrote ligand-bound PDB: $out_pdb"
quit
