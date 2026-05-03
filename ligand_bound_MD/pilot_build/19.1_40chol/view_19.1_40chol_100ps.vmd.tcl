# ============================================================
# VMD visualization script
# System: TRKB 40chol + ligand 19.1 / L003
# Trajectory: 100 ps R2 restrained MD pilot
# Purpose:
#   1. Load PSF + DCD
#   2. Fix apparent PBC wrapping artifacts
#   3. Show protein, ligand, key residues, nearby membrane
#   4. Export 0/25/50/75/100 ps snapshots
# ============================================================

# ----------------------------
# User settings
# ----------------------------
set base "C:/TRKB_WP2/ligand_bound_MD/pilot_build/19.1_40chol"

set psf_file "$base/outputs/TRKB_40chol_L003_19.1_bound_cleaned.psf"
set dcd_file "$base/md100ps_R2/outputs/TRKB_40chol_L003_19.1_100ps_R2.dcd"
set gro_file "C:/TRKB_WP2/TRKB_40CHOL/gromacs/step5_input.gro"

set outdir "$base/vmd_visualization_100ps"
set ligand_sel "resname L003 or segid LIG"
set keyres_sel "protein and (resid 13 or resid 17 or resid 20)"

# Create output directory if needed
file mkdir $outdir

# ----------------------------
# Basic display settings
# ----------------------------
display resetview
display projection Orthographic
display depthcue off
display shadows on
display ambientocclusion on
axes location Off
color Display Background white

# ----------------------------
# Load trajectory
# ----------------------------
puts "Loading PSF:"
puts $psf_file
mol new $psf_file type psf waitfor all

puts "Loading DCD:"
puts $dcd_file
mol addfile $dcd_file type dcd waitfor all

set molid [molinfo top]
set nframes [molinfo $molid get numframes]
puts "Loaded frames: $nframes"

# ----------------------------
# Read box size from GRO file
# GRO box is in nm; VMD uses Angstrom
# ----------------------------
set fh [open $gro_file r]
set gro_lines [split [read $fh] "\n"]
close $fh

# Last non-empty line
set lastline ""
foreach line $gro_lines {
    if {[string trim $line] ne ""} {
        set lastline [string trim $line]
    }
}

set fields [regexp -all -inline {\S+} $lastline]
set lx_nm [lindex $fields 0]
set ly_nm [lindex $fields 1]
set lz_nm [lindex $fields 2]

set lx_A [expr {$lx_nm * 10.0}]
set ly_A [expr {$ly_nm * 10.0}]
set lz_A [expr {$lz_nm * 10.0}]

puts "GRO box nm: $lx_nm $ly_nm $lz_nm"
puts "VMD box A:  $lx_A $ly_A $lz_A"

# ----------------------------
# PBC repair
# This is the important part for the problem:
# "two chains separated / protein not in membrane"
# ----------------------------
package require pbctools

pbc set [list $lx_A $ly_A $lz_A 90 90 90] -all

# Unwrap trajectory, join broken fragments, then wrap around protein+ligand
puts "Applying PBC unwrap/join/wrap..."
pbc unwrap -all
pbc join fragment -all
pbc wrap -all -compound fragment -center com -centersel "protein or $ligand_sel" -all

# Draw box for reference
pbc box -color black -width 2

# ----------------------------
# Remove default representation
# ----------------------------
mol delrep 0 $molid

# ----------------------------
# Representation 1: protein cartoon
# ----------------------------
mol representation NewCartoon 0.35 10 4
mol color Chain
mol selection {protein}
mol material Opaque
mol addrep $molid

# ----------------------------
# Representation 2: ligand L003
# ----------------------------
mol representation Licorice 0.28 14 14
mol color ColorID 4
mol selection $ligand_sel
mol material Opaque
mol addrep $molid

# ----------------------------
# Representation 3: key residues TYR13 / VAL17 / SER20
# If this selection fails, check actual resid in VMD.
# ----------------------------
mol representation Licorice 0.28 14 14
mol color ColorID 1
mol selection $keyres_sel
mol material Opaque
mol addrep $molid

# ----------------------------
# Representation 4: protein residues within 5 A of ligand
# ----------------------------
mol representation Licorice 0.20 10 10
mol color ColorID 9
mol selection "same residue as (protein within 5 of ($ligand_sel))"
mol material Opaque
mol addrep $molid

# ----------------------------
# Representation 5: nearby membrane lipids/cholesterol around protein/ligand
# Keep this light to avoid blocking the ligand view.
# ----------------------------
mol representation Lines 1.0
mol color ColorID 8
mol selection "same residue as ((resname POPC CHL1 CHOL CLR) and within 12 of (protein or $ligand_sel))"
mol material Transparent
mol addrep $molid

# ----------------------------
# Representation 6: hide water/ions by not displaying them
# If you want water, add another rep manually.
# ----------------------------

# ----------------------------
# Set view around ligand and key residues
# ----------------------------
set center_sel [atomselect $molid "($ligand_sel) or ($keyres_sel)"]
set center [measure center $center_sel]
$center_sel delete

puts "Centering view around ligand and key residues..."
translate to $center
scale by 1.4
display update

# ----------------------------
# Helper: render snapshot at chosen frame
# ----------------------------
proc render_frame {molid frame label outdir ligand_sel keyres_sel} {
    animate goto $frame
    display update

    # Recenter at each frame
    set s [atomselect $molid "($ligand_sel) or ($keyres_sel)" frame $frame]
    set c [measure center $s]
    $s delete

    # Reset view lightly and center
    display update

    set outfile "$outdir/snapshot_${label}.tga"
    puts "Rendering frame $frame as $outfile"
    render TachyonInternal $outfile
}

# ----------------------------
# Export snapshots at 0, 25, 50, 75, 100 percent of trajectory
# If DCD has 50 frames or 100 frames, this still works.
# ----------------------------
set f0 0
set f25 [expr {int(($nframes - 1) * 0.25)}]
set f50 [expr {int(($nframes - 1) * 0.50)}]
set f75 [expr {int(($nframes - 1) * 0.75)}]
set f100 [expr {$nframes - 1}]

render_frame $molid $f0   "000ps" $outdir $ligand_sel $keyres_sel
render_frame $molid $f25  "025ps" $outdir $ligand_sel $keyres_sel
render_frame $molid $f50  "050ps" $outdir $ligand_sel $keyres_sel
render_frame $molid $f75  "075ps" $outdir $ligand_sel $keyres_sel
render_frame $molid $f100 "100ps" $outdir $ligand_sel $keyres_sel

# ----------------------------
# Save current VMD state
# ----------------------------
save_state "$outdir/view_19.1_40chol_100ps_state.vmd"

puts "============================================================"
puts "VMD visualization setup complete."
puts "Output directory:"
puts $outdir
puts ""
puts "Snapshots saved as .tga files:"
puts "snapshot_000ps.tga"
puts "snapshot_025ps.tga"
puts "snapshot_050ps.tga"
puts "snapshot_075ps.tga"
puts "snapshot_100ps.tga"
puts ""
puts "If the protein still looks split, run in Tk Console:"
puts "pbc unwrap -all"
puts "pbc join fragment -all"
puts "pbc wrap -all -compound fragment -center com -centersel \"protein or resname L003\" -all"
puts "============================================================"