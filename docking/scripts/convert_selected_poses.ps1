$root = "C:\TRKB_WP2\docking"
$outDir = "$root\top_pose_review"

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$items = @(
    # Candidate 20% CHOL
    @{group="candidate"; id="6.1";  system="20chol"},
    @{group="candidate"; id="2.3";  system="20chol"},
    @{group="candidate"; id="18.1"; system="20chol"},
    @{group="candidate"; id="13.4"; system="20chol"},
    @{group="candidate"; id="11.4"; system="20chol"},
    @{group="candidate"; id="11.1"; system="20chol"},

    # Candidate 40% CHOL
    @{group="candidate"; id="8.1";  system="40chol"},
    @{group="candidate"; id="3.2";  system="40chol"},
    @{group="candidate"; id="20.2"; system="40chol"},
    @{group="candidate"; id="10.4"; system="40chol"},
    @{group="candidate"; id="13.4"; system="40chol"},
    @{group="candidate"; id="2.2";  system="40chol"},
    @{group="candidate"; id="19.2"; system="40chol"},
    @{group="candidate"; id="6.2";  system="40chol"},

    # Reference 20% CHOL
    @{group="reference"; id="FLX_R";  system="20chol"},
    @{group="reference"; id="FLX_S";  system="20chol"},
    @{group="reference"; id="HNK_RR"; system="20chol"},
    @{group="reference"; id="HNK_SS"; system="20chol"},
    @{group="reference"; id="DPH";    system="20chol"},

    # Reference 40% CHOL
    @{group="reference"; id="FLX_R";  system="40chol"},
    @{group="reference"; id="FLX_S";  system="40chol"},
    @{group="reference"; id="IMI";    system="40chol"},
    @{group="reference"; id="HNK_RR"; system="40chol"},
    @{group="reference"; id="HNK_SS"; system="40chol"},
    @{group="reference"; id="CPZ";    system="40chol"},
    @{group="reference"; id="ISO";    system="40chol"}
)

foreach ($item in $items) {
    $group = $item.group
    $id = $item.id
    $system = $item.system

    if ($group -eq "candidate") {
        $inFile = "$root\results_batch\$system\$id`_$system`_out.pdbqt"
    } else {
        $inFile = "$root\reference\results\$system\$id`_$system`_out.pdbqt"
    }

    $outFile = "$outDir\$group`_$id`_$system`_pose.pdb"

    if (Test-Path $inFile) {
        Write-Host "Converting:" $inFile
        obabel $inFile -O $outFile
    } else {
        Write-Host "MISSING:" $inFile
    }
}

Write-Host "Done. Output folder:"
Write-Host $outDir