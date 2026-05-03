$root = "C:\TRKB_WP2\docking"

$ligandDir = "$root\reference\pdbqt"

$rec20 = "$root\receptors\20chol\TRKB_20chol_receptor.pdbqt"
$rec40 = "$root\receptors\40chol\TRKB_40chol_receptor.pdbqt"

$box20File = "$root\configs\TRKB_20chol_box.txt"
$box40File = "$root\configs\TRKB_40chol_box.txt"

$out20 = "$root\reference\results\20chol"
$out40 = "$root\reference\results\40chol"
$logDir = "$root\reference\logs"

New-Item -ItemType Directory -Force -Path $out20 | Out-Null
New-Item -ItemType Directory -Force -Path $out40 | Out-Null
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$condaEnv = "C:\Users\14566\miniconda3\envs\trkb_openmm"
$vinaCandidates = @(
    "$condaEnv\Library\bin\vina.exe",
    "$condaEnv\Scripts\vina.exe",
    "$condaEnv\bin\vina.exe"
)

$vinaExe = $null
foreach ($candidate in $vinaCandidates) {
    if (Test-Path $candidate) {
        $vinaExe = $candidate
        break
    }
}

if (-not $vinaExe) {
    throw "Could not find vina.exe in trkb_openmm."
}

Write-Host "Using vina:" $vinaExe

function Read-BoxFile {
    param([string]$boxFile)

    $box = @{}
    Get-Content $boxFile | ForEach-Object {
        if ($_ -match "^\s*(\w+)\s*=\s*([-0-9.]+)") {
            $box[$matches[1]] = $matches[2]
        }
    }
    return $box
}

function Run-Vina {
    param(
        [string]$receptor,
        [string]$ligand,
        [hashtable]$box,
        [string]$outFile,
        [string]$logFile
    )

    $vinaArgs = @(
        "--receptor", $receptor,
        "--ligand", $ligand,
        "--center_x", $box["center_x"],
        "--center_y", $box["center_y"],
        "--center_z", $box["center_z"],
        "--size_x", $box["size_x"],
        "--size_y", $box["size_y"],
        "--size_z", $box["size_z"],
        "--exhaustiveness", "8",
        "--num_modes", "10",
        "--energy_range", "5",
        "--out", $outFile
    )

    & $script:vinaExe @vinaArgs *> $logFile
    return $LASTEXITCODE
}

$box20 = Read-BoxFile $box20File
$box40 = Read-BoxFile $box40File

$ligands = Get-ChildItem $ligandDir -Filter *.pdbqt | Sort-Object Name

Write-Host "Found reference ligands:" $ligands.Count

foreach ($lig in $ligands) {
    $id = [System.IO.Path]::GetFileNameWithoutExtension($lig.Name)

    Write-Host "Docking reference 20chol:" $id
    $exit20 = Run-Vina `
        -receptor $rec20 `
        -ligand $lig.FullName `
        -box $box20 `
        -outFile "$out20\$id`_20chol_out.pdbqt" `
        -logFile "$logDir\$id`_20chol.log"

    if ($exit20 -ne 0) {
        Write-Host "WARNING: 20chol failed for $id, exit code $exit20"
    }

    Write-Host "Docking reference 40chol:" $id
    $exit40 = Run-Vina `
        -receptor $rec40 `
        -ligand $lig.FullName `
        -box $box40 `
        -outFile "$out40\$id`_40chol_out.pdbqt" `
        -logFile "$logDir\$id`_40chol.log"

    if ($exit40 -ne 0) {
        Write-Host "WARNING: 40chol failed for $id, exit code $exit40"
    }
}

Write-Host "Reference docking finished."