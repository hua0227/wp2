$pymolExe = "D:\software\Pymol\Scripts\pymol.exe"

$root = "C:\TRKB_WP2\docking"
$pmlDir = "$root\figure_best\pml"
$viewDir = "$root\figure_best\views"

if (-not (Test-Path $pymolExe)) {
    throw "找不到 PyMOL，请修改脚本第一行 `$pymolExe 为你的 PyMOLWin.exe 路径。当前路径：$pymolExe"
}

New-Item -ItemType Directory -Force -Path $viewDir | Out-Null

$pmlFiles = Get-ChildItem $pmlDir -Filter *.pml | Sort-Object Name

Write-Host "Using PyMOL:" $pymolExe
Write-Host "PML files:" $pmlFiles.Count

foreach ($pml in $pmlFiles) {
    Write-Host "Rendering:" $pml.Name
    & $pymolExe -cq $pml.FullName
}

Write-Host "Rendering done."
Write-Host "View images should be in:"
Write-Host $viewDir