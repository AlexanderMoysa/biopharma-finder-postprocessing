# Publish v0.1.0 archive to the existing repository.
# Requires Git for Windows and authenticated GitHub git access.
param(
  [string]$Archive = "$env:USERPROFILE\Downloads\biopharma-finder-postprocessing-v0.1.0.zip",
  [string]$Repo = "https://github.com/AlexanderMoysa/biopharma-finder-postprocessing.git"
)
$ErrorActionPreference = 'Stop'
if (-not (Get-Command git -ErrorAction SilentlyContinue)) { throw 'Install Git for Windows: https://git-scm.com/download/win' }
if (-not (Test-Path -LiteralPath $Archive -PathType Leaf)) { throw "Archive not found: $Archive" }
$work = Join-Path $env:TEMP ("bpf_publish_" + [guid]::NewGuid().ToString("N"))
$unpacked = Join-Path $work 'unpacked'
$checkout = Join-Path $work 'checkout'
New-Item -ItemType Directory -Path $unpacked -Force | Out-Null
try {
  Expand-Archive -LiteralPath $Archive -DestinationPath $unpacked -Force
  $source = Join-Path $unpacked "biopharma-finder-postprocessing"
  if (-not (Test-Path -LiteralPath (Join-Path $source "pyproject.toml"))) { throw "Unexpected archive structure" }
  & git clone $Repo $checkout
  if ($LASTEXITCODE -ne 0) { throw "Clone failed" }
  foreach ($item in Get-ChildItem -LiteralPath $source -Force) {
    Copy-Item -LiteralPath $item.FullName -Destination $checkout -Recurse -Force
  }
  & git -C $checkout add -A
  if ($LASTEXITCODE -ne 0) { throw "git add failed" }
  & git -C $checkout diff --cached --quiet
  if ($LASTEXITCODE -eq 0) { Write-Host "Nothing to commit"; return }
  & git -C $checkout commit -m "Import modular BioPharma Finder post-processing v0.1.0"
  if ($LASTEXITCODE -ne 0) { throw "Configure git author identity with git config user.name and user.email" }
  & git -C $checkout push origin HEAD:main
  if ($LASTEXITCODE -ne 0) { throw "Push failed: check GitHub authentication" }
  Write-Host "Published code to https://github.com/AlexanderMoysa/biopharma-finder-postprocessing"
} finally {
  if (Test-Path -LiteralPath $work) { Remove-Item -LiteralPath $work -Force -Recurse }
}
