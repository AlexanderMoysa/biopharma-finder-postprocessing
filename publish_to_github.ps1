<#
Optional publishing helper. Run from the extracted project root only after
reviewing the code. Requires Git and authenticated GitHub CLI (`gh auth login`).
Defaults to a PRIVATE repository under the signed-in user's GitHub account.
#>
param([string]$RepoName = "biopharma-finder-postprocessing")
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw "Git is not installed. Install Git for Windows first."
}
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw "GitHub CLI is not installed. Install it from https://cli.github.com/ and run: gh auth login"
}
& gh auth status
if ($LASTEXITCODE -ne 0) { throw "Log into GitHub CLI: gh auth login" }
if (-not (Test-Path ".git")) {
    & git init -b main
    if ($LASTEXITCODE -ne 0) { throw "git init failed" }
}
& git add -A
if ($LASTEXITCODE -ne 0) { throw "git add failed" }
& git diff --cached --quiet
if ($LASTEXITCODE -ne 0) {
    & git commit -m "Initial modular BioPharma Finder post-processing scaffold"
    if ($LASTEXITCODE -ne 0) { throw "git commit failed; configure git user.name and user.email" }
}
& gh repo create $RepoName --private --source . --remote origin --push
if ($LASTEXITCODE -ne 0) { throw "GitHub repository creation/push failed" }
Write-Host "Created private repository and pushed initial code."
