# Run the EDF pipeline for every subject folder from a terminal (no MATLAB GUI).
#
#   .\run_all_subjects.ps1                               # all folders with EDFs
#   .\run_all_subjects.ps1 -Subjects 001-s,005-s         # only these
#   .\run_all_subjects.ps1 -OutputRoot D:\out -CasesFile cases.csv
#
# Output is also saved to <OutputRoot>\run_all_subjects_<timestamp>.log
param(
    [string]$DataRoot   = (Split-Path -Parent $PSScriptRoot),
    [string]$OutputRoot = '',
    [string[]]$Subjects = @(),
    [string]$CasesFile  = '',
    [string]$RecordingLog = '',
    [string]$Matlab     = 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe'
)

if (-not $RecordingLog) { $RecordingLog = Join-Path $DataRoot 'EEG_recording_log.xlsx' }
if (-not $OutputRoot) { $OutputRoot = Join-Path $DataRoot 'pipeline_output' }
New-Item -ItemType Directory -Force $OutputRoot | Out-Null
$log = Join-Path $OutputRoot ("run_all_subjects_{0}.log" -f (Get-Date -Format 'yyyyMMdd_HHmmss'))

function q([string]$s) { "'" + $s.Replace("'", "''") + "'" }
$Subjects = $Subjects | ForEach-Object { $_ -split ',' } | Where-Object { $_ }   # -File passes "a,b" as one string
$subj = '{' + (($Subjects | ForEach-Object { q $_ }) -join ',') + '}'
$repo = $PSScriptRoot
$cmd  = "cd($(q $repo)); run_all_subjects($(q $DataRoot), $(q $OutputRoot), $subj, $(q $CasesFile), $(q $RecordingLog))"

Write-Host "MATLAB: $cmd"
Write-Host "Log:    $log"
& $Matlab -batch $cmd 2>&1 | Tee-Object -FilePath $log
exit $LASTEXITCODE
