# Run the EDF campaign (src/campaign_config.m) from a terminal (no MATLAB GUI).
#
#   .\run_all_subjects.ps1                                   # every animal folder in the Excel
#   .\run_all_subjects.ps1 -FilesCsv edf_list.csv            # only the EDFs listed (animal_id, edf_file)
#   .\run_all_subjects.ps1 -Subjects 004-s,005-s             # only these animals
#   .\run_all_subjects.ps1 -OutputRoot D:\out -FilesCsv edf_list.csv
#
# Default DataRoot: edfDir of run_all_subjects.m (the EDF folder).
# Default OutputRoot: <parent folder of this repo>\pipeline_output_campaign_<yyyyMMdd>
# (pipeline_output_campaign_preliminar_<yyyyMMdd> with -FilesCsv).
# Output is also saved to <OutputRoot>\run_all_subjects_<timestamp>.log
param(
    [string]$DataRoot   = '',   # '' -> edfDir of run_all_subjects.m
    [string]$OutputRoot = '',
    [string[]]$Subjects = @(),
    [string]$CasesFile  = '',
    [string]$RecordingLog = '',
    [string]$FilesCsv   = '',
    [string]$Matlab     = 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe'
)

if (-not $RecordingLog -and $DataRoot) { $RecordingLog = Join-Path $DataRoot 'EEG_recording_log.xlsx' }   # '' -> <edfDir>\EEG_recording_log.xlsx (run_all_subjects.m)
if ($FilesCsv) { $FilesCsv = (Resolve-Path $FilesCsv).Path }
if (-not $OutputRoot) {
    $stamp = Get-Date -Format 'yyyyMMdd'
    if ($FilesCsv) { $OutputRoot = Join-Path (Split-Path -Parent $PSScriptRoot) "pipeline_output_campaign_preliminar_$stamp" }
    else           { $OutputRoot = Join-Path (Split-Path -Parent $PSScriptRoot) "pipeline_output_campaign_$stamp" }
}
New-Item -ItemType Directory -Force $OutputRoot | Out-Null
$log = Join-Path $OutputRoot ("run_all_subjects_{0}.log" -f (Get-Date -Format 'yyyyMMdd_HHmmss'))

function q([string]$s) { "'" + $s.Replace("'", "''") + "'" }
$Subjects = $Subjects | ForEach-Object { $_ -split ',' } | Where-Object { $_ }   # -File passes "a,b" as one string
$subj = '{' + (($Subjects | ForEach-Object { q $_ }) -join ',') + '}'
$repo = $PSScriptRoot
$cmd  = "cd($(q $repo)); run_all_subjects($(q $DataRoot), $(q $OutputRoot), $subj, $(q $CasesFile), $(q $RecordingLog), $(q $FilesCsv))"

Write-Host "MATLAB: $cmd"
Write-Host "Log:    $log"
& $Matlab -batch $cmd 2>&1 | Tee-Object -FilePath $log
exit $LASTEXITCODE
