# EDF pipeline (Natus)

MATLAB pipeline that takes Natus EDF(+C/+D) exports and produces seizure
and interictal-spike (IID) detections, ready to compare against a manual
Natus review. It is a restructured, methodology-preserving port of the
original Intan-era scripts at the repo root (`complete_pipeline_seizures.m`,
`IID_detection_FINAL.m`, `edf_txt_conversion_Samara.m`, ...), which are
kept untouched as reference. See `KNOWN_ISSUES.md` for the methodological
caveats that were deliberately carried over unchanged.

## Flow

```
EDF (Natus)
  |
  |  src/edf_import.m            one txt per selected channel (continuous,
  v                               NaN-padded at real gaps) + one gaps CSV
01_txt/{subject}_{yyyyMMdd}_{HHMMSS}_{region}.txt
01_txt/{subject}_{yyyyMMdd}_{HHMMSS}_gaps.csv
  |
  |  src/clean_lfp.m             amplitude-threshold outlier removal,
  v                               per gap-free block
02_clean/{subject}_{yyyyMMdd}_{HHMMSS}_{region}_clean.txt

  |  src/detect_seizures.m       Hilbert-envelope detector, per gap-free
  v                               block, global threshold
03_seizures/..._seizures.mat / .fig / .png

  |  src/detect_iid.m            spike/polyspike/burst detector, exclusion
  v                               zones = seizures + gaps + manual
04_iid/..._IID_results.mat / .fig / .png

  |  src/run_pipeline_edf.m      consolidates every file/channel above
  v
05_summaries/ (9 CSVs + pipeline_summary.xlsx, see below)
logs/run_<timestamp>.log
config_used.mat / config_used.json
```

Every stage is a plain function that takes a config struct and the
previous stage's output and returns a struct/table -- none of them depend
on global state, and each can be called on a single file to debug it in
isolation (see **Invocation** below). `run_pipeline_edf.m` is the only
piece that walks a folder, writes the final CSVs, and is resilient to a
single file or channel failing (logged, not fatal -- see `qc_report.csv`).

### Design decisions worth knowing before you read the code

- **No splitting at discontinuities.** Unlike `edf_txt_conversion_Samara.m`,
  a gap (real acquisition discontinuity) does **not** start a new file.
  Each channel's txt is one continuous, NaN-padded timeline from the first
  to the last sample of the EDF; `t_rel = (i-1)/fs` for every sample index
  `i`, gaps included. The companion `*_gaps.csv` records where those NaN
  spans are (see schema below).
- **A gap is different from jitter.** `cfg.edf.gap_tol_s` (default 1e-3 s)
  is the threshold used to even flag a record-time deviation as a
  candidate (same rule as `edf_txt_conversion_Samara.m`); `cfg.edf.gap_min_s`
  (default 0.5 s) is the threshold that separates a real gap from
  sub-sample timestamp jitter. Jitter is logged, never turned into NaN.
- **Sample and time conservation are exact, not approximate.** `edf_import.m`
  places each record's samples contiguously and only advances the write
  cursor by `round(gap_duration_s * fs)` NaN samples at a *real* gap
  (never at jitter). This guarantees `n_valid_samples == n_records *
  samples_per_record` and `valid_duration_s + sum(gaps) == total_duration_s`
  exactly, by construction, not within a tolerance.
- **Filtering, peak detection, and event grouping never cross a gap.**
  `detect_seizures.m` and `detect_iid.m` run their bandpass/Hilbert/
  findpeaks per contiguous valid block (`mask_to_segments` +
  `filter_by_blocks`), and additionally force a new complex/burst boundary
  whenever consecutive spikes/complexes fall in different blocks -- so no
  seizure, complex, or burst can span a real gap even for parameter values
  where the ISI/gap-time threshold alone would have merged them. This was
  verified with a synthetic test (the sample EDF has no real gaps to
  exercise the path with real data -- see **Testing status** below).
- **Every threshold that should be "global" stays global.** The seizure
  energy threshold and the IID baseline are single statistics computed
  over the concatenation of all valid blocks (seizures) or all included
  samples (IID), never recomputed per block.
- **IID exclusion zones are automatic**, built from detected seizures
  (+`cfg.iid.exclusion_buffer_s`), the recording's own gaps, and optional
  manual zones (`cfg.iid.exclusion_zones_manual`) -- see `iid_summary.n_exclusion_zones`
  and the `reason` column inside each channel's own exclusion table
  (`iid_results.exclusion_table`, not written to a CSV by itself; the
  final per-event tables already reflect it).
- **Channel-to-region mapping is a placeholder.** `pipeline_config.m`
  ships with `cfg.edf.channels.mode = 'list'` and the raw EDF labels found
  in the shipped test file (`A7C1`, `A7C3`) used as-is for `region` (no
  invented anatomy). Switch to `cfg.edf.channels.mode = 'map'` and edit
  `cfg.edf.channels.map` with the real label -> region names once you've
  confirmed the montage for a given study.
- **Per-file mapping from a recording log.** `cfg.edf.channels.mode = 'log'`
  reads `cfg.edf.channels.log_file` (`EEG_recording_log.xlsx`: `Filename`,
  `Animal ID`, `Port`, `HPCr_channel`, `HPCl_channel`) and, for each EDF,
  picks the row matching (`cfg.edf.subject_id`, EDF filename) -- never the
  filename alone, since animals recorded together share one filename. Port +
  channel gives the EDF label (`A1` + `C2` -> `EEG A1C2`, or `A1C2` in older
  exports); regions are named `cfg.edf.channels.log_regions` (`HPCr`, `HPCl`)
  regardless of port, so runs stay comparable when an animal changes port.
  An EDF with no log row, or whose logged channels aren't in the file, fails
  `edf_import` and is reported in `qc_report.csv`. `run_all_subjects.m` /
  `run_all_subjects.ps1` run every subject folder this way.

### Testing status

Validated end to end against `097-s/*.EDF` (continuous, no real gaps):
sample/time conservation, block-wise vs. monolithic numerical equivalence
(acceptance criterion 5), no detected event overlapping an exclusion zone,
no empty identity columns, and absolute-time agreement across CSVs
(criteria 4 and 7) all pass. The gap-handling code path itself (cursor
placement, jitter/anomaly classification, complex/burst block-boundary
guards) was verified with synthetic unit tests, since the only EDF
available has no real discontinuities. **Acceptance criterion 6** (a real
EDF+D file processing cleanly, with its gaps CSV matching
`edf_txt_conversion_Samara.m`'s cut points) is still unverified against
real data -- run it against a genuinely gapped EDF before trusting that
path in production.

### Resumability (scoped)

`cfg.general.overwrite = false` makes every stage refuse to clobber an
existing output file. `run_pipeline_edf.m` additionally *skips recomputing*
stage 2 (`clean_lfp`) when its output already exists, because that's the
stage where avoiding recomputation is worth it and its output filename is
a trivial, low-risk one-line rule (`<raw_basename>_clean.txt`). Stages 1
(`edf_import`, the expensive `edfread` call), 3, and 4 still run every
time; they just won't overwrite what's already on disk. This is a
deliberate scope limitation, not a bug.

## Header field dictionary

Every txt file (`01_txt/*.txt` and `02_clean/*_clean.txt`) starts with
`# key = value` lines, parsed generically by `load_lfp_txt.m` (any field,
not just a hardcoded few) via `parse_header.m`. All values are read back
as text; callers convert what they need (see `load_lfp_txt.m`'s handling
of `fs`).

| Field | Written by | Meaning |
|---|---|---|
| `subject_id` | edf_import | From `cfg.edf.subject_id`, or derived from the EDF filename (logged as a warning) if left empty |
| `fs` | edf_import | Sampling frequency (Hz) for *this channel*: `NumSamples(channel) / DataRecordDuration` |
| `time_unit` | edf_import | Always `seconds` |
| `region` | edf_import | Resolved region name (see channel mapping above) |
| `channel_label` | edf_import | The original EDF signal label (e.g. `A7C1`) |
| `source_file` | edf_import | EDF filename (not full path) |
| `source_format` | edf_import | `EDF+C` or `EDF+D`, from the EDF header's reserved field |
| `session_start_datetime` | edf_import | ISO 8601 with UTC offset, e.g. `2026-08-19T12:20:22.734+02:00` -- the instant of *sample 1*, which can differ from the EDF header's nominal StartTime by the first record's own onset |
| `session_start_unix` | edf_import | POSIX timestamp of the same instant, full `double` precision (`%.17g`) -- the authoritative source `load_lfp_txt.m` reconstructs `session_start` from; the ISO string is for humans |
| `timezone` | edf_import | From `cfg.general.timezone` (default `Europe/Paris`) |
| `n_samples` | edf_import | Length of the continuous (NaN-padded) vector written to this file |
| `n_valid_samples` | edf_import | Samples that actually came from an EDF record (`n_samples - sum(gap samples)`) |
| `n_gaps` | edf_import | Number of real gaps in this file (see the sibling `*_gaps.csv`) |
| `total_duration_s` / `valid_duration_s` | edf_import | `n_samples/fs` and `n_valid_samples/fs` |
| `units` | edf_import | `microvolts`, unless the channel's EDF PhysicalDimensions wasn't a recognized voltage unit (only reachable via `cfg.edf.channels.mode='all'`) |
| `columns` | edf_import | `amplitude_microvolts` (data section is one amplitude value per line, no time column) |
| `clean_source_file` | clean_lfp | The raw txt this clean file was built from |
| `threshold_type` | clean_lfp | `FIXED (user config)`, `FIXED` (auto-switched), or `AUTOMATIC (k=...)` |
| `k_factor`, `lower_threshold_uV`, `upper_threshold_uV` | clean_lfp | Outlier thresholds actually used |
| `auto_switch`, `auto_switch_reason` | clean_lfp | Whether the >0.1% preliminary-outlier rule forced fixed thresholds, and why |
| `n_outliers`, `outlier_pct` | clean_lfp | Outliers found (relative to valid, i.e. non-gap, samples) |
| `signal_range_original_uV`, `signal_range_clean_uV` | clean_lfp | `[min, max]` before/after cleaning (valid samples only) |

All other fields present in a raw txt (everything above `clean_source_file`)
are copied through into the clean txt unchanged -- subject/region/session
identity never gets lost at this step.

## CSV schemas (`05_summaries/`)

All nine outputs carry full identification (`subject_id`, `region` where
applicable, `source_file`) and both relative (`_s`) and absolute (`_abs`)
times where a time is meaningful.

**seizures_events.csv** -- one row per seizure
`subject_id, region, session_start, source_file, seizure_id, start_s, end_s, duration_s, start_abs, end_abs, block_id, adjacent_to_gap`

**seizures_summary.csv** -- one row per file/channel
`subject_id, region, session_start, source_file, total_duration_min, valid_duration_min, n_gaps, gap_duration_min, n_seizures, total_seizure_time_s, pct_time_in_seizure, mean_duration_s, min_duration_s, max_duration_s, median_energy, threshold_value, pct_above_thr, n_segments, n_rejected, bandpass_low, bandpass_high, power_exponent, window_s, median_factor, min_seizure_duration`

**iid_events.csv** -- one row per spike complex, classified
`subject_id, region, session_start, source_file, complex_id, start_s, end_s, start_abs, end_abs, duration_ms, n_spikes, classification ('single'|'polyspike'), max_amplitude_mV, mean_amplitude_mV, in_burst, burst_id`

**iid_summary.csv** -- one row per file/channel
`subject_id, region, session_start, source_file, total_duration_min, analyzed_duration_min, excluded_duration_min, n_exclusion_zones, baseline_uV, lower_threshold_uV, upper_threshold_uV, total_peaks, total_complexes, n_single, n_polyspike, pct_polyspike, complexes_per_min, single_per_min, polyspikes_per_min, n_bursts, bursts_per_hour, mean_spikes_per_polyspike`

**iid_bursts.csv** -- one row per burst
`subject_id, region, start_s, end_s, start_abs, end_abs, duration_s, n_complexes`

**gaps_summary.csv** -- every gap, every file, consolidated
`gap_id, start_s, end_s, duration_s, start_abs, end_abs, prev_record_idx, next_record_idx, source_file`

**qc_report.csv** -- one row per file/channel
`subject_id, region, source_file, session_start, stages_completed, n_errors, error_messages, n_blocks_rejected_short, outlier_pct, nan_pct, n_warnings, warning_messages`

**pipeline_summary.xlsx** -- the same seven tables above (everything
except `natus_review_sheet`) as sheets named `seizures_events`,
`seizures_summary`, `iid_events`, `iid_summary`, `iid_bursts`, `gaps`, `qc`.

**natus_review_sheet.csv** -- built for the manual Natus comparison:
every seizure, IID burst, and gap merged and sorted by absolute time
`abs_time, clock_time, event_type ('seizure'|'iid_burst'|'gap'), duration_s, region, subject_id, natus_confirmed`
(`natus_confirmed` is left blank for the reviewer to fill in by hand).

### Re-reading these CSVs safely

Plain `readtable()` is **not** safe on any file above: MATLAB's automatic
column-type detection turns a `subject_id` like `"097"` into the number
`97`, and a file with zero rows (e.g. zero gaps) comes back with its
datetime columns typed as a plain `struct` instead of `datetime`, because
there is nothing in the file for `readtable` to infer a format from. Use
`read_pipeline_csv.m` instead, which forces every column to the type the
pipeline actually writes:

```matlab
se = read_pipeline_csv('pipeline_output/05_summaries/seizures_events.csv', 'seizures_events');
```

`kind` (the 2nd argument) is one of `seizures_events`, `seizures_summary`,
`iid_events`, `iid_summary`, `iid_bursts`, `gaps`, `qc` -- matching the
CSV you're reading. A 3rd argument sets the TimeZone to reconstruct on
`start_abs`/`end_abs`/`session_start` (default `Europe/Paris`, must match
whatever `cfg.general.timezone` the run actually used -- plain CSV text
carries no timezone of its own).

## Multiple subjects / multiple recording days

`run_pipeline_edf.m` uses a single `cfg.edf.subject_id` for every EDF it
processes in one call. If a folder mixes EDFs from different mice, either
name the EDF files so `cfg.edf.subject_id = ''` can derive a correct,
distinct ID per file from each filename, or -- the more reliable option --
run `run_pipeline_edf` once per subject (its own input folder, its own
`cfg.edf.subject_id`, its own `cfg.paths.output_root`) and combine the
results afterward with `merge_pipeline_runs.m`:

```matlab
r1 = run_pipeline_edf('data/mouse_097', set_subject(pipeline_config(), '097', 'out/097'));
r2 = run_pipeline_edf('data/mouse_098', set_subject(pipeline_config(), '098', 'out/098'));

merged = merge_pipeline_runs({'out/097', 'out/098'}, 'out/merged');
merged.seizures_events           % both subjects, correctly identified
merged.paths.natus_review_sheet  % one combined review sheet, chronological

function cfg = set_subject(cfg, id, out_dir)
    cfg.edf.subject_id = id;
    cfg.paths.output_root = out_dir;
end
```

Recordings of the *same* mouse on different days need no special
handling: `session_start` (from the EDF's own header) and the output
filenames already differentiate sessions, so they can all go through one
`run_pipeline_edf` call with one fixed `subject_id`.

## Invocation

```matlab
addpath('src');
addpath('src/utils');

cfg = pipeline_config();
cfg.edf.subject_id = '097';                       % '' would derive it from the EDF filename
cfg.paths.output_root = fullfile(pwd, 'pipeline_output');
% cfg.edf.channels.mode = 'map';                  % once you know the real region names:
% cfg.edf.channels.map  = containers.Map({'A7C1','A7C3'}, {'HPCleft','HPCright'});

result = run_pipeline_edf('097-s', cfg);           % folder of .edf files

% Everything above is also in memory, not just on disk:
result.seizures_events
result.iid_summary
result.paths.natus_review_sheet                    % path to the CSV
```

To debug a single file/channel without running the whole batch:

```matlab
cfg = pipeline_config();
cfg.edf.subject_id = '097';
cfg.edf.output_dir = 'tmp/01_txt';
[manifest, gaps, meta] = edf_import('097-s/some_file.EDF', cfg);

raw = load_lfp_txt(manifest.txt_file{1});
cfg.clean.output_dir = 'tmp/02_clean';
clean_result = clean_lfp(raw, cfg);
clean = load_lfp_txt(clean_result.txt_file);

cfg.seizure.output_dir = 'tmp/03_seizures';
seizures = detect_seizures(clean, cfg);

cfg.iid.output_dir = 'tmp/04_iid';
iid = detect_iid(clean, seizures.seizures, cfg);    % omit/[] the 2nd arg to skip seizure-based exclusion
```

## Repository layout

```
src/
  pipeline_config.m          all parameters, one struct, traceable to source scripts
  edf_import.m                EDF -> continuous txt per channel + gaps CSV
  load_lfp_txt.m               generic txt reader (replaces load_LFP_intan_txt.m)
  clean_lfp.m                  outlier removal -> *_clean.txt
  detect_seizures.m            seizure detection
  detect_iid.m                  interictal spike/polyspike/burst detection
  run_pipeline_edf.m           batch orchestrator over a folder of EDFs
  merge_pipeline_runs.m         combine several run_pipeline_edf.m runs (e.g. multiple subjects)
  utils/
    write_lfp_txt.m, parse_header.m      generic txt header read/write
    mask_to_segments.m, filter_by_blocks.m   gap-respecting block utilities
    rel_to_abs_time.m                     relative seconds -> absolute datetime
    ensure_findpeaks_signal_toolbox.m     forces Signal Toolbox findpeaks over Chronux
    read_pipeline_csv.m                   safe reader for the 05_summaries/ CSVs
    write_all_summaries.m, build_natus_review_sheet.m   shared writer, used by both
      run_pipeline_edf.m and merge_pipeline_runs.m
    apply_tz.m, vertcat_or_empty.m, empty_*_table.m     table-schema plumbing shared by
      the writer and read_pipeline_csv.m (single source of truth, see their headers)
KNOWN_ISSUES.md               methodological caveats carried over unchanged, on purpose
```

The Intan-era scripts at the repo root and `functions/dual_convert_txt_dat/`
are unmodified reference material (format/logic), not part of this
pipeline's call graph.
