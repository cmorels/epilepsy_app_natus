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

## Review band and event categories (robust branch)

### Review band

`detect_seizures_robust.m`'s final filter is no longer a single cut at
`ll_threshold`: each candidate gets an `ll_status`

| `ll_status` | rule |
|---|---|
| `accepted` | `ll_ratio >= cfg.seizure_robust.ll_accept` (default 1.90) |
| `in_band` | `ll_reject <= ll_ratio < ll_accept` -- kept, needs review |
| `rejected` | `ll_ratio < cfg.seizure_robust.ll_reject` (default 1.60) |

`.seizures` holds accepted + in_band rows, `.rejected_events` the rejected
ones, `.candidates` all of them (`kept = ll_status ~= 'rejected'`). The
line-length panel shows both lines with the band shaded. Legacy rows
(`detect_seizures.m`, untouched) are always `ll_status = 'accepted'`: the
band is a robust-branch concept only. **Backward compatibility:** set
`ll_accept` or `ll_reject` to `[]`/NaN and the band is OFF -- every output
is byte-identical to the binary pipeline; `ll_accept = ll_reject = 1.75`
keeps the machinery on with a zero-width band and reproduces the binary
classification.

### Events and categories (`src/bilateral_events.m`)

Robust detections of the same **(subject_id, session_start)** -- never the
EDF file, which in log mode holds several animals -- are grouped across
channels (overlap or gap < `cfg.seizure_robust.bilateral_tol_s`,
transitive; two rows of the same channel are never merged, they share the
`event_id` with `fragmented = true`). Each event gets a new **`event_id`**
(1..K per animal + session, chronological, across all categories;
`seizure_id` keeps its meaning and values as the channel-local index) and a
category, counting channels:

| category | rule |
|---|---|
| `Crisis` | accepted in >= 2 channels |
| `Candidates` | accepted in exactly 1 channel |
| `Candidates_in_band` | accepted in none, in_band in >= 1 |

**Warning:** `Candidates_in_band` is material the binary pipeline
discarded outright -- the band does not reorder what existed, it adds new
material to review. It is never mixed into Crisis counts or rates.

### Hemisphere mapping

The hemisphere comes ONLY from the **region** name (`cfg.edf.channels.mode =
'log'` names them `HPCr` = right, `HPCl` = left), through
`cfg.output.hemisphere.left_regions` / `right_regions` (case-insensitive,
trimmed; default `{'HPCl','HPCleft','HPC_L'}` / `{'HPCr','HPCright','HPC_R'}`).
Never from the channel label (A5C2 / A7C1 ... change between recordings and
animals) nor from processing order. A region in neither list is processed
anyway (`hemisphere = 'unknown'`) and gets an extra joint-figure column on
the right, with a log warning.

### Figures

```
03_seizures/                       (per-channel .mat and panoramas stay here)
  Crisis/              individual/   joint/
  Candidates/          individual/   joint/
  Candidates_in_band/  individual/   joint/
individual: {subject}_{session}_event{ID:02d}_{region}.png/.fig   (every robust channel, detected or not)
joint     : {subject}_{session}_event{ID:02d}_joint.png/.fig      (cfg.output.joint_figures)
```

Joint figure = 4 rows (voltage, band-passed, energy + threshold, line
length + band) x N channels (left hemisphere, right, then unknown), same
time window and linked x on every panel, **same y-limits across columns
within each row** (energy row in log scale on every column if the maxima
differ by > 100x), the event's reference window shaded per column by that
channel's `ll_status` (or "sin detección"), column headers with region,
hemisphere, channel, status, `ll_ratio` and distance to the nearest
threshold, and the absolute start time in the title. Individual and joint
figures are drawn by the same `utils/plot_trace_*.m` functions. Stale files
of the same subject + session (event gone or re-categorized) are deleted;
with `overwrite = false` existing correct files are not redrawn.

### CSVs

- `seizures_events.csv`: `ll_status` + `event_id, category, n_accepted,
  n_in_band, n_channels_in_group, fragmented, hemisphere,
  accepted_in_regions, ref_start_s, ref_end_s, figure_individual_path,
  figure_joint_path` (a legacy-only run adds `ll_status` only).
- `seizures_summary.csv`: pre-existing columns count **accepted rows only**;
  appended `n_events_crisis, n_events_candidates, n_events_candidates_in_band,
  n_accepted_this_channel, n_in_band_this_channel`, and time in seizure
  twice: `*_crisis` (Crisis only) and `*_crisis_candidates` (Crisis +
  Candidates). `Candidates_in_band` never enters a rate.
- `natus_review_sheet.csv`: one row per (event_id, channel), grouped by
  event, with `category, ll_status, ll_ratio, review_priority` (1 =
  Candidates, 2 = Candidates_in_band, 3 = Crisis), `figure_joint_path`;
  `natus_confirmed` is the reviewer's verdict column.
- IID: in_band rows do not exclude time unless
  `cfg.bilateral.exclude_in_band_from_iid = true` (default false).
- `merge_pipeline_runs.m` propagates everything; the event key across runs
  is (`subject_id`, `session_start`, `event_id`), never renumbered.

The review band and `cfg.bilateral.rescue_mode ~= 'off'` cannot be combined
yet (the run stops with a clear message); `bilateral_events.m` marks the
extension point where rescued/imputed rows would plug in.

## Bilateral reconciliation (reconciliación bilateral)

Without it, each channel (hemisphere) is detected and numbered on its
own: the same seizure can pass the filters in one hemisphere and fall just
below them in the other (implant quality, attenuation, sitting near a
threshold), so the two channels of one animal report different counts and
"seizure 3" of one has nothing to do with "seizure 3" of the other.

`src/bilateral_reconcile.m` runs after detection on every channel of one
recording and before any seizure CSV, figure or IID exclusion zone:
**a seizure is dropped only if every channel rejected it; if at least one
channel accepted it, it is reported in all of them**, under one shared
`seizure_id` (chronological, 1..K per animal + session).

**It applies to the robust branch only** (channels whose case resolved to
`seizure_mode = 'robust'`, i.e. `attenuated` / `line` / `both`). Legacy
channels (`detect_seizures.m`, case `normal`) are never reconciled nor
used as a reference: they keep their own figures and per-channel ids, and
their rows are tagged `detection_status = 'not_reconciled'`.

```matlab
cfg.bilateral.rescue_mode = 'rescue_and_impute'; % RECOMMENDED. 'off' (default) | 'rescue' (diagnosis only)
cfg.bilateral.match_tol_s = 5;                   % detections of different channels closer than this = one event
cfg.bilateral.exclude_rescued_from_iid = true;   % IID exclusion zones also cover rescued/imputed rows
cfg.bilateral.require_same_subject = true;       % never group different animals
```

Per event and channel, `detection_status` says how that row got there:

| status | meaning | window |
|---|---|---|
| `accepted` | this channel's detector accepted it on its own merit | the detector's own |
| `rescued` | this channel had it as a candidate rejected by the line-length filter (`kept=false` in `detect_seizures_robust.m`'s `.candidates`, `rejected_by='ll_ratio'`) -- the one with the highest `ll_ratio` if several; a candidate is rescued into one event at most | the candidate's OWN limits in this channel |
| `imputed` | no rejected candidate of this channel overlapped the event | the event's reference window (union of the accepted rows), clipped to one valid block |
| `imputed`, `rejected_by='detection_failed'` | this (robust) channel's clean or detection stage FAILED -- the failure is still reported, once per event; metrics are NaN (nothing to measure) and no figure is drawn; the channel counts in `n_channels_total` | the event's reference window |
| `not_reconciled` | legacy channel, or reconciliation failed for that recording | the detector's own |

`ll_ratio`, `peak_energy_ratio`, `hf_ratio_db`, `envelope_cv` and the
duration are always computed on the row's own window in its own channel,
never copied (a rescued `ll_ratio=1.62` against the 1.75 cut says the event
was *close* in that hemisphere, not absent). In `'rescue'` mode, events with
no candidate get no row, so channels can end up with different counts --
that mode is for diagnosis, not normal use.

**Unit of analysis changes.** After reconciliation the unit is the *event
per animal*, not the *detection per channel*. That is why "accepted" and
"reported" columns coexist and must not be mixed:

- `seizures_events.csv`: `seizure_id` becomes the shared id;
  `channel_seizure_id` keeps the detector's own per-channel index (NaN on
  rescued/imputed rows); appended columns `detection_status`,
  `accepted_in_n_channels`, `accepted_in_regions`, `n_channels_total`,
  `is_bilateral_accepted` (accepted in >= 2 channels), `fragmented` (this
  channel contributes more than one row to the same event -- kept as
  separate rows, never merged), `ref_start_s`/`ref_end_s` (reference
  window, same on every channel), `rejected_by`. Whether an event was
  bilateral or unilateral is answered by `accepted_in_n_channels`, never
  by the row count.
- `seizures_summary.csv`: every pre-existing column still describes what
  THIS channel accepted (so `n_seizures` = `n_seizures_accepted`); appended
  `n_seizures_accepted`, `n_seizures_reported`, `n_events_reported`,
  `n_rescued`, `n_imputed`, and time-in-seizure twice:
  `total_seizure_time_s_accepted` / `_reported`,
  `pct_time_in_seizure_accepted` / `_reported`.
- `natus_review_sheet.csv`: one row per (seizure_id, channel), grouped by
  event, with `seizure_id`, `detection_status`, `accepted_in_regions`,
  `source_file`.
- `03_seizures/`: figures are drawn by `utils/save_bilateral_seizure_figures.m`
  (the detectors' own per-channel figures are discarded):
  `{base}_seizure{seizure_id:02d}.png/.fig`, titled e.g.
  `Crisis 3 (HPCr, accepted)` / `Crisis 3 (HPCl, rescued, ll_ratio=1.62 < 1.75)` /
  `Crisis 3 (HPCl, imputed: no candidate in this channel)`; red = accepted,
  orange = rescued, blue = imputed (legend on every figure), reference
  window dotted; rescued/imputed outlined dashed on the panorama.
- IID: with `exclude_rescued_from_iid = true` the exclusion zones include
  rescued/imputed rows (their `source_id` is the shared `seizure_id`).
- `merge_pipeline_runs.m` propagates every column and never renumbers:
  across runs the key is (`subject_id`, `source_file`/`session_start`,
  `seizure_id`).

With `rescue_mode = 'off'` none of this runs and every output is
byte-identical to the pipeline before this stage existed. See
KNOWN_ISSUES.md for why systematic rescue can mask genuinely unilateral
seizures.

## Invocation

```matlab
addpath('src');
addpath('src/utils');

cfg = pipeline_config();
cfg.edf.subject_id = '097';                       % '' would derive it from the EDF filename
cfg.paths.output_root = fullfile(pwd, 'pipeline_output');
% cfg.edf.channels.mode = 'map';                  % once you know the real region names:
% cfg.edf.channels.map  = containers.Map({'A7C1','A7C3'}, {'HPCleft','HPCright'});
cfg.bilateral.rescue_mode = 'rescue_and_impute';  % RECOMMENDED: one shared seizure_id per event across hemispheres
                                                  % (default 'off' = per-channel, unchanged; see "Bilateral reconciliation")

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
