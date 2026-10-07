# Known methodological issues

This file tracks methodological problems in the ported detection algorithms
that are **intentionally not fixed** during the EDF-pipeline restructuring.
The goal of that work is a faithful, modular port of the existing methods
(see `src/pipeline_config.m` for the exact parameter values, each traceable
to the original script it came from), not a scientific revision. Each item
below is a candidate for a future, separately-reviewed change, module by
module, now that the code is parametrized enough to allow that.

## Outlier removal (`clean_lfp.m`, ported from `complete_pipeline_seizures.m`)

- **Fixed thresholds ignore per-animal/per-session baseline.** ±2500 µV is
  applied identically regardless of electrode, impedance, or animal.
- **`k_factor` auto-switch heuristic conflates artifacts and real events.**
  The ">0.1% outliers -> switch to fixed thresholds" rule is a coarse proxy
  for "these are probably interictal spikes, not artifacts"; it has no
  physiological basis and can go either way silently.
- **MAD = 0 fallback (`noise_uV = 1`) is arbitrary.** It prevents a
  division-by-zero-like degenerate threshold but was not derived from data.

## Seizure detection (`detect_seizures.m`, ported from `complete_pipeline_seizures.m`)

- **Min-max normalization before band-passing is sensitive to outliers
  that survive cleaning** (or to the single largest artifact in the
  recording), which rescales the whole detection metric for that file.
- **No 50 Hz (or 60 Hz) notch filter.** Line noise is not rejected before
  the [5-75 Hz] bandpass, which overlaps mains frequency in many regions.
- **Threshold = median(energy) × 10 is a heuristic, not a statistically
  derived cutoff**, and is sensitive to how much of the recording is
  already seizure activity (the median shifts with seizure burden).
- **Fixed 1 s edge trim** discards a constant amount regardless of filter
  transient length at the actual sampling rate.
- **No upper bound on seizure duration**, so two adjacent events separated
  by a brief dip below threshold, or a prolonged artifact, can be recorded
  as one very long "seizure".

## Review band (`detect_seizures_robust.m`, `ll_accept` / `ll_reject`)

- **1.90 and 1.60 come from four seizures of ONE animal (005).** Three
  video-confirmed seizures exceed 3.4 in both hemispheres; the one of
  31/03 22:29:16 sits at 1.87 and 1.96, against false positives at 1.71 and
  1.60 -- the single 1.75 cut decided that seizure on a ~2 % margin in one
  channel. The band edges are placed around those few values, not derived
  from a distribution; recalibrate on more animals before trusting them.
- **The band adds review material, it does not re-sort it.**
  `Candidates_in_band` events were discarded outright by the binary cut;
  counting them as seizures would inflate every rate. Summaries keep them
  apart and out of every default rate.
- **event_id depends on the run.** It is re-assigned chronologically on
  every run; if detections change between runs, the same number can point
  to a different event. With `overwrite = false` an existing figure file of
  the same event_id in the right folder is kept as is.

## IID / spike detection (`detect_iid.m`, ported from `IID_detection_FINAL.m`)

- **Baseline search is quantized** (100:5:130 µV steps) and clamped at
  130 µV; both the step size and the ceiling are arbitrary round numbers,
  not derived from the noise distribution.
- **`MinPeakProminence` (0.2 mV) is an absolute value**, not scaled to the
  file's own baseline/noise level, so it behaves inconsistently across
  recordings with different impedance or amplification.
- **Complexes have no maximum duration cap.** A run of spikes with each
  gap <= 150 ms can chain indefinitely into one "complex" no matter how
  long it lasts.
- **No 50 Hz notch** before the [15-70 Hz] bandpass, same caveat as
  seizure detection.
- **Exclusion-zone buffer (5 s) and burst grouping window (5 s) are fixed**
  and not related to any measured property of the recording (e.g. seizure
  post-ictal suppression length, which varies by animal/seizure severity).

## Bilateral reconciliation (`bilateral_reconcile.m`)

- **Systematic rescue can mask genuinely unilateral seizures.** With
  `rescue_mode = 'rescue_and_impute'` every event accepted in one
  hemisphere is *reported* in all of them, so row counts, per-channel
  `n_seizures_reported` and `*_reported` time-in-seizure are identical
  across hemispheres by construction -- including for a seizure that was
  truly focal. Laterality must be read from `accepted_in_n_channels` /
  `is_bilateral_accepted` / `accepted_in_regions` (and each row's own
  `ll_ratio`), and per-channel rates from the `*_accepted` columns, never
  from the `*_reported` ones.
- **An imputed row is not a measurement of a seizure in that channel.** It
  is the other channel's window, measured here; its metrics say what that
  window looks like in this hemisphere, not that a seizure happened there.
- **Only line-length rejections can be rescued.** Crossings the robust
  branch dropped as too short after merging (< `min_duration_s`) are not
  exposed by `detect_seizures_robust.m`, so such an event is imputed, not
  rescued, even if the channel had a near-miss there.
- **A failed channel is reported, not measured.** When a robust channel's
  detection fails (e.g. 005 `20260326 mTor` HPCr: flat signal, 99.98 %
  outliers), each event of that recording gets an `imputed` row there with
  `rejected_by='detection_failed'` and NaN metrics. Such rows must be
  excluded from any per-channel metric analysis; they only say "this
  hemisphere could not be assessed for this event".
- **Robust branch only.** Legacy channels are not reconciled; in a
  recording where one channel is robust and another legacy (different
  cases per channel), the legacy one neither receives nor contributes
  events, and its ids stay per-channel.
- **`match_tol_s = 5` s is a round number**, not derived from measured
  inter-hemispheric propagation delays; a larger value merges more nearby
  but distinct events, a smaller one splits one event in two.

## Campaign settings (`campaign_config.m`, `compute_run_reference.m`)

- **The robust branch now runs on every channel**, so its parameters
  (calibrated on 4 seizures of ONE animal, 005 -- see above) apply to the
  whole cohort, not only to a few problem channels.
- **The gain reference depends on the run.** It is the median amplitude of
  the non-attenuated channels of the EDFs processed together, so a
  preliminary run (a few EDFs) and the full campaign give attenuated
  channels different gains, hence different cleaning and IID results for
  them. Seizure detection is barely affected (relative thresholds).
- **The reference falls back silently in quality.** A region with no
  non-attenuated channel in the run takes the all-regions median
  (`reference_source = 'fallback'`); a run with none at all applies no
  gain (`reference_source = 'none'`). Both are flagged in
  `qc_report.csv`'s warnings, never an error.
- **Unknown Attenuation values count as attenuated.** Anything in the
  Excel's Attenuation column other than no / none / empty / severe / mild
  / "unsure, perhaps mild" (e.g. a typo) gets gain, with a warning.
- **Attenuation is per (animal, EDF), not per channel.** The Excel column
  does not say which channel is attenuated, so both channels get the gain.
- **The notch reaches IID too.** It is applied when cleaning, before both
  detectors; removing 49-51 Hz slightly lowers the amplitude of spikes and
  of the IID baseline wherever mains noise was present -- intended, but it
  makes IID numbers not directly comparable with un-notched runs.
- **`evaluate_detections.m` is not reconciliation-aware** (see README
  "Not implemented").

## Cross-cutting

- **No cross-validation against Natus's own visual/automatic scoring** is
  built into the detectors themselves; `natus_review_sheet.csv` (planned,
  not yet built) is meant to make that comparison possible downstream, but
  the detectors do not use Natus output as ground truth or feedback.
