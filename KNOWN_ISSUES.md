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

## Cross-cutting

- **No cross-validation against Natus's own visual/automatic scoring** is
  built into the detectors themselves; `natus_review_sheet.csv` (planned,
  not yet built) is meant to make that comparison possible downstream, but
  the detectors do not use Natus output as ground truth or feedback.
