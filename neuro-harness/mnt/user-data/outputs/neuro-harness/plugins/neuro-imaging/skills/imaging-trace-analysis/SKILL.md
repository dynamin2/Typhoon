---
name: imaging-trace-analysis
description: Use when analyzing presynaptic optical recordings — iGluSnFR glutamate traces, pHluorin (SypH, vGlut-pHluorin) exocytosis/endocytosis traces, GCaMP or FRCaMPi calcium traces, dual Ca2+/glutamate recordings, stimulus-train responses, event detection, or latency measurement. Covers baseline/ΔF-F0 handling, bleach correction, train depression fitting, latency quantization limits, and the checks to run before interpreting a kinetic parameter.
---

# Presynaptic imaging trace analysis

## Before any fitting, report the acquisition constraints

State the frame interval and the number of stimuli/frames. Every kinetic
parameter is bounded by them:

- A time constant shorter than ~2 frame intervals is not measured, it is
  extrapolated. Say so rather than reporting it.
- Latency measured at frame interval Δt is **quantized** at Δt. A latency
  distribution on a 20 ms grid has at most a handful of distinct values in the
  fast range — density estimates and mixture fits on such data must account for
  the discreteness, and any claimed component separation smaller than Δt is
  not resolvable.

Dual-channel recordings trade temporal resolution for simultaneity. When the
frame interval is the limiting factor because two indicators are recorded
together, note this in the results — it is a design constraint, not a
shortcoming to be hidden.

## ΔF/F₀ and bleaching

- Define F₀ explicitly (which frames, mean or median) and keep it identical
  across compared conditions.
- Fit and remove photobleaching from a pre-stimulus and, where available,
  post-response baseline segment. Report the bleach time constant. If bleaching
  is comparable to the response decay you are fitting, the decay estimate is
  confounded — say so.
- Never normalize to a peak that is itself the quantity being compared.

## Event detection

Report the detection threshold in units of baseline SD, the minimum event
separation, and the resulting false-positive estimate from a signal-free
segment. Apply the identical detector to every condition. If a condition
required different settings, that is a result to disclose, not a preprocessing
choice to bury.

## Stimulus-train depression

For pulse-by-pulse amplitude decay across a train:

- Fit per-pulse amplitudes, not the raw cumulative trace.
- Report whether a single exponential is adequate before introducing a second;
  apply `mixture-model-guardrails` Rule 1 logic to the model comparison.
- Give the decay constant in **pulse number and in seconds**, and state the
  stimulation frequency. τ in pulses and τ in time are different claims about
  mechanism.
- Distinguish depression of release from indicator saturation or
  reuptake/clearance kinetics. State which the data can and cannot separate.

## Latency between two indicators

- Define the timing landmark on each channel identically (e.g. both at half-rise
  or both at onset), and state which.
- Physically impossible latencies (downstream signal preceding upstream) mark
  detection failures, not fast events — exclude them and report how many.
- Trials with no detectable signal in the reference channel cannot yield a
  latency. Handle them as a separate declared category with its own n; do not
  merge them into the latency distribution and do not silently drop them, since
  their frequency is itself informative.

## Always state what the measurement cannot distinguish

Close every kinetic analysis with one sentence naming the mechanisms the data
do not separate. This is a required part of the output, not a caveat to omit.
