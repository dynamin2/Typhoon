---
name: mixture-model-guardrails
description: Use whenever fitting Gaussian mixture models, selecting the number of components, or making a claim that a distribution has multiple populations — amplitude histograms, latency distributions, quantal analysis, event-size distributions. Also use when interpreting AIC/BIC output, comparing conditions on component counts, or when asked whether a treatment "induces a new population". Enforces null calibration before any multi-component claim.
---

# Mixture-model guardrails

The failure mode this exists to prevent: reporting that a condition has *k*
components because AIC or BIC picked *k*, when the information criterion had no
calibrated null to be compared against and the sample size could not have
resolved *k* in the first place.

## Rule 1 — an information criterion alone never establishes a component count

AIC/BIC differences are not a test. Before reporting that a distribution has
more than one component, run a **parametric bootstrap null calibration**:

1. Fit the *k*-component and (*k*+1)-component models to the real data; record
   the observed ΔBIC (or ΔAIC, or the log-likelihood ratio).
2. Simulate ≥1000 datasets **of the same n** from the fitted *k*-component
   model (the null).
3. Refit both models to each simulated dataset and build the null distribution
   of Δ.
4. Report the empirical p-value: the fraction of null Δ at least as extreme as
   observed.

Report the calibrated p-value, not the raw ΔBIC, as the basis for the claim.
Use the same calibration procedure across every dataset in a manuscript — an
amplitude histogram and a latency distribution in the same paper must not be
held to different standards.

## Rule 2 — check that n could resolve the claim

Two-component resolution in these distributions typically needs several hundred
events; a few dozen cannot support it regardless of what BIC says. If n is
below what the calibration in Rule 1 shows is resolvable, state that the data
are **underpowered to distinguish 1 from 2 components** and stop. Do not report
the fitted two-component parameters as if they were estimates.

## Rule 3 — figure and model selection must agree, and disagreement is a finding

If the selected component count contradicts what the histogram visibly shows —
or if a condition that looks more heterogeneous is assigned *fewer* components
than the control — do not narrate around it. Halt and report the contradiction
explicitly, then check, in this order:

- Were the two conditions filtered by the same criteria? Asymmetric inclusion
  (e.g. dropping zero or negative values in one group but not the other) will
  produce exactly this artifact.
- Is the fit converging to a degenerate solution (a component collapsing onto
  a few points, near-zero variance)? Report convergence status.
- Is the binning of the displayed histogram doing the visual work rather than
  the data?

## Rule 4 — asymmetric filtering must be declared

For every group, report: total events, how many were removed, and by which
criterion. Groups compared against each other must pass through the same
filter. If they cannot (e.g. a physically-motivated exclusion applies to only
one condition), say so in the results text — it constrains the comparison.

## Rule 5 — cross-condition component claims

To claim that a treatment *induces a component absent in control*, fitting each
condition separately and comparing component counts is not sufficient. Fit the
constrained model in which the shared component's parameters are tied to the
control values, and test the added component against that constrained null with
the Rule 1 bootstrap.

## Reporting template

Always report, for each fit: n, number of components compared, initialization
and number of restarts, convergence status, observed Δ criterion, bootstrap
null p-value with number of simulations, and the resulting claim in one
sentence with its limits.
