---
name: experiment-provenance
description: Use at the end of any analysis that produces a figure, statistic, or results document, and whenever re-running or revising a previous analysis. Writes a provenance record linking each output artifact to the exact input file, script version, parameters, seeds, and environment that produced it, so results can be reproduced or defended months later.
---

# Experiment provenance

Every analysis run writes a record. Without it, a figure six months old cannot
be defended in review.

## Where it goes

Alongside the outputs, as `PROVENANCE.md` plus a machine-readable
`provenance.json` in the same directory.

## What it must contain

- **Inputs**: absolute path, file size, and SHA-256 of every input data file.
  The hash is what makes "we reran it on the same data" checkable.
- **Code**: script path, git commit if under version control, SHA-256 of the
  script if not. Note explicitly if the working tree was dirty.
- **Parameters**: every non-default parameter actually used, including random
  seeds. Seeds are mandatory for anything involving bootstrap, GMM
  initialization, or resampling.
- **Environment**: Python/MATLAB version and versions of numpy, scipy,
  scikit-learn, matplotlib, or the MATLAB toolboxes used.
- **Outputs**: each produced file with its SHA-256.
- **Exclusions**: records in, records out, and the criterion. This duplicates
  what the results document says, deliberately.
- **Runtime**: timestamp and wall time.

## Rules

1. Never overwrite a previous run's outputs. Write into a timestamped
   subdirectory so earlier versions of a figure remain recoverable.
2. If an analysis is rerun with changed parameters, write a new record and note
   which record it supersedes and why.
3. If an input cannot be hashed (still open, on a network share), say so in the
   record rather than omitting the field.
4. A results claim that cannot be traced to a recorded run is not reportable.
   When asked about a number whose provenance record is missing, say the
   provenance is missing rather than recomputing and presenting the new value
   as the original.

## Helper

`scripts/write_provenance.py` collects most of this automatically; call it at
the end of an analysis script. It is a convenience, not a substitute — the
exclusion criteria and superseding notes are written by hand.
