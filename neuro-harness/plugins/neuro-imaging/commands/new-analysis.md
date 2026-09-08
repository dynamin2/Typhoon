---
description: Start a provenance-tracked analysis run on an imaging dataset
---

Start a new analysis run.

1. Ask which input data file(s) to use if not given. Never assume a path.
2. Create a timestamped output directory `results/<UTC timestamp>/`.
3. Before any fitting, report: n per group, the acquisition frame interval, how
   many records fall into each declared exclusion category, and whether the
   groups passed through identical filters. Stop and flag any asymmetry.
4. Run the requested analysis, applying `imaging-trace-analysis`, and
   `mixture-model-guardrails` for any component-count or population claim.
5. Write all figures per `lab-output-conventions` (SVG + PDF + PNG + source
   data xlsx + figure source file) and the results to a .docx.
6. Write the provenance record per `experiment-provenance`.
7. Close with one sentence stating what the data cannot distinguish.
