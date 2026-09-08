---
name: lab-output-conventions
description: Apply the Chang lab's standing output rules whenever producing figures, analysis results, or analysis scripts from imaging or electrophysiology data. Use for any task that generates a plot, a statistics table, or a results document — including quick exploratory plots. Covers required file formats (SVG/PDF/PNG), source-data export for GraphPad, Word results documents, and figure-editability requirements for Illustrator.
---

# Lab output conventions

These are standing requirements. Apply them without being asked. If a request
conflicts with one, say so and ask — do not silently drop it.

## Every figure

Save **all three** formats, same basename, in the same output directory:

| Format | Purpose | Requirement |
|---|---|---|
| `.svg` | Illustrator editing | Text must stay as text, **not** outlines. No rasterized panels. |
| `.pdf` | Vector archive / submission | Fonts embedded (Type 42 / TrueType). |
| `.png` | Slides, quick viewing | 300 dpi minimum. |

Also save the **figure source file** of the plotting environment itself:
`.fig` for MATLAB, pickled figure or the generating `.py` for Python. The point
is that a figure can be reopened and edited months later, not just re-rendered.

Matplotlib preamble that satisfies the above:

```python
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'    # keep text editable in Illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['figure.dpi']   = 300
```

## Every figure needs its source data

For each figure, write an Excel workbook `<figure_basename>_sourcedata.xlsx`
with one sheet per panel, containing the exact plotted values (not the raw
upstream file). Column headers must be self-describing. This is so the figure
can be rebuilt in GraphPad or given to a journal without re-running anything.

Histograms: export bin edges, counts, and the fitted curve evaluated on a dense
grid — not just the raw values.

## Analysis results go into a Word document

Numerical results, model-selection tables, and test statistics are written to a
`.docx`, not left as terminal output. Include: n per group, what was excluded
and why, the exact model compared, and the software/version used.

## Data loading

Never hard-code an input path. Analysis scripts take the Excel/data file as an
argument or a clearly-marked variable at the top of the file. The lab re-runs
these scripts on new files constantly.

## Units and axis conventions

State units on every axis. For imaging traces, mark the acquisition frame
interval on the figure or in the legend when latency or kinetics are plotted —
quantization by frame interval is a real constraint on interpretation, not a
detail.
