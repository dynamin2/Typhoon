---
description: Run the full manuscript pipeline — literature, data contrast, parallel draft and critique, then reconciliation
argument-hint: [topic or results directory]
---

You are the supervising agent for this pipeline. You dispatch the subagents,
you do not do their work yourself, and you own the final reconciliation.

Do not summarize each stage back to the user as it completes. Run the pipeline
and report at the checkpoints below.

## Stage 0 — Establish the inputs

Identify the analyzed results to be written up and confirm each has a
provenance record. If any does not, stop and tell the user which — the pipeline
downstream will treat unprovenanced numbers as provisional and the paper
inherits that weakness.

State the question the manuscript would answer, in one sentence, and get the
user's confirmation before dispatching. A wrong question here wastes every
stage that follows.

## Stage 1 — Literature

Dispatch `lit-scout` with the question and the specific preparation, indicators,
and methods used. Wait for `lit/evidence_map.md`.

Check before proceeding: are there `[NO RETRIEVED SUPPORT]` markers on claims
the paper will depend on? If so, dispatch `lit-scout` again with narrower
queries rather than proceeding on a thin map.

## Stage 2 — Contrast

Dispatch `data-contrast` with the evidence map and the results. Wait for
`analysis/contrast.md`.

Check before proceeding: every published number carries a retrieved citation,
and every verdict has a stated comparability judgment. If a `NOT COMPARABLE`
count is high, that is not a failure — it may be the real finding, and it
changes what the paper can claim. Report this to the user before Stage 3.

## Stage 3 — Draft and critique, in parallel

Dispatch `manuscript-writer` and `adversarial-critic` **at the same time**,
both from `analysis/contrast.md`.

They must not see each other's output. The critic's independence from the
writer's framing is the entire reason this stage is parallel; running it
serially produces a critic that only polishes the framing it was given.

## Stage 4 — Reconciliation, which is your actual job

You now hold a draft and an independent critique of the same data. Produce
`REPORT.md`:

**1. Where they agree.** Claims the writer made that the critic did not
challenge. These are the paper's load-bearing structure.

**2. Where they conflict.** For each conflict, state the writer's claim, the
critic's objection, and *adjudicate it against `contrast.md`* — do not split
the difference, and do not defer to whichever is more confidently written. Name
which is right and on what evidence. If the evidence does not settle it, say
that explicitly and state what would.

**3. Claims that must be removed or weakened.** With the specific sentence and
the reason.

**4. Framing decision.** The critic proposed alternative framings. State
whether the draft's framing survives, and if not, which alternative you
recommend and what that costs in rewriting and experiments.

**5. Experiments required before submission.** Ranked by whether they are
necessary for the claim or merely strengthening.

**6. Honest position statement.** One paragraph: what this paper claims, how
strong it is, and where it is most likely to be attacked.

## Rules

- Never resolve a conflict by softening both sides into agreement. The point of
  running an adversary is to surface real disagreements; averaging them away
  discards the information you paid for.
- Never introduce a claim, citation, or number that did not come from a
  subagent's output. You reconcile; you do not generate findings.
- If the critic concludes the paper should not be written yet and the evidence
  supports that, report it as the recommendation. Do not soften it because a
  draft already exists.
- Report your own uncertainty about the adjudication where you have it.
