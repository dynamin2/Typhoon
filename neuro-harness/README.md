# neuro-harness

A Claude Code harness for presynaptic / synaptic neuroscience research.

This is not a collection of prompts. It is a set of standing constraints that
apply to every session without being restated: the lab's output conventions,
statistical guardrails for population and component-count claims, trace-analysis
discipline for optical recordings, and a provenance record for every result.

## Layout

```
neuro-harness/
├── .claude-plugin/
│   └── marketplace.json           # lists the plugins in this repo
└── plugins/
    ├── neuro-imaging/
    │   ├── .claude-plugin/plugin.json
    │   ├── skills/
    │   │   ├── lab-output-conventions/     # SVG+PDF+PNG, source data, .docx
    │   │   ├── mixture-model-guardrails/   # bootstrap null before any k>1 claim
    │   │   ├── imaging-trace-analysis/     # ΔF/F₀, bleaching, latency, trains
    │   │   └── experiment-provenance/      # hashes, seeds, environment
    │   ├── commands/new-analysis.md        # /neuro-imaging:new-analysis
    │   ├── hooks/hooks.json                # advisory figure-format check
    │   ├── scripts/write_provenance.py
    │   └── .mcp.json
    └── neuro-lit/
        ├── .claude-plugin/plugin.json
        ├── skills/citation-discipline/
        └── .mcp.json
```

## Install

```bash
# 1. push this directory to your own GitHub repo first, then:
/plugin marketplace add YOURNAME/neuro-harness
/plugin install neuro-imaging
/plugin install neuro-lit
```

For local development without pushing, point the marketplace at the directory:

```bash
/plugin marketplace add /absolute/path/to/neuro-harness
```

Verify:

```bash
/plugin           # lists installed plugins
/mcp              # confirms MCP servers connected
```

Skills load automatically when their `description` matches what you are doing;
they are not invoked by name.

## Extending it

The three things worth adding next, in order of payoff:

1. **Wrap your existing pipelines.** Move `Ca_Glu_alignment_GMM_v4.py`,
   `MuTiGaUsS_AIC_v4.m`, and the iGluSnFR train scripts into
   `plugins/neuro-imaging/scripts/`, and add a skill per pipeline stating what
   it takes, what it emits, and its known limitations. This is what turns the
   harness from conventions into an actual analysis environment.
2. **A construct-design skill** encoding the lab's cloning conventions —
   Gibson overlap rules, linker choices, which luminal loop tolerates insertion
   in which SCAMP family member, AAV titer expectations.
3. **A manuscript skill** for the lab's writing conventions, if the existing
   editor skills do not already cover it.

## Caveats

- The `hooks.json` schema and the `${CLAUDE_PLUGIN_ROOT}` variable should be
  checked against the current Claude Code plugin reference before you rely on
  the hook firing; the skills, commands and MCP config do not depend on it.
- The `biomcp` MCP entry assumes `uvx` is on PATH and that the package name
  matches the current release. Run `/mcp` after install to confirm it connects.
- Nothing here validates science. The guardrails make certain errors loud; they
  do not make results correct.
