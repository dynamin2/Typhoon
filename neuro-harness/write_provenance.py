#!/usr/bin/env python3
"""Write a provenance record for an analysis run.

Usage (from inside an analysis script):

    from write_provenance import write_provenance
    write_provenance(
        outdir="results/2026-09-07T14-30",
        inputs=["data/N_CTL_d4.xlsx", "data/N_TTX_d.xlsx"],
        script=__file__,
        params={"n_bootstrap": 1000, "seed": 20260907, "k_max": 3},
        outputs=["fig2_latency.svg", "fig2_latency.pdf", "fig2_latency.png"],
        exclusions={"CTL": {"in": 102, "out": 49, "criterion": "positive latency only"}},
    )

Or from the shell to hash a set of files:

    python write_provenance.py file1.xlsx file2.xlsx
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: str | Path, _chunk: int = 1 << 20) -> str | None:
    p = Path(path)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    try:
        with p.open("rb") as fh:
            for block in iter(lambda: fh.read(_chunk), b""):
                h.update(block)
    except OSError:
        return None
    return h.hexdigest()


def _file_record(path: str | Path) -> dict:
    p = Path(path)
    digest = sha256(p)
    return {
        "path": str(p.resolve()) if p.exists() else str(p),
        "exists": p.exists(),
        "size_bytes": p.stat().st_size if p.is_file() else None,
        "sha256": digest,
        "note": None if digest else "could not hash (missing, open, or unreadable)",
    }


def _git_state(script: str | Path) -> dict:
    d = Path(script).resolve().parent
    def run(*args):
        try:
            return subprocess.run(
                args, cwd=d, capture_output=True, text=True, timeout=10
            ).stdout.strip() or None
        except Exception:
            return None
    commit = run("git", "rev-parse", "HEAD")
    if not commit:
        return {"under_version_control": False}
    return {
        "under_version_control": True,
        "commit": commit,
        "branch": run("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(run("git", "status", "--porcelain")),
    }


def _packages() -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("numpy", "scipy", "sklearn", "matplotlib", "pandas", "openpyxl"):
        try:
            out[name] = __import__(name).__version__
        except Exception:
            pass
    return out


def write_provenance(
    outdir: str | Path,
    inputs: list,
    script: str,
    params: dict,
    outputs: list | None = None,
    exclusions: dict | None = None,
    supersedes: str | None = None,
    notes: str | None = None,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if "seed" not in params:
        print(
            "WARNING: no random seed recorded. Any bootstrap, GMM "
            "initialization, or resampling in this run is not reproducible.",
            file=sys.stderr,
        )

    record = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "working_directory": os.getcwd(),
        "script": _file_record(script) | {"git": _git_state(script)},
        "inputs": [_file_record(p) for p in inputs],
        "parameters": params,
        "environment": _packages(),
        "outputs": [_file_record(outdir / p) for p in (outputs or [])],
        "exclusions": exclusions or {},
        "supersedes": supersedes,
        "notes": notes,
    }

    (outdir / "provenance.json").write_text(json.dumps(record, indent=2))

    lines = [f"# Provenance\n", f"Run: {record['timestamp_utc']}\n"]
    lines.append(f"## Script\n\n`{record['script']['path']}`\n")
    g = record["script"]["git"]
    if g.get("under_version_control"):
        dirty = "  **(uncommitted changes present)**" if g["dirty"] else ""
        lines.append(f"commit `{g['commit'][:12]}` on `{g['branch']}`{dirty}\n")
    else:
        lines.append(f"not under version control; sha256 `{record['script']['sha256']}`\n")

    lines.append("\n## Inputs\n")
    for f in record["inputs"]:
        lines.append(f"- `{f['path']}` — sha256 `{f['sha256']}`\n")

    lines.append("\n## Parameters\n\n```json\n")
    lines.append(json.dumps(params, indent=2) + "\n```\n")

    if exclusions:
        lines.append("\n## Exclusions\n")
        for grp, e in exclusions.items():
            lines.append(
                f"- **{grp}**: {e.get('in')} in → {e.get('out')} out "
                f"({e.get('criterion')})\n"
            )

    if record["outputs"]:
        lines.append("\n## Outputs\n")
        for f in record["outputs"]:
            lines.append(f"- `{Path(f['path']).name}` — sha256 `{f['sha256']}`\n")

    if supersedes:
        lines.append(f"\n## Supersedes\n\n{supersedes}\n")
    if notes:
        lines.append(f"\n## Notes\n\n{notes}\n")

    (outdir / "PROVENANCE.md").write_text("".join(lines))
    return outdir / "PROVENANCE.md"


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(json.dumps(_file_record(arg), indent=2))
