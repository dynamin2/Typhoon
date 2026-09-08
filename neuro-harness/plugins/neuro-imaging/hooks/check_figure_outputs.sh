#!/usr/bin/env bash
# Advisory check: warn when a recently written figure lacks the lab's full
# format set. Prints to stderr; never blocks.
set -uo pipefail
root="${CLAUDE_PROJECT_DIR:-$PWD}"
missing=0
while IFS= read -r -d '' svg; do
  base="${svg%.svg}"
  for ext in pdf png; do
    if [ ! -f "${base}.${ext}" ]; then
      echo "figure missing ${ext}: ${base}.${ext}" >&2; missing=1
    fi
  done
  if ! ls "${base}"*_sourcedata.xlsx >/dev/null 2>&1; then
    echo "figure missing source data: ${base}_sourcedata.xlsx" >&2; missing=1
  fi
done < <(find "$root" -name '*.svg' -newermt '-2 hours' -print0 2>/dev/null)
if [ "$missing" -eq 1 ]; then
  echo "lab-output-conventions not satisfied for the above." >&2
fi
exit 0
