#!/usr/bin/env python3
"""
AAVLINK CDS Split-Site Finder v2.0

Based on: Lin et al. (2026) Cell 189, 969–986.
"AAVLINK: A potent DNA-recombination method for large cargo delivery in gene therapy"

Design principles extracted from the paper:
1. Split junction MUST conform to exonic AG/G consensus (Sibley et al., 2016, Nat Rev Genet)
   - 5' fragment ends with ...AG
   - 3' fragment starts with G...
   - This enables proper RNA splicing after Cre-mediated DNA recombination

2. Fragment sizes must fit within AAV packaging capacity (~4.7 kb per vector)
   - 5' CDS vector overhead: ~550 bp (ITRs + promoter + SD + lox)
   - 3' CDS vector w/ Cre (Dual, Triple II): ~2300 bp (ITRs + promoter + Cre + SD + lox + SA + polyA)
   - 3' CDS vector w/o Cre (Triple I): ~550 bp
   - Body CDS vector (Triple II): ~550 bp

3. System type selection:
   - Dual:      CDS < ~6.2 kb (gene split into 2 vectors; Cre in 3' vector)
   - Triple I:  CDS 6.2–8.0 kb (gene split into 2 vectors; Cre in separate 3rd vector)
   - Triple II: CDS > 8.0 kb (gene split into 3 vectors; Cre in 3' vector)

4. Codon-phase preference at split junction:
   - Phase 2 (AG are 2nd–3rd nt of preceding codon): 47% of validated sites
   - Phase 0 (split at codon boundary): 35%
   - Phase 1 (AG straddles two codons): 18%

5. Synonymous mutations can create AG|G junctions without changing protein sequence

6. More balanced fragment lengths are preferred (max/min ratio ≤ 3 required, ≤ 2 preferred)

Features:
- CLI mode (headless) and GUI mode (tkinter)
- Automatic system type recommendation
- Composite scoring with weighted criteria
- Synonymous mutation suggestions for additional candidate sites
- AAV packaging feasibility check
- Export results as text or CSV

Usage:
  GUI:  python AAVLINK_SplitFinder.py
  CLI:  python AAVLINK_SplitFinder.py --cli --seq ATGCCC... --cds_start 1 --cds_end 5418
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import math
import os
import re
import sys
import textwrap
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

# ─────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────

CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Build reverse table: amino acid -> set of codons
AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# AAV vector overhead estimates (bp) — derived from validated constructs in Table S1
AAV_CAPACITY = 4700  # bp, total AAV packaging capacity

# Overhead for each vector type (ITRs + regulatory elements)
OVERHEAD_5CDS = 550        # 5' CDS vector: ITR(145×2) + Promoter(~200) + SD(~30) + lox(34)
OVERHEAD_3CDS_WITH_CRE = 2300  # 3' CDS w/ Cre: ITR(145×2) + Promoter(~200) + Cre(1032) + SD + lox + SA + polyA
OVERHEAD_3CDS_NO_CRE = 550     # 3' CDS w/o Cre (Triple I): ITR(145×2) + SA + polyA + lox
OVERHEAD_BODY_CDS = 550        # Body CDS (Triple II): ITR(145×2) + SA + SD + lox×2

MAX_CDS_5 = AAV_CAPACITY - OVERHEAD_5CDS                  # ~4150 bp
MAX_CDS_3_WITH_CRE = AAV_CAPACITY - OVERHEAD_3CDS_WITH_CRE  # ~2400 bp
MAX_CDS_3_NO_CRE = AAV_CAPACITY - OVERHEAD_3CDS_NO_CRE      # ~4150 bp
MAX_CDS_BODY = AAV_CAPACITY - OVERHEAD_BODY_CDS              # ~4150 bp

# System type thresholds (bp of CDS)
DUAL_MAX_CDS = 6200
TRIPLE_I_MAX_CDS = 8000

# Scoring weights
W_BALANCE = 0.30       # Fragment length balance
W_PACKAGING = 0.25     # AAV packaging feasibility
W_PHASE = 0.20         # Codon phase preference
W_NATURAL = 0.15       # Natural AG|G (vs synonymous mutation)
W_MIDPOINT = 0.10      # Proximity to ideal split position

# Codon phase scores (higher = more preferred, based on empirical data from Table S1)
PHASE_SCORE = {0: 0.70, 1: 0.36, 2: 1.00}


# ─────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────

@dataclass
class SplitCandidate:
    """A single candidate split position."""
    c_pos: int                    # CDS position (1-based); 3' fragment starts here
    is_natural: bool              # True if native AG|G; False if requires synonymous mutation
    synonymous_mutations: List[str] = field(default_factory=list)  # e.g., ["c.2200 A→G (Arg→Arg)"]
    codon_phase: int = 0          # 0, 1, or 2
    left_codon: str = ""          # codon straddling the junction (5' side)
    right_codon: str = ""         # codon straddling the junction (3' side)
    left_aa: str = ""
    right_aa: str = ""
    nt_context: str = ""          # flanking nt sequence with | at junction
    aa_context: str = ""          # flanking aa sequence


@dataclass
class SplitCombination:
    """A valid combination of split points for a given fragment count."""
    candidates: Tuple[SplitCandidate, ...]
    fragment_lengths: List[int]
    system_type: str              # "Dual", "Triple I", "Triple II"
    score: float = 0.0
    score_details: Dict[str, float] = field(default_factory=dict)
    packaging_ok: bool = True


# ─────────────────────────────────────────────────────
# Core sequence utilities
# ─────────────────────────────────────────────────────

def clean_sequence(seq: str) -> str:
    """Remove whitespace and convert to uppercase."""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def validate_dna(seq: str) -> None:
    if not seq:
        raise ValueError("Sequence is empty.")
    invalid = sorted(set(re.sub(r"[ACGTN]", "", seq)))
    if invalid:
        raise ValueError(f"Invalid characters found: {', '.join(invalid)}")


def translate_cds(cds: str) -> str:
    """Translate CDS to amino acid sequence."""
    aa = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def get_codon_at(cds: str, pos_1b: int) -> Tuple[str, int, int]:
    """Return (codon, codon_start_1b, phase) for a given CDS position."""
    pos_0b = pos_1b - 1
    codon_start_0b = (pos_0b // 3) * 3
    phase = pos_0b % 3
    codon = cds[codon_start_0b:codon_start_0b + 3]
    return codon, codon_start_0b + 1, phase


# ─────────────────────────────────────────────────────
# AG|G junction detection
# ─────────────────────────────────────────────────────

def find_natural_agg_sites(cds: str, search_start: int, search_end: int) -> List[SplitCandidate]:
    """
    Find positions where the native CDS has AG|G junction.

    c.n means: 5' fragment = cds[0:n-1], 3' fragment = cds[n-1:]
    Condition: cds[n-3:n-1] == "AG" and cds[n-1] == "G"
    """
    candidates = []
    cds_len = len(cds)
    search_start = max(3, search_start)
    search_end = min(cds_len, search_end)

    for n in range(search_start, search_end + 1):
        if n < 3 or n > cds_len:
            continue
        left_2nt = cds[n - 3:n - 1]   # last 2 nt of 5' fragment
        right_1nt = cds[n - 1]         # first nt of 3' fragment
        if left_2nt == "AG" and right_1nt == "G":
            cand = _build_candidate(cds, n, is_natural=True)
            candidates.append(cand)

    return candidates


def find_synonymous_agg_sites(cds: str, search_start: int, search_end: int) -> List[SplitCandidate]:
    """
    Find positions where synonymous mutations could create AG|G junction.

    For each position n (3' fragment start), we need:
      cds[n-3] == A, cds[n-2] == G, cds[n-1] == G

    If any of these don't match, check if a synonymous mutation at that
    position can change the nucleotide without altering the amino acid.
    """
    candidates = []
    cds_len = len(cds)
    search_start = max(3, search_start)
    search_end = min(cds_len, search_end)

    for n in range(search_start, search_end + 1):
        if n < 3 or n > cds_len:
            continue

        # Already a natural site? Skip
        if cds[n - 3:n - 1] == "AG" and cds[n - 1] == "G":
            continue

        # Positions to check (0-based): n-3, n-2, n-1
        # Required nucleotides:          A,   G,   G
        target_nts = [("A", n - 3), ("G", n - 2), ("G", n - 1)]

        mutations_needed = []
        feasible = True

        for required_nt, pos_0b in target_nts:
            if pos_0b < 0 or pos_0b >= cds_len:
                feasible = False
                break
            current_nt = cds[pos_0b]
            if current_nt == required_nt:
                continue

            # Check if synonymous mutation is possible
            codon_start_0b = (pos_0b // 3) * 3
            if codon_start_0b + 3 > cds_len:
                feasible = False
                break

            original_codon = cds[codon_start_0b:codon_start_0b + 3]
            original_aa = CODON_TABLE.get(original_codon, "X")
            if original_aa in ("*", "X"):
                feasible = False
                break

            pos_in_codon = pos_0b - codon_start_0b
            mutant_codon = list(original_codon)
            mutant_codon[pos_in_codon] = required_nt
            mutant_codon = "".join(mutant_codon)
            mutant_aa = CODON_TABLE.get(mutant_codon, "X")

            if mutant_aa == original_aa:
                c_pos = pos_0b + 1  # 1-based CDS position
                mutations_needed.append(
                    f"c.{c_pos} {current_nt}→{required_nt} "
                    f"({original_codon}→{mutant_codon}, {original_aa} silent)"
                )
            else:
                feasible = False
                break

        if feasible and mutations_needed:
            cand = _build_candidate(cds, n, is_natural=False)
            cand.synonymous_mutations = mutations_needed
            candidates.append(cand)

    return candidates


def _build_candidate(cds: str, n: int, is_natural: bool) -> SplitCandidate:
    """Build a SplitCandidate with context information."""
    cds_len = len(cds)

    # Codon phase: the G at position n is at phase (n-1) % 3
    phase = (n - 1) % 3

    # Flanking nucleotide context (20 nt each side)
    flank = 20
    left_start = max(0, n - 1 - flank)
    right_end = min(cds_len, n - 1 + flank)
    left_seq = cds[left_start:n - 1]
    right_seq = cds[n - 1:right_end]
    nt_context = f"{left_seq.lower()}|{right_seq.lower()}"
    # Highlight the AG|G junction in uppercase
    if len(left_seq) >= 2:
        nt_context = (
            left_seq[:-2].lower()
            + left_seq[-2:].upper()
            + "|"
            + right_seq[0:1].upper()
            + right_seq[1:].lower()
        )

    # Flanking amino acid context
    aa_context = _get_aa_context(cds, n, flank_nt=20)

    # Codons at junction
    left_codon = ""
    right_codon = ""
    left_aa = ""
    right_aa = ""
    if n >= 4:
        lc, _, _ = get_codon_at(cds, n - 1)
        left_codon = lc
        left_aa = CODON_TABLE.get(lc, "X")
    if n <= cds_len:
        rc, _, _ = get_codon_at(cds, n)
        right_codon = rc
        right_aa = CODON_TABLE.get(rc, "X")

    return SplitCandidate(
        c_pos=n,
        is_natural=is_natural,
        codon_phase=phase,
        left_codon=left_codon,
        right_codon=right_codon,
        left_aa=left_aa,
        right_aa=right_aa,
        nt_context=nt_context,
        aa_context=aa_context,
    )


def _get_aa_context(cds: str, split_start: int, flank_nt: int = 20) -> str:
    """Generate amino acid context string around split point."""
    split0 = split_start - 1
    start_codon_idx = max(0, (max(0, split0 - flank_nt) // 3))
    end_codon_idx = min(len(cds) // 3, ((split0 + flank_nt) // 3) + 1)

    split_codon_idx = split0 // 3
    phase_in_codon = split0 % 3

    aa_list = []
    for ci in range(start_codon_idx, end_codon_idx):
        codon = cds[ci * 3:(ci + 1) * 3]
        aa = CODON_TABLE.get(codon, "X")
        if ci == split_codon_idx:
            if phase_in_codon == 0:
                aa_list.append(f"|{aa}")
            else:
                aa_list.append(f"{aa}(split@{phase_in_codon + 1}/3)")
        else:
            aa_list.append(aa)

    return "".join(aa_list)


# ─────────────────────────────────────────────────────
# System type recommendation
# ─────────────────────────────────────────────────────

def recommend_system(cds_len: int) -> str:
    """Recommend AAVLINK system type based on CDS length."""
    if cds_len <= DUAL_MAX_CDS:
        return "Dual"
    elif cds_len <= TRIPLE_I_MAX_CDS:
        return "Triple I"
    else:
        return "Triple II"


def get_fragment_count(system_type: str) -> int:
    """Return required number of gene fragments (split points = fragments - 1)."""
    if system_type == "Dual":
        return 2
    elif system_type == "Triple I":
        return 2   # Gene in 2 parts; Cre in separate vector
    else:  # Triple II
        return 3   # Gene in 3 parts


def get_max_fragment_sizes(system_type: str) -> List[int]:
    """Return max CDS size per fragment slot, ordered [5', (body,) 3']."""
    if system_type == "Dual":
        return [MAX_CDS_5, MAX_CDS_3_WITH_CRE]
    elif system_type == "Triple I":
        return [MAX_CDS_5, MAX_CDS_3_NO_CRE]
    else:  # Triple II
        return [MAX_CDS_5, MAX_CDS_BODY, MAX_CDS_3_WITH_CRE]


# ─────────────────────────────────────────────────────
# Scoring
# ─────────────────────────────────────────────────────

def score_combination(
    combo: SplitCombination,
    cds_len: int,
) -> float:
    """
    Compute a composite score (0–1, higher = better) for a split combination.
    """
    frag_lens = combo.fragment_lengths
    max_sizes = get_max_fragment_sizes(combo.system_type)

    # --- Balance score ---
    ratio = max(frag_lens) / max(1, min(frag_lens))
    if ratio <= 1.0:
        balance_score = 1.0
    elif ratio <= 2.0:
        balance_score = 1.0 - 0.5 * (ratio - 1.0)
    elif ratio <= 3.0:
        balance_score = 0.5 - 0.5 * (ratio - 2.0)
    else:
        balance_score = 0.0

    # --- Packaging score ---
    pkg_score = 1.0
    all_fit = True
    for i, (frag, maxs) in enumerate(zip(frag_lens, max_sizes)):
        if frag > maxs:
            overshoot = (frag - maxs) / maxs
            pkg_score -= min(0.5, overshoot)
            all_fit = False
        else:
            margin = (maxs - frag) / maxs
            pkg_score_i = min(1.0, 0.5 + margin)  # bonus for headroom
    pkg_score = max(0.0, pkg_score)
    combo.packaging_ok = all_fit

    # --- Phase score (average across split points) ---
    phase_scores = [PHASE_SCORE.get(c.codon_phase, 0.5) for c in combo.candidates]
    phase_score = sum(phase_scores) / len(phase_scores) if phase_scores else 0.5

    # --- Natural AG|G score ---
    natural_count = sum(1 for c in combo.candidates if c.is_natural)
    natural_score = natural_count / len(combo.candidates) if combo.candidates else 0.0

    # --- Midpoint proximity score ---
    # Ideal split: equal-size fragments
    n_frags = len(frag_lens)
    ideal_size = cds_len / n_frags
    deviations = [abs(f - ideal_size) / cds_len for f in frag_lens]
    midpoint_score = max(0.0, 1.0 - 2.0 * max(deviations))

    # --- Composite ---
    total = (
        W_BALANCE * balance_score
        + W_PACKAGING * pkg_score
        + W_PHASE * phase_score
        + W_NATURAL * natural_score
        + W_MIDPOINT * midpoint_score
    )

    combo.score = round(total, 4)
    combo.score_details = {
        "balance": round(balance_score, 3),
        "packaging": round(pkg_score, 3),
        "phase": round(phase_score, 3),
        "natural_AGG": round(natural_score, 3),
        "midpoint": round(midpoint_score, 3),
    }
    return total


# ─────────────────────────────────────────────────────
# Fragment length calculation
# ─────────────────────────────────────────────────────

def compute_fragment_lengths(cds_len: int, split_positions: Tuple[int, ...]) -> List[int]:
    """
    Given CDS length and sorted split positions (c.n = 3' fragment start),
    compute fragment lengths.
    """
    lengths = []
    prev = 1
    for sp in split_positions:
        lengths.append(sp - prev)
        prev = sp
    lengths.append(cds_len - prev + 1)
    return lengths


# ─────────────────────────────────────────────────────
# Main analysis pipeline
# ─────────────────────────────────────────────────────

def analyze_cds(
    cds: str,
    system_type: Optional[str] = None,
    search_start: int = 1,
    search_end: Optional[int] = None,
    include_synonymous: bool = True,
    max_ratio: float = 3.0,
    top_n: int = 20,
) -> Dict:
    """
    Main analysis function.

    Parameters
    ----------
    cds : str
        CDS sequence (ATG...stop).
    system_type : str, optional
        "Dual", "Triple I", or "Triple II". Auto-detected if None.
    search_start, search_end : int
        CDS coordinate range to search (1-based).
    include_synonymous : bool
        Whether to include synonymous-mutation-derived AG|G sites.
    max_ratio : float
        Maximum allowed max/min fragment length ratio.
    top_n : int
        Number of top-scoring combinations to return.

    Returns
    -------
    dict with keys: system_type, candidates, combinations, summary
    """
    cds_len = len(cds)
    if search_end is None:
        search_end = cds_len

    # Auto-detect system type
    if system_type is None:
        system_type = recommend_system(cds_len)

    frag_count = get_fragment_count(system_type)
    num_splits = frag_count - 1

    # Find candidate split sites
    natural_sites = find_natural_agg_sites(cds, search_start, search_end)
    syn_sites = []
    if include_synonymous:
        syn_sites = find_synonymous_agg_sites(cds, search_start, search_end)

    all_candidates = natural_sites + syn_sites
    all_candidates.sort(key=lambda c: c.c_pos)

    # Remove duplicates (same position)
    seen = set()
    unique_candidates = []
    for c in all_candidates:
        if c.c_pos not in seen:
            seen.add(c.c_pos)
            unique_candidates.append(c)

    # Build lookup
    cand_map = {c.c_pos: c for c in unique_candidates}
    cand_positions = [c.c_pos for c in unique_candidates]

    # Enumerate valid combinations
    valid_combos: List[SplitCombination] = []

    for combo_positions in itertools.combinations(cand_positions, num_splits):
        frag_lens = compute_fragment_lengths(cds_len, combo_positions)

        # Check minimum fragment size (avoid trivially small fragments)
        if min(frag_lens) < 100:
            continue

        # Check ratio
        ratio = max(frag_lens) / min(frag_lens)
        if ratio > max_ratio:
            continue

        cands_tuple = tuple(cand_map[p] for p in combo_positions)
        sc = SplitCombination(
            candidates=cands_tuple,
            fragment_lengths=frag_lens,
            system_type=system_type,
        )
        score_combination(sc, cds_len)
        valid_combos.append(sc)

    # Sort by score descending
    valid_combos.sort(key=lambda x: x.score, reverse=True)
    top_combos = valid_combos[:top_n]

    # Collect all unique valid split points
    valid_points_set: Set[int] = set()
    for sc in valid_combos:
        for c in sc.candidates:
            valid_points_set.add(c.c_pos)

    summary = {
        "cds_length": cds_len,
        "system_type": system_type,
        "fragment_count": frag_count,
        "search_range": (search_start, search_end),
        "natural_AG|G_sites": len(natural_sites),
        "synonymous_AG|G_sites": len(syn_sites),
        "total_candidates": len(unique_candidates),
        "valid_combinations": len(valid_combos),
        "max_fragment_sizes": get_max_fragment_sizes(system_type),
    }

    return {
        "system_type": system_type,
        "candidates": unique_candidates,
        "combinations": top_combos,
        "all_valid_count": len(valid_combos),
        "summary": summary,
    }


# ─────────────────────────────────────────────────────
# Formatting / reporting
# ─────────────────────────────────────────────────────

def format_report(result: Dict, cds_start_genomic: int = 1) -> str:
    """Format analysis results as a human-readable text report."""
    lines = []
    s = result["summary"]

    lines.append("=" * 100)
    lines.append("  AAVLINK Split-Site Finder v2.0 — Analysis Report")
    lines.append("  Based on Lin et al. (2026) Cell 189, 969–986")
    lines.append("=" * 100)
    lines.append("")

    # Summary
    lines.append("[INPUT SUMMARY]")
    lines.append(f"  CDS length:           {s['cds_length']} bp")
    lines.append(f"  Recommended system:   {s['system_type']} "
                 f"({'gene in ' + str(s['fragment_count']) + ' fragment(s)'})")
    lines.append(f"  Search range (CDS):   c.{s['search_range'][0]}–c.{s['search_range'][1]}")
    lines.append(f"  Max CDS per vector:   {', '.join(str(x) + ' bp' for x in s['max_fragment_sizes'])}")
    lines.append("")

    lines.append("[CANDIDATE SPLIT SITES]")
    lines.append(f"  Natural AG|G sites:    {s['natural_AG|G_sites']}")
    lines.append(f"  Synonymous AG|G sites: {s['synonymous_AG|G_sites']}")
    lines.append(f"  Total candidates:      {s['total_candidates']}")
    lines.append(f"  Valid combinations:    {result['all_valid_count']}")
    lines.append("")

    # ASCII map
    if result["candidates"]:
        lines.append("[SPLIT SITE MAP]")
        lines.append(_make_ascii_map(s["cds_length"], result["candidates"]))
        lines.append("")

    # Top combinations
    combos = result["combinations"]
    if not combos:
        lines.append("  No valid split combinations found.")
        lines.append("  Try expanding the search range, enabling synonymous mutations,")
        lines.append("  or relaxing the max/min ratio constraint.")
        return "\n".join(lines)

    lines.append(f"[TOP {len(combos)} COMBINATIONS (ranked by composite score)]")
    lines.append("")

    for rank, combo in enumerate(combos, 1):
        pkg_flag = "✓" if combo.packaging_ok else "⚠ EXCEEDS AAV CAPACITY"
        lines.append("─" * 100)
        split_labels = ", ".join(f"c.{c.c_pos}" for c in combo.candidates)
        lines.append(f"  Rank #{rank}   Score: {combo.score:.3f}   "
                     f"Split: {split_labels}   Packaging: {pkg_flag}")

        frag_strs = []
        for i, fl in enumerate(combo.fragment_lengths):
            max_sz = get_max_fragment_sizes(combo.system_type)
            if i < len(max_sz):
                pct = fl / max_sz[i] * 100
                frag_strs.append(f"{fl} bp ({pct:.0f}% of max)")
            else:
                frag_strs.append(f"{fl} bp")
        lines.append(f"  Fragments: {' | '.join(frag_strs)}")

        # Score breakdown
        sd = combo.score_details
        lines.append(f"  Scores: balance={sd['balance']:.2f}  packaging={sd['packaging']:.2f}  "
                     f"phase={sd['phase']:.2f}  natural={sd['natural_AGG']:.2f}  "
                     f"midpoint={sd['midpoint']:.2f}")
        lines.append("")

        for j, cand in enumerate(combo.candidates, 1):
            genomic_pos = cds_start_genomic + cand.c_pos - 1
            nat_label = "natural" if cand.is_natural else "synonymous mutation"
            phase_label = f"phase {cand.codon_phase}"
            lines.append(f"  [Split {j}]  c.{cand.c_pos}  "
                         f"(genomic: {genomic_pos})  "
                         f"[{nat_label}, {phase_label}]")
            lines.append(f"    nt:  ...{cand.nt_context}...")
            lines.append(f"    aa:  {cand.aa_context}")
            if cand.synonymous_mutations:
                for mut in cand.synonymous_mutations:
                    lines.append(f"    mutation: {mut}")
            lines.append("")

    lines.append("=" * 100)
    lines.append("Notes:")
    lines.append("  - Phase 2 = AG are nt 2–3 of preceding codon (most common in validated genes)")
    lines.append("  - Phase 0 = split at codon boundary")
    lines.append("  - Phase 1 = AG straddles two codons")
    lines.append("  - Packaging '✓' means all fragments fit within AAV capacity estimates")
    lines.append("  - Score weights: balance=0.30, packaging=0.25, phase=0.20, "
                 "natural=0.15, midpoint=0.10")
    lines.append("=" * 100)

    return "\n".join(lines)


def _make_ascii_map(cds_len: int, candidates: List[SplitCandidate], width: int = 90) -> str:
    """Generate an ASCII map of candidate positions along the CDS."""
    bar = ["-"] * width
    markers = [" "] * width

    nat_positions = []
    syn_positions = []

    for c in candidates:
        pos = min(width - 1, max(0, round((c.c_pos - 1) / max(1, cds_len - 1) * (width - 1))))
        if c.is_natural:
            bar[pos] = "|"
            markers[pos] = "↑"
            nat_positions.append(c.c_pos)
        else:
            if bar[pos] != "|":
                bar[pos] = ":"
            if markers[pos] != "↑":
                markers[pos] = "·"
            syn_positions.append(c.c_pos)

    parts = []
    parts.append(f"  CDS: {cds_len} bp   | = natural AG|G   : = synonymous mutation")
    parts.append(f"  5' {''.join(bar)} 3'")
    parts.append(f"     {''.join(markers)}")
    parts.append(f"  Natural sites ({len(nat_positions)}): "
                 + ", ".join(f"c.{p}" for p in nat_positions[:30])
                 + ("..." if len(nat_positions) > 30 else ""))
    if syn_positions:
        parts.append(f"  Synonymous sites ({len(syn_positions)}): "
                     + ", ".join(f"c.{p}" for p in syn_positions[:20])
                     + ("..." if len(syn_positions) > 20 else ""))

    return "\n".join(parts)


def format_csv(result: Dict, cds_start_genomic: int = 1) -> str:
    """Format top combinations as CSV."""
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow([
        "Rank", "Score", "System", "Packaging_OK",
        "Split_Points", "Fragment_Lengths",
        "Balance", "Packaging", "Phase", "Natural", "Midpoint",
        "Mutations_Needed"
    ])

    for rank, combo in enumerate(result["combinations"], 1):
        splits = "; ".join(f"c.{c.c_pos}" for c in combo.candidates)
        frags = "; ".join(str(f) for f in combo.fragment_lengths)
        muts = "; ".join(
            m for c in combo.candidates for m in c.synonymous_mutations
        ) or "None"
        sd = combo.score_details
        writer.writerow([
            rank, combo.score, combo.system_type, combo.packaging_ok,
            splits, frags,
            sd["balance"], sd["packaging"], sd["phase"], sd["natural_AGG"], sd["midpoint"],
            muts
        ])

    return output.getvalue()


# ─────────────────────────────────────────────────────
# GUI (tkinter)
# ─────────────────────────────────────────────────────

def run_gui():
    """Launch the tkinter GUI."""
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    class App:
        def __init__(self, root: tk.Tk) -> None:
            self.root = root
            self.root.title("AAVLINK Split-Site Finder v2.0")
            self.root.geometry("1100x820")
            self.root.minsize(900, 700)
            self._last_result = None
            self._build_ui()

        def _build_ui(self) -> None:
            main = ttk.Frame(self.root, padding=8)
            main.pack(fill="both", expand=True)

            outer = ttk.Panedwindow(main, orient="vertical")
            outer.pack(fill="both", expand=True)

            top = ttk.Frame(outer)
            bottom = ttk.Frame(outer)
            outer.add(top, weight=3)
            outer.add(bottom, weight=5)

            self._build_input(top)
            self._build_output(bottom)

        def _build_input(self, parent: ttk.Frame) -> None:
            hp = ttk.Panedwindow(parent, orient="horizontal")
            hp.pack(fill="both", expand=True)

            inp = ttk.LabelFrame(hp, text="Input", padding=10)
            guide = ttk.LabelFrame(hp, text="Guide", padding=10)
            hp.add(inp, weight=5)
            hp.add(guide, weight=3)

            # Sequence input
            ttk.Label(inp, text="DNA Sequence (CDS or full genomic)").grid(
                row=0, column=0, columnspan=6, sticky="w")
            self.seq_text = scrolledtext.ScrolledText(inp, width=80, height=8, wrap="word")
            self.seq_text.grid(row=1, column=0, columnspan=6, sticky="nsew", pady=(2, 8))

            # Parameters
            r = 2
            ttk.Label(inp, text="CDS start").grid(row=r, column=0, sticky="w")
            self.cds_start_var = tk.StringVar(value="1")
            ttk.Entry(inp, textvariable=self.cds_start_var, width=10).grid(row=r, column=1, padx=4)

            ttk.Label(inp, text="CDS end").grid(row=r, column=2, sticky="w")
            self.cds_end_var = tk.StringVar()
            ttk.Entry(inp, textvariable=self.cds_end_var, width=10).grid(row=r, column=3, padx=4)

            r = 3
            ttk.Label(inp, text="Search start (CDS)").grid(row=r, column=0, sticky="w", pady=(6, 0))
            self.search_start_var = tk.StringVar()
            ttk.Entry(inp, textvariable=self.search_start_var, width=10).grid(
                row=r, column=1, padx=4, pady=(6, 0))

            ttk.Label(inp, text="Search end (CDS)").grid(row=r, column=2, sticky="w", pady=(6, 0))
            self.search_end_var = tk.StringVar()
            ttk.Entry(inp, textvariable=self.search_end_var, width=10).grid(
                row=r, column=3, padx=4, pady=(6, 0))

            r = 4
            ttk.Label(inp, text="System type").grid(row=r, column=0, sticky="w", pady=(6, 0))
            self.system_var = tk.StringVar(value="Auto")
            ttk.Combobox(
                inp, textvariable=self.system_var,
                values=["Auto", "Dual", "Triple I", "Triple II"],
                width=10, state="readonly"
            ).grid(row=r, column=1, padx=4, pady=(6, 0))

            self.syn_var = tk.BooleanVar(value=True)
            ttk.Checkbutton(inp, text="Include synonymous mutations",
                           variable=self.syn_var).grid(
                row=r, column=2, columnspan=2, sticky="w", pady=(6, 0))

            r = 5
            btn = ttk.Frame(inp)
            btn.grid(row=r, column=0, columnspan=6, sticky="w", pady=(10, 0))
            ttk.Button(btn, text="Analyze", command=self.analyze).pack(side="left", padx=(0, 6))
            ttk.Button(btn, text="Example (SCN1A)", command=self.load_example).pack(side="left", padx=(0, 6))
            ttk.Button(btn, text="Clear", command=self.clear_all).pack(side="left", padx=(0, 6))
            ttk.Button(btn, text="Save TXT", command=self.save_txt).pack(side="left", padx=(0, 6))
            ttk.Button(btn, text="Save CSV", command=self.save_csv).pack(side="left")

            inp.columnconfigure(5, weight=1)
            inp.rowconfigure(1, weight=1)

            # Guide text
            guide_content = (
                "AAVLINK Split-Site Finder v2.0\n"
                "Based on Lin et al. (2026) Cell\n\n"
                "Split junction: ...AG|G...\n"
                "  5' CDS ends with AG\n"
                "  3' CDS starts with G\n\n"
                "System types:\n"
                "  Dual:      CDS ≤ 6.2 kb\n"
                "  Triple I:  6.2–8.0 kb\n"
                "  Triple II: > 8.0 kb\n\n"
                "Scoring criteria:\n"
                "  Fragment balance (30%)\n"
                "  AAV packaging fit (25%)\n"
                "  Codon phase (20%)\n"
                "  Natural AG|G (15%)\n"
                "  Midpoint proximity (10%)\n\n"
                "Codon phases:\n"
                "  Phase 2 (most common): AG are\n"
                "    nt 2–3 of preceding codon\n"
                "  Phase 0: split at codon boundary\n"
                "  Phase 1: AG straddles codons"
            )
            gt = scrolledtext.ScrolledText(guide, width=32, height=10, wrap="word")
            gt.pack(fill="both", expand=True)
            gt.insert("1.0", guide_content)
            gt.configure(state="disabled", font=("Arial", 10))

        def _build_output(self, parent: ttk.Frame) -> None:
            self.output_text = scrolledtext.ScrolledText(
                parent, width=120, height=24, wrap="none",
                font=("Consolas", 10))
            self.output_text.pack(fill="both", expand=True, padx=4, pady=4)

        def _read_params(self):
            full_seq = clean_sequence(self.seq_text.get("1.0", "end"))
            validate_dna(full_seq)

            cds_start = int(self.cds_start_var.get() or "1")
            cds_end_str = self.cds_end_var.get()
            cds_end = int(cds_end_str) if cds_end_str else len(full_seq)

            cds = full_seq[cds_start - 1:cds_end]
            if len(cds) < 100:
                raise ValueError("CDS too short (< 100 bp).")

            ss = self.search_start_var.get()
            se = self.search_end_var.get()
            search_start = int(ss) if ss else 1
            search_end = int(se) if se else len(cds)

            sys_type = self.system_var.get()
            if sys_type == "Auto":
                sys_type = None

            return cds, cds_start, sys_type, search_start, search_end

        def analyze(self) -> None:
            try:
                cds, cds_start, sys_type, ss, se = self._read_params()
                result = analyze_cds(
                    cds=cds,
                    system_type=sys_type,
                    search_start=ss,
                    search_end=se,
                    include_synonymous=self.syn_var.get(),
                )
                report = format_report(result, cds_start_genomic=cds_start)
                self._last_result = result
                self._last_cds_start = cds_start

                self.output_text.delete("1.0", "end")
                self.output_text.insert("1.0", report)
            except Exception as e:
                messagebox.showerror("Error", str(e))

        def clear_all(self) -> None:
            self.seq_text.delete("1.0", "end")
            self.cds_start_var.set("1")
            self.cds_end_var.set("")
            self.search_start_var.set("")
            self.search_end_var.set("")
            self.system_var.set("Auto")
            self.syn_var.set(True)
            self.output_text.delete("1.0", "end")
            self._last_result = None

        def load_example(self) -> None:
            """Load SCN1A CDS start (first 200 nt of a truncated example for demo)."""
            # Truncated SCN1A-like example (not real sequence, for demonstration)
            example = (
                "ATGGCAGCGCAGCCAGAGGCAGCGGCAGCTCCTGATGTGGCACCGGATCCTCTGAGAACACAAGG"
                "GATGCTCAAGGACAGCAGCTGCACCTTGAGCAGGGGCAGCCTTGGGCAGATGTCGCCTCCTAATC"
                "TGGTCACCCTGCCCAGGAACATCGTCAACCTGGCCCAGGTCATCAGCGCCATCAAGCCCCTGGAG"
                "GATGGAAGCAAGCTGGTGGTCTCTAAGGCCTACAGCTTCTTCCTCAAGCTCAACATTGAGACCATG"
                "GGCGTGCTCCAGGATGTCCTCACCATCATCCTGTCCATGATGATCGTCACCATGATCCTGGAGAAG"
                "AAGATGTTCACCACCATCTTCCTCTTCATCATCGTCTTCGGCATCTTTGGCAACATGCTGGTCATG"
                "CTGGACAAGAATGACCAGGCCAGCCGCCTGGTGCCCACCCTCATCACCATCTCCACCATGGGCTGG"
                "AACCTGGTCTTCACCGGGATCTTCACCGACTTCCCCAAGGTCATCGTCCTGAAGAGCCTGGAGAGC"
            )
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", example)
            self.cds_start_var.set("1")
            self.cds_end_var.set(str(len(clean_sequence(example))))
            self.search_start_var.set("")
            self.search_end_var.set("")
            self.system_var.set("Auto")
            self.output_text.delete("1.0", "end")
            self.output_text.insert("1.0", "Example loaded. Click 'Analyze' to run.")

        def save_txt(self) -> None:
            content = self.output_text.get("1.0", "end").strip()
            if not content:
                messagebox.showinfo("Info", "No results to save.")
                return
            path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text", "*.txt"), ("All", "*.*")])
            if path:
                with open(path, "w", encoding="utf-8") as f:
                    f.write(content)
                messagebox.showinfo("Saved", f"Results saved to {path}")

        def save_csv(self) -> None:
            if self._last_result is None:
                messagebox.showinfo("Info", "Run analysis first.")
                return
            path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=[("CSV", "*.csv"), ("All", "*.*")])
            if path:
                csv_text = format_csv(self._last_result, self._last_cds_start)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(csv_text)
                messagebox.showinfo("Saved", f"CSV saved to {path}")

    root = tk.Tk()
    try:
        root.call("tk", "scaling", 1.05)
    except Exception:
        pass
    App(root)
    root.mainloop()


# ─────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────

def run_cli(args):
    """Run analysis in command-line mode."""
    # Read sequence
    if args.seq_file:
        with open(args.seq_file, "r") as f:
            raw = f.read()
        # Simple FASTA handling
        lines = raw.strip().split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        seq = clean_sequence("".join(seq_lines))
    elif args.seq:
        seq = clean_sequence(args.seq)
    else:
        print("Error: provide --seq or --seq_file", file=sys.stderr)
        sys.exit(1)

    validate_dna(seq)

    cds_start = args.cds_start
    cds_end = args.cds_end or len(seq)
    cds = seq[cds_start - 1:cds_end]

    sys_type = None if args.system == "Auto" else args.system

    result = analyze_cds(
        cds=cds,
        system_type=sys_type,
        search_start=args.search_start or 1,
        search_end=args.search_end or len(cds),
        include_synonymous=not args.no_synonymous,
        max_ratio=args.max_ratio,
        top_n=args.top_n,
    )

    report = format_report(result, cds_start_genomic=cds_start)
    print(report)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"\nReport saved to {args.output}")

    if args.csv_output:
        csv_text = format_csv(result, cds_start_genomic=cds_start)
        with open(args.csv_output, "w", encoding="utf-8") as f:
            f.write(csv_text)
        print(f"CSV saved to {args.csv_output}")


# ─────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="AAVLINK CDS Split-Site Finder v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              GUI mode:   python AAVLINK_SplitFinder.py
              CLI mode:   python AAVLINK_SplitFinder.py --cli --seq_file SCN1A.fa --cds_start 1
              CLI inline: python AAVLINK_SplitFinder.py --cli --seq ATGCCC... --cds_end 6030
        """)
    )

    parser.add_argument("--cli", action="store_true", help="Run in command-line mode (no GUI)")
    parser.add_argument("--seq", type=str, help="CDS/genomic sequence string")
    parser.add_argument("--seq_file", type=str, help="FASTA file with sequence")
    parser.add_argument("--cds_start", type=int, default=1, help="CDS start in sequence (1-based)")
    parser.add_argument("--cds_end", type=int, default=None, help="CDS end in sequence")
    parser.add_argument("--search_start", type=int, default=None, help="Search range start (CDS coords)")
    parser.add_argument("--search_end", type=int, default=None, help="Search range end (CDS coords)")
    parser.add_argument("--system", type=str, default="Auto",
                       choices=["Auto", "Dual", "Triple I", "Triple II"])
    parser.add_argument("--no_synonymous", action="store_true",
                       help="Disable synonymous mutation search")
    parser.add_argument("--max_ratio", type=float, default=3.0,
                       help="Max fragment length ratio (default: 3.0)")
    parser.add_argument("--top_n", type=int, default=20, help="Number of top results")
    parser.add_argument("--output", type=str, help="Save report to file")
    parser.add_argument("--csv_output", type=str, help="Save CSV to file")

    args = parser.parse_args()

    if args.cli:
        run_cli(args)
    else:
        run_gui()


if __name__ == "__main__":
    main()
