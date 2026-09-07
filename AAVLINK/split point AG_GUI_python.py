#!/usr/bin/env python3
"""
AAVLINK CDS Split Finder GUI

기능:
1. 전체 시퀀스를 입력받고 그 안에서 CDS 시작/끝을 지정
2. 전체 시퀀스 기준으로 split 검색 범위를 지정
3. fragment 수는 2개 또는 3개
4. split point는 5' fragment가 AG로 끝나고, 3' fragment가 G로 시작하는 지점만 선택
5. split 후 모든 fragment 길이에 대해 max/min <= 3
6. split point는 coding DNA numbering (예: c.2107)로 표기
   - c.n = n번째 염기부터 3' fragment가 시작
7. split point 앞뒤 nt 20개와 그에 대응하는 a.a. 표시
8. 전체 CDS 선 위에 조건을 만족하는 usable split point만 모두 표시

GUI:
- tkinter 기반
- 입력값 유지
- 결과를 텍스트창에 출력
- 예시/초기화/저장 기능 포함
- 상단(입력/가이드) / 하단(map/결과) 2분할
"""

from __future__ import annotations

import itertools
import re
import textwrap
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Dict, List, Set, Tuple


CODON_TABLE = {
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


def clean_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper()


def validate_dna(seq: str) -> None:
    if not seq:
        raise ValueError("시퀀스가 비어 있습니다.")
    invalid = sorted(set(re.sub(r"[ACGTN]", "", seq.upper())))
    if invalid:
        raise ValueError(
            f"유효하지 않은 문자가 포함되어 있습니다: {', '.join(invalid)}\n"
            "허용 문자는 A, C, G, T, N 입니다."
        )


def translate_dna(seq: str) -> str:
    aa = []
    usable_len = len(seq) - (len(seq) % 3)
    for i in range(0, usable_len, 3):
        codon = seq[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def extract_subseq(seq: str, start_1b: int, end_1b: int) -> str:
    if start_1b < 1 or end_1b > len(seq):
        raise ValueError("지정 범위가 시퀀스 길이를 벗어났습니다.")
    if start_1b > end_1b:
        raise ValueError("시작 위치는 끝 위치보다 작거나 같아야 합니다.")
    return seq[start_1b - 1:end_1b]


def genome_to_cds_coord(genomic_pos_1b: int, cds_start_1b: int, cds_end_1b: int) -> int:
    if genomic_pos_1b < cds_start_1b or genomic_pos_1b > cds_end_1b:
        raise ValueError("입력한 위치가 CDS 범위를 벗어났습니다.")
    return genomic_pos_1b - cds_start_1b + 1


def cds_to_genome_coord(c_pos: int, cds_start_1b: int) -> int:
    return cds_start_1b + c_pos - 1


def find_candidate_split_starts(cds: str, search_start_c: int, search_end_c: int) -> List[int]:
    """
    c.n = CDS의 n번째 염기부터 3' fragment가 시작

    즉 split은 n번째 염기 바로 앞 경계이며,
    조건은 다음과 같다.
    - 5' fragment 끝 2nt가 AG
    - 3' fragment 첫 nt가 G (즉 c.n이 G)

    예:
        ...AG|G...
             ^
            c.n
    """
    candidates = []
    cds_len = len(cds)

    if search_start_c < 3:
        search_start_c = 3
    if search_end_c > cds_len:
        search_end_c = cds_len

    for n in range(search_start_c, search_end_c + 1):
        right_first = cds[n - 1]
        left_last2 = cds[n - 3:n - 1]
        if left_last2 == "AG" and right_first == "G":
            candidates.append(n)

    return candidates


def fragment_lengths(cds_len: int, split_starts: Tuple[int, ...]) -> List[int]:
    """
    split_start c.n 은 n번째부터 다음 fragment 시작을 뜻함.
    예:
      split_start = c.7
      -> 첫 fragment 길이 = 6 nt (c.1~c.6)
      -> 둘째 fragment는 c.7부터 시작
    """
    lengths = []
    prev_start = 1

    for sp in split_starts:
        lengths.append(sp - prev_start)
        prev_start = sp

    lengths.append(cds_len - prev_start + 1)
    return lengths


def ratio_ok(lengths: List[int], max_ratio: float = 3.0) -> bool:
    min_len = min(lengths)
    max_len = max(lengths)
    if min_len == 0:
        return False
    return (max_len / min_len) <= max_ratio


def get_nt_context(cds: str, split_start: int, flank: int = 20) -> str:
    split0 = split_start - 1
    left_start0 = max(0, split0 - flank)
    right_end0 = min(len(cds), split0 + flank)

    left_seq = cds[left_start0:split0]
    right_seq = cds[split0:right_end0]
    return f"{left_seq}|{right_seq}"


def get_aa_context(cds: str, split_start: int, flank_nt: int = 20) -> str:
    split0 = split_start - 1

    start_nt_1b = max(1, split_start - flank_nt)
    end_nt_1b = min(len(cds), split_start + flank_nt - 1)

    frame0_start = start_nt_1b - 1
    frame0_end = end_nt_1b - 1

    codon_start0 = (frame0_start // 3) * 3
    codon_end0 = ((frame0_end // 3) + 1) * 3
    codon_end0 = min(codon_end0, len(cds))

    sub_nt = cds[codon_start0:codon_end0]
    sub_aa = translate_dna(sub_nt)

    local_left_aa_count = max(0, split0 // 3 - codon_start0 // 3)
    split_inside = split0 % 3

    aa_prefix = sub_aa[:local_left_aa_count]
    aa_suffix = sub_aa[local_left_aa_count:]

    if split_inside == 0:
        marker = "|"
    elif split_inside == 1:
        marker = "[split, codon+1]|"
    else:
        marker = "[split, codon+2]|"

    return f"{aa_prefix}{marker}{aa_suffix}"


def _distribute_labels(width: int, label_positions: List[Tuple[int, str]]) -> List[str]:
    rows: List[List[str]] = []
    occupied_ranges: List[List[Tuple[int, int]]] = []

    for pos, label in label_positions:
        label_len = len(label)
        start = max(0, min(width - label_len, pos - label_len // 2))
        end = start + label_len - 1

        placed = False
        for row_idx, ranges in enumerate(occupied_ranges):
            overlap = any(not (end < s or start > e) for s, e in ranges)
            if not overlap:
                rows[row_idx][start:start + label_len] = list(label)
                ranges.append((start, end))
                placed = True
                break

        if not placed:
            new_row = [" "] * width
            new_row[start:start + label_len] = list(label)
            rows.append(new_row)
            occupied_ranges.append([(start, end)])

    return ["".join(r) for r in rows]


def format_valid_points_block(valid_points: List[int], width: int = 90) -> str:
    if not valid_points:
        return "Valid split point(s): 없음"

    label = "Valid split point(s): "
    joined = ", ".join(f"c.{x}" for x in valid_points)

    wrapped = textwrap.wrap(
        joined,
        width=max(20, width - len(label)),
        break_long_words=False,
        break_on_hyphens=False,
    )

    if not wrapped:
        return label

    lines = [label + wrapped[0]]
    indent = " " * len(label)
    for line in wrapped[1:]:
        lines.append(indent + line)
    return "\n".join(lines)


def make_ascii_map_all_valid_points(cds_len: int, valid_points: List[int], width: int = 100) -> str:
    width = max(width, 60)
    base_line = ["-"] * width
    arrow_line = [" "] * width
    plotted_positions: List[Tuple[int, str]] = []

    for sp in valid_points:
        pos = min(width - 1, max(0, round(((sp - 1) / max(1, cds_len - 1)) * (width - 1))))
        base_line[pos] = "|"
        arrow_line[pos] = "↑"
        plotted_positions.append((pos, f"c.{sp}"))

    label_lines = _distribute_labels(width, plotted_positions)

    map_parts = [
        f"CDS length: {cds_len} bp",
        f"5' {''.join(base_line)} 3'",
        f"   {''.join(arrow_line)}",
    ]
    for line in label_lines:
        map_parts.append(f"   {line}")

    return "\n".join(map_parts)


def enumerate_valid_splits(cds: str, fragment_count: int, candidate_points: List[int]) -> List[Dict]:
    if fragment_count not in (2, 3):
        raise ValueError("fragment 수는 2 또는 3만 가능합니다.")

    cds_len = len(cds)
    results = []
    num_splits = fragment_count - 1

    for combo in itertools.combinations(candidate_points, num_splits):
        lengths = fragment_lengths(cds_len, combo)
        if ratio_ok(lengths, max_ratio=3.0):
            split_infos = []
            for sp in combo:
                split_infos.append({
                    "split_point": sp,
                    "nt_context": get_nt_context(cds, sp, flank=20),
                    "aa_context": get_aa_context(cds, sp, flank_nt=20),
                })

            results.append({
                "split_points": combo,
                "fragment_lengths": lengths,
                "split_infos": split_infos,
            })

    return results


def collect_points_used_in_valid_results(results: List[Dict]) -> List[int]:
    used: Set[int] = set()
    for result in results:
        for sp in result["split_points"]:
            used.add(sp)
    return sorted(used)


def format_results(results: List[Dict], cds_start: int) -> str:
    if not results:
        return "조건을 만족하는 split 조합이 없습니다."

    out = []
    for idx, result in enumerate(results, start=1):
        out.append("=" * 100)
        out.append(f"[결과 {idx}]")
        out.append(f"Split point: {', '.join(f'c.{x}' for x in result['split_points'])}")
        out.append(f"Fragment length(nt): {', '.join(map(str, result['fragment_lengths']))}")
        out.append("")

        for j, info in enumerate(result["split_infos"], start=1):
            genome_pos = cds_to_genome_coord(info["split_point"], cds_start)
            out.append(f"  [split {j}] c.{info['split_point']} (전체 시퀀스 좌표: {genome_pos})")
            out.append(f"  - nt  : {info['nt_context']}")
            out.append(f"  - a.a.: {info['aa_context']}")
            out.append("")

    return "\n".join(out)


def run_analysis(
    full_seq: str,
    cds_start: int,
    cds_end: int,
    search_start_genomic: int,
    search_end_genomic: int,
    fragment_count: int,
) -> Tuple[str, str]:
    cds = extract_subseq(full_seq, cds_start, cds_end)

    search_start_c = genome_to_cds_coord(search_start_genomic, cds_start, cds_end)
    search_end_c = genome_to_cds_coord(search_end_genomic, cds_start, cds_end)

    candidate_points = find_candidate_split_starts(
        cds=cds,
        search_start_c=search_start_c,
        search_end_c=search_end_c,
    )

    results = enumerate_valid_splits(
        cds=cds,
        fragment_count=fragment_count,
        candidate_points=candidate_points,
    )

    valid_points_for_map = collect_points_used_in_valid_results(results)

    fixed_map_text = "[전체 valid split point map]\n\n"
    if valid_points_for_map:
        fixed_map_text += make_ascii_map_all_valid_points(len(cds), valid_points_for_map)
        fixed_map_text += "\n\n"
        fixed_map_text += format_valid_points_block(valid_points_for_map, width=90)
    else:
        fixed_map_text += "max/min <= 3 조건까지 만족하는 split point가 없습니다."

    out = []
    out.append("=" * 100)
    out.append("[현재 입력값 요약]")
    out.append(f"전체 시퀀스 길이: {len(full_seq)} bp")
    out.append(f"CDS 범위: {cds_start}..{cds_end} (길이 {len(cds)} bp)")
    out.append(f"split 검색 범위: {search_start_genomic}..{search_end_genomic} (전체 시퀀스 기준)")
    out.append(f"fragment 수: {fragment_count}")
    out.append("")
    out.append(f"candidate split point 수 (AG|G): {len(candidate_points)}")
    out.append(f"valid split point 수 (max/min <= 3 포함): {len(valid_points_for_map)}")
    out.append("")
    out.append("[검색 결과]")
    out.append(format_results(results, cds_start))

    return fixed_map_text, "\n".join(out)


class App:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("AAVLINK CDS Split Finder GUI")
        self.root.geometry("1000x760")
        self.root.minsize(850, 650)
        self._build_ui()

    def _build_ui(self) -> None:
        main = ttk.Frame(self.root, padding=8)
        main.pack(fill="both", expand=True)

        outer_paned = ttk.Panedwindow(main, orient="vertical")
        outer_paned.pack(fill="both", expand=True)

        top_container = ttk.Frame(outer_paned)
        bottom_container = ttk.Frame(outer_paned)

        outer_paned.add(top_container, weight=3)
        outer_paned.add(bottom_container, weight=4)

        self._build_top_area(top_container)
        self._build_bottom_area(bottom_container)

    def _build_top_area(self, parent: ttk.Frame) -> None:
        top_paned = ttk.Panedwindow(parent, orient="horizontal")
        top_paned.pack(fill="both", expand=True)

        input_frame = ttk.LabelFrame(top_paned, text="입력", padding=10)
        guide_frame = ttk.LabelFrame(top_paned, text="입력 가이드", padding=10)

        top_paned.add(input_frame, weight=5)
        top_paned.add(guide_frame, weight=3)

        ttk.Label(input_frame, text="전체 DNA sequence").grid(row=0, column=0, sticky="w")
        self.seq_text = scrolledtext.ScrolledText(input_frame, width=80, height=10, wrap="word")
        self.seq_text.grid(row=1, column=0, columnspan=6, sticky="nsew", pady=(5, 10))

        ttk.Label(input_frame, text="CDS 시작").grid(row=2, column=0, sticky="w")
        self.cds_start_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.cds_start_var, width=10).grid(row=2, column=1, sticky="w", padx=(5, 15))

        ttk.Label(input_frame, text="CDS 끝").grid(row=2, column=2, sticky="w")
        self.cds_end_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.cds_end_var, width=10).grid(row=2, column=3, sticky="w", padx=(5, 15))

        ttk.Label(input_frame, text="split 범위 시작").grid(row=3, column=0, sticky="w", pady=(10, 0))
        self.search_start_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.search_start_var, width=10).grid(row=3, column=1, sticky="w", padx=(5, 15), pady=(10, 0))

        ttk.Label(input_frame, text="split 범위 끝").grid(row=3, column=2, sticky="w", pady=(10, 0))
        self.search_end_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.search_end_var, width=10).grid(row=3, column=3, sticky="w", padx=(5, 15), pady=(10, 0))

        ttk.Label(input_frame, text="fragment 수").grid(row=3, column=4, sticky="w", pady=(10, 0))
        self.fragment_var = tk.StringVar(value="2")
        ttk.Combobox(
            input_frame,
            textvariable=self.fragment_var,
            values=["2", "3"],
            width=8,
            state="readonly"
        ).grid(row=3, column=5, sticky="w", padx=(5, 0), pady=(10, 0))

        button_frame = ttk.Frame(input_frame)
        button_frame.grid(row=4, column=0, columnspan=6, sticky="w", pady=(15, 0))

        ttk.Button(button_frame, text="분석 실행", command=self.analyze).pack(side="left", padx=(0, 8))
        ttk.Button(button_frame, text="예시 입력", command=self.load_example).pack(side="left", padx=(0, 8))
        ttk.Button(button_frame, text="초기화", command=self.clear_all).pack(side="left", padx=(0, 8))
        ttk.Button(button_frame, text="결과 저장", command=self.save_result).pack(side="left")

        input_frame.columnconfigure(0, weight=0)
        input_frame.columnconfigure(1, weight=0)
        input_frame.columnconfigure(2, weight=0)
        input_frame.columnconfigure(3, weight=0)
        input_frame.columnconfigure(4, weight=0)
        input_frame.columnconfigure(5, weight=1)
        input_frame.rowconfigure(1, weight=1)

        help_text = (
            "- split 검색 범위는 전체 시퀀스 기준으로 입력하지만, 결과 표기는 CDS 기준 c.n 으로 표시됩니다.\n\n"
            "- c.n 은 CDS의 n번째 염기부터 3' fragment가 시작된다는 뜻입니다.\n\n"
            "- split motif는 AG|G\n"
            "  (5' CDS가 AG로 끝나고 3' CDS가 G로 시작) 입니다.\n\n"
            "- [split, codon+1], [split, codon+2]는 codon 중간 split을 뜻합니다.\n"
            "  예: codon+2 = AG|G\n\n"
            "- fragment 길이 조건은 max/min <= 3 입니다."
        )

        self.guide_text = scrolledtext.ScrolledText(guide_frame, width=35, height=10, wrap="word")
        self.guide_text.pack(fill="both", expand=True)
        self.guide_text.insert("1.0", help_text)
        self.guide_text.configure(state="disabled", font=("Arial", 10))

    def _build_bottom_area(self, parent: ttk.Frame) -> None:
        bottom_paned = ttk.Panedwindow(parent, orient="vertical")
        bottom_paned.pack(fill="both", expand=True)

        map_frame = ttk.LabelFrame(bottom_paned, text="Split point map", padding=10)
        output_frame = ttk.LabelFrame(bottom_paned, text="결과", padding=10)

        bottom_paned.add(map_frame, weight=2)
        bottom_paned.add(output_frame, weight=5)

        self.map_text = scrolledtext.ScrolledText(map_frame, width=120, height=10, wrap="word")
        self.map_text.pack(fill="both", expand=True)
        self.map_text.configure(font=("Consolas", 10))

        self.output_text = scrolledtext.ScrolledText(output_frame, width=120, height=18, wrap="none")
        self.output_text.pack(fill="both", expand=True)
        self.output_text.configure(font=("Consolas", 10))

    def _read_inputs(self) -> Tuple[str, int, int, int, int, int]:
        full_seq = clean_sequence(self.seq_text.get("1.0", "end"))
        validate_dna(full_seq)

        cds_start = int(self.cds_start_var.get())
        cds_end = int(self.cds_end_var.get())
        search_start = int(self.search_start_var.get())
        search_end = int(self.search_end_var.get())
        fragment_count = int(self.fragment_var.get())

        _ = extract_subseq(full_seq, cds_start, cds_end)

        if search_start > search_end:
            raise ValueError("split 범위 시작은 끝보다 작거나 같아야 합니다.")
        if search_start < cds_start or search_end > cds_end:
            raise ValueError("split 검색 범위는 CDS 내부여야 합니다.")
        if fragment_count not in (2, 3):
            raise ValueError("fragment 수는 2 또는 3만 가능합니다.")

        return full_seq, cds_start, cds_end, search_start, search_end, fragment_count

    def analyze(self) -> None:
        try:
            full_seq, cds_start, cds_end, search_start, search_end, fragment_count = self._read_inputs()
            map_result, result = run_analysis(
                full_seq=full_seq,
                cds_start=cds_start,
                cds_end=cds_end,
                search_start_genomic=search_start,
                search_end_genomic=search_end,
                fragment_count=fragment_count,
            )

            self.map_text.delete("1.0", "end")
            self.map_text.insert("1.0", map_result)

            self.output_text.delete("1.0", "end")
            self.output_text.insert("1.0", result)

        except Exception as e:
            messagebox.showerror("오류", str(e))

    def clear_all(self) -> None:
        self.seq_text.delete("1.0", "end")
        self.cds_start_var.set("")
        self.cds_end_var.set("")
        self.search_start_var.set("")
        self.search_end_var.set("")
        self.fragment_var.set("2")
        self.map_text.delete("1.0", "end")
        self.output_text.delete("1.0", "end")

    def load_example(self) -> None:
        example_seq = (
            "ATGGCCGAGGCTCAGGGTGAGGCTGAGGAGGAGGCTGAGGGTGAGGCTAAGGGTGAG"
            "GCTCAGGGTGAGGCTGAGGAGGAGGCTGAGGGTGAGGCTCAGGGTGAGGCTGAGGAG"
        )
        self.seq_text.delete("1.0", "end")
        self.seq_text.insert("1.0", example_seq)
        self.cds_start_var.set("1")
        self.cds_end_var.set(str(len(example_seq)))
        self.search_start_var.set("1")
        self.search_end_var.set(str(len(example_seq)))
        self.fragment_var.set("2")
        self.map_text.delete("1.0", "end")
        self.output_text.delete("1.0", "end")
        self.output_text.insert("1.0", "예시 입력이 채워졌습니다. '분석 실행'을 눌러 주세요.")

    def save_result(self) -> None:
        map_content = self.map_text.get("1.0", "end").strip()
        result_content = self.output_text.get("1.0", "end").strip()
        content = "\n\n".join(part for part in [map_content, result_content] if part)

        if not content:
            messagebox.showinfo("안내", "저장할 결과가 없습니다.")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="결과 저장"
        )
        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(content)
            messagebox.showinfo("저장 완료", f"결과를 저장했습니다.\n{file_path}")
        except Exception as e:
            messagebox.showerror("저장 실패", str(e))


def main() -> None:
    root = tk.Tk()
    try:
        root.call("tk", "scaling", 1.05)
    except Exception:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()