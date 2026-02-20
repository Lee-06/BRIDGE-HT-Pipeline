#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse
from pathlib import Path

DELIMS = set([":", ",", ")", ";"])


def load_mapping_tsv(path: Path) -> list[tuple[str, str]]:
    """
    Read 2-column TSV mapping: old<TAB>new
    Keeps order (sequential application).
    """
    pairs: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"Mapping TSV invalid at line {line_no}: expected 2 columns (tab-separated)."
                )
            old = parts[0].strip()
            new = parts[1].strip()
            if not old:
                raise ValueError(f"Mapping TSV invalid at line {line_no}: empty 'old' value.")
            pairs.append((old, new))
    return pairs


def should_transform_label(label: str) -> bool:
    """
    Avoid touching numeric internal supports like '34' or '100'.
    Transform leaf labels (usually contain letters, underscores, pipes, dots, etc.).
    """
    s = label.strip()
    if not s:
        return False
    if s.replace(".", "", 1).isdigit():
        return False
    return True


def transform_text(s: str, mapping: list[tuple[str, str]], remove_assembly_scaffolds: bool) -> str:
    """
    Apply sequential replacements.
    """
    if remove_assembly_scaffolds:
        s = s.replace("_AssemblyScaffolds.fasta", "")
    for old, new in mapping:
        s = s.replace(old, new)
    return s


# ---------- NEWICK ----------
def rewrite_newick_labels(newick: str, mapping: list[tuple[str, str]], remove_assembly_scaffolds: bool) -> str:
    """
    Minimal Newick label rewriter:
    - rewrites labels that appear after '(' or ',' and before ':'/','/')'/';'
    - DOES NOT rewrite internal node labels right after ')'
    """
    out = []
    i = 0
    n = len(newick)

    while i < n:
        ch = newick[i]

        if ch in "(,":
            out.append(ch)
            i += 1

            start = i
            while i < n and newick[i] not in DELIMS:
                i += 1
            label = newick[start:i]

            if label and should_transform_label(label):
                label = transform_text(label, mapping, remove_assembly_scaffolds)

            out.append(label)
            continue

        out.append(ch)
        i += 1

    return "".join(out)


def process_tree_file(in_path: Path, out_path: Path, mapping: list[tuple[str, str]], remove_assembly_scaffolds: bool) -> None:
    txt = in_path.read_text(encoding="utf-8")
    new_txt = rewrite_newick_labels(txt, mapping, remove_assembly_scaffolds=remove_assembly_scaffolds)
    out_path.write_text(new_txt, encoding="utf-8")


# ---------- FASTA ----------
def rewrite_fasta_headers(fasta_text: str, mapping: list[tuple[str, str]], remove_assembly_scaffolds: bool) -> str:
    """
    Rewrites only FASTA headers:
    - line startswith '>'
    - transforms only the identifier token (up to first whitespace)
    - keeps the rest of the header line unchanged
    - does NOT touch sequence lines
    """
    out_lines: list[str] = []
    for line in fasta_text.splitlines(keepends=False):
        if line.startswith(">"):
            header = line[1:]  # without '>'
            if header.strip() == "":
                out_lines.append(line)
                continue

            # split on first whitespace (identifier vs description)
            parts = header.split(None, 1)
            ident = parts[0]
            rest = parts[1] if len(parts) == 2 else ""

            new_ident = transform_text(ident, mapping, remove_assembly_scaffolds)

            if rest:
                out_lines.append(f">{new_ident} {rest}")
            else:
                out_lines.append(f">{new_ident}")
        else:
            out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if fasta_text.endswith("\n") else "")


def process_fasta_file(in_path: Path, out_path: Path, mapping: list[tuple[str, str]], remove_assembly_scaffolds: bool) -> None:
    txt = in_path.read_text(encoding="utf-8")
    new_txt = rewrite_fasta_headers(txt, mapping, remove_assembly_scaffolds=remove_assembly_scaffolds)
    out_path.write_text(new_txt, encoding="utf-8")


def guess_mode_from_pattern(pattern: str) -> str:
    p = pattern.lower()
    if p.endswith(".treefile") or p.endswith(".nwk") or p.endswith(".newick") or "tree" in p:
        return "tree"
    if p.endswith(".fa") or p.endswith(".fna") or p.endswith(".faa") or p.endswith(".fasta") or "fasta" in p:
        return "fasta"
    return "fasta"  # default sensible


def main():
    p = argparse.ArgumentParser(
        description="Rename IDs in Newick tree files OR FASTA headers using a TSV mapping; writes results to a new folder."
    )
    p.add_argument("-i", "--input-dir", required=True, help="Folder containing input files")
    p.add_argument("-m", "--mapping-tsv", required=True, help="TSV mapping: old<TAB>new")
    p.add_argument("-o", "--output-dir", required=True, help="Output folder for rewritten files")
    p.add_argument("--pattern", required=True, help="Glob pattern to select files (e.g. *.treefile or *.fasta)")
    p.add_argument("--mode", choices=["auto", "tree", "fasta"], default="auto",
                   help="How to process files: tree (Newick), fasta (FASTA headers), or auto (guess from pattern).")
    p.add_argument("--no-remove-assembly-scaffolds", action="store_true",
                   help="Disable removal of '_AssemblyScaffolds.fasta' in labels/headers")
    args = p.parse_args()

    in_dir = Path(args.input_dir).expanduser().resolve()
    out_dir = Path(args.output_dir).expanduser().resolve()
    mapping_path = Path(args.mapping_tsv).expanduser().resolve()

    if not in_dir.is_dir():
        raise SystemExit(f"Input dir not found: {in_dir}")
    if not mapping_path.is_file():
        raise SystemExit(f"Mapping TSV not found: {mapping_path}")

    mapping = load_mapping_tsv(mapping_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(in_dir.glob(args.pattern))
    if not files:
        raise SystemExit(f"No files matched pattern '{args.pattern}' in {in_dir}")

    remove_assembly_scaffolds = not args.no_remove_assembly_scaffolds

    mode = args.mode
    if mode == "auto":
        mode = guess_mode_from_pattern(args.pattern)

    for f in files:
        out_path = out_dir / f.name
        if mode == "tree":
            process_tree_file(f, out_path, mapping, remove_assembly_scaffolds)
        else:
            process_fasta_file(f, out_path, mapping, remove_assembly_scaffolds)

    print(f"Done. Rewritten {len(files)} file(s) into: {out_dir} (mode={mode})")


if __name__ == "__main__":
    main()
