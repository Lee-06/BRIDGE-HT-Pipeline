#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO


DEFAULT_KEYWORDS = [
    "ribosomal", "ATPase", "18s", "28s", "5s", "rrna", "rdna",
    "ribonucleoprotein", "translation elongation factor", "elongation factor",
    "mitochondrial", "mitochondrion", "cytochrome",
    "cox1", "cox2", "cox3",
    "nad", "nadh dehydrogenase",
    "atp6", "atp8", "atp9", "atp synthase",
    "chloroplast", "plastid", "rubisco", "rbcl",
    "dna polymerase", "rna polymerase",
    "helicase", "topoisomerase", "histone",
    "ligase", "exonuclease", "primase",
    "actin", "tubulin", "kinesin", "dynein", "myosin",
    "heat shock protein", "hsp", "chaperone",
    "ubiquitin", "kinase",
]


def eprint(*args):
    print(*args, file=sys.stderr)


def unique_path(path: Path) -> Path:
    """If path exists, create path_v2, path_v3, ... (no overwrite)."""
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    i = 2
    while True:
        candidate = parent / f"{stem}_v{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def write_default_keyword_file(outdir: Path) -> Path:
    """
    Write default keyword file into outdir without overwriting.
    Comma-separated keywords.
    """
    kw_path = unique_path(outdir / "default_housekeeping_keywords.csv")
    kw_path.write_text(",".join(DEFAULT_KEYWORDS) + "\n", encoding="utf-8")
    print(f"[INFO] Default keyword file written: {kw_path}")
    return kw_path


def load_keywords_csv(path: Path) -> List[str]:
    """
    Load comma-separated keywords from file (multi-line allowed).
    Keep original case; we do case-insensitive matching later.
    """
    txt = path.read_text(encoding="utf-8", errors="replace")
    kws: List[str] = []
    for part in txt.replace("\n", ",").split(","):
        kw = part.strip()
        if kw:
            kws.append(kw)

    # unique, keep order (case-insensitive uniqueness)
    seen = set()
    out: List[str] = []
    for k in kws:
        key = k.casefold()
        if key not in seen:
            out.append(k)
            seen.add(key)
    return out


def normalize_query_id(q: str) -> str:
    """
    EggNOG query IDs often end with _0, _1, etc. Remove ONLY trailing _<digits>.
    """
    if "_" not in q:
        return q
    left, right = q.rsplit("_", 1)
    if right.isdigit():
        return left
    return q


def parse_eggnog_annotations(annotations: Path) -> Dict[str, Dict[str, str]]:
    """
    Return dict keyed by normalized query id with fields:
      description, preferred_name, pfams, cog_category, seed_ortholog
    """
    header_line = None
    with annotations.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#query\t"):
                header_line = line.strip()[1:]  # remove leading '#'
                break
    if header_line is None:
        raise SystemExit("[ERROR] Could not find '#query' header line in annotations file.")

    cols = header_line.split("\t")

    df = pd.read_csv(
        annotations,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=str
    )

    info: Dict[str, Dict[str, str]] = {}
    for _, row in df.iterrows():
        q = str(row.get("query", "")).strip()
        if not q:
            continue
        qn = normalize_query_id(q)

        # Keep first hit per query (good enough for filtering/reporting)
        if qn not in info:
            info[qn] = {
                "query_raw": q,
                "description": str(row.get("Description", "")).strip(),
                "preferred_name": str(row.get("Preferred_name", "")).strip(),
                "pfams": str(row.get("PFAMs", "")).strip(),
                "cog_category": str(row.get("COG_category", "")).strip(),
                "seed_ortholog": str(row.get("seed_ortholog", "")).strip(),
            }

    return info


def _compile_whole_word_patterns(keywords: List[str]) -> List[Tuple[str, re.Pattern]]:
    """
    Compile regex patterns that match keywords as WHOLE WORDS (or whole phrases),
    case-insensitive.

    Rule:
    - For each keyword, we match it only when it's not embedded inside a larger
      alphanumeric/underscore token.
    - This avoids false positives like 'nad' matching inside 'canada' or 'nadh'.
    - Works for phrases too: 'atp synthase' must appear as a phrase with boundaries
      around the full phrase.

    Boundary definition:
    - left boundary: start of string OR a character that is NOT [A-Za-z0-9_]
    - right boundary: end of string OR a character that is NOT [A-Za-z0-9_]

    Using a custom boundary rather than \b because \b treats '_' as a word char,
    and we want consistent behavior for IDs / PFAM-like tokens.
    """
    compiled: List[Tuple[str, re.Pattern]] = []
    for kw in keywords:
        kw_strip = kw.strip()
        if not kw_strip:
            continue

        # Escape regex meta-chars in the keyword
        escaped = re.escape(kw_strip)

        # Custom "whole token/phrase" boundaries
        pattern = re.compile(rf"(^|[^A-Za-z0-9_]){escaped}($|[^A-Za-z0-9_])", flags=re.IGNORECASE)
        compiled.append((kw_strip, pattern))
    return compiled


def find_matches_for_filtering(
    eggnog_info: Dict[str, Dict[str, str]],
    keywords: List[str]
) -> Dict[str, Tuple[str, str, str]]:
    """
    Filtering uses Description + Preferred_name + PFAMs.
    Whole-word/whole-phrase, case-insensitive matching to reduce false positives.
    Returns:
      id -> (matched_keyword, description, preferred_name)
    """
    matches: Dict[str, Tuple[str, str, str]] = {}
    patterns = _compile_whole_word_patterns(keywords)

    for qid, meta in eggnog_info.items():
        desc = meta.get("description", "")
        pref = meta.get("preferred_name", "")
        pfams = meta.get("pfams", "")

        haystack = f"{desc} {pref} {pfams}"

        hit_kw = None
        for kw_raw, pat in patterns:
            if pat.search(haystack):
                hit_kw = kw_raw
                break

        if hit_kw:
            matches[qid] = (hit_kw, desc, pref)

    return matches


def main():
    ap = argparse.ArgumentParser(
        description="Step 8: Filter highly conserved genes using EggNOG annotations (keywords). "
                    "Matches are case-insensitive and whole-word/whole-phrase to reduce false positives. "
                    "Search includes Description, Preferred_name and PFAMs."
    )
    ap.add_argument("--annotations", required=True, type=Path,
                    help="EggNOG .emapper.annotations from Step 7")
    ap.add_argument("--fasta-in", required=True, type=Path,
                    help="Clustered FASTA from Step 6 (e.g., Result_HT/ht_clusters.fasta)")

    ap.add_argument("--outdir", default=".", type=Path,
                    help="Output directory (default: current directory)")

    ap.add_argument("--fasta-out", default="ht_clusters.housekeeping_filtered.fasta", type=str,
                    help="Output FASTA filename (in --outdir). If exists, _v2 suffix is used.")
    ap.add_argument("--summary", default="housekeeping_filter_summary.tsv", type=str,
                    help="Summary TSV filename (in --outdir). If exists, _v2 suffix is used.")

    ap.add_argument("--keyword-file", default=None, type=Path,
                    help="User keyword file (comma-separated). If provided, ONLY this file is used.")
    ap.add_argument("--no-filter", action="store_true",
                    help="Disable filtering (script runs but keeps everything; warns about runtime/false positives).")

    args = ap.parse_args()

    if not args.annotations.exists():
        sys.exit(f"[ERROR] annotations not found: {args.annotations}")
    if not args.fasta_in.exists():
        sys.exit(f"[ERROR] fasta-in not found: {args.fasta_in}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Keyword handling
    if args.no_filter:
        eprint(
            "\n[WARNING] Housekeeping filtering is DISABLED (--no-filter).\n"
            "          Potential true HT candidates will be retained,\n"
            "          BUT false positives and downstream runtime may increase significantly.\n"
        )
        keywords: List[str] = []
        keyword_source = "NO_FILTER"
    else:
        if args.keyword_file is not None:
            if not args.keyword_file.exists():
                sys.exit(f"[ERROR] keyword-file not found: {args.keyword_file}")
            keyword_path = args.keyword_file
            print(f"[INFO] Using USER keyword file ONLY: {keyword_path}")
        else:
            keyword_path = write_default_keyword_file(args.outdir)
            print(f"[INFO] Using DEFAULT keyword file: {keyword_path}")

        keywords = load_keywords_csv(keyword_path)
        keyword_source = str(keyword_path)

        if not keywords:
            eprint(
                "\n[WARNING] Keyword list is empty. Nothing will be filtered.\n"
                "          (You can still continue, but expect many conserved genes.)\n"
            )

    # Parse EggNOG annotations
    eggnog_info = parse_eggnog_annotations(args.annotations)

    # Determine removals (whole-word matching on Description+Preferred_name+PFAMs)
    matches = find_matches_for_filtering(eggnog_info, keywords) if keywords else {}
    ids_to_remove = set(matches.keys())

    # Output paths (no overwrite)
    fasta_out_path = unique_path(args.outdir / args.fasta_out)
    summary_path = unique_path(args.outdir / args.summary)

    kept = 0
    removed = 0

    with fasta_out_path.open("w", encoding="utf-8") as fout, summary_path.open("w", encoding="utf-8") as s:
        s.write("id\tstatus\treason\tmatched_keyword\tdescription\tpreferred_name\tpfams\tkeyword_source\n")

        for rec in SeqIO.parse(str(args.fasta_in), "fasta"):
            rid = rec.id
            meta = eggnog_info.get(rid, {})
            desc = meta.get("description", "")
            pref = meta.get("preferred_name", "")
            pfams = meta.get("pfams", "")

            if rid in ids_to_remove:
                kw, _, _ = matches[rid]
                s.write(f"{rid}\tremoved\tEggNOG_keyword_wholeword\t{kw}\t{desc}\t{pref}\t{pfams}\t{keyword_source}\n")
                removed += 1
            else:
                SeqIO.write(rec, fout, "fasta")
                s.write(f"{rid}\tkept\t-\t-\t{desc}\t{pref}\t{pfams}\t{keyword_source}\n")
                kept += 1

    print("[SUCCESS] Step 8 complete.")
    print(f"[INFO] Output FASTA  : {fasta_out_path}")
    print(f"[INFO] Summary TSV   : {summary_path}")
    print(f"[INFO] Kept          : {kept}")
    print(f"[INFO] Removed       : {removed}")


if __name__ == "__main__":
    main()
