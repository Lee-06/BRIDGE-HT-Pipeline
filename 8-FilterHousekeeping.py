#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO


DEFAULT_KEYWORDS = [
    "18s", "28s", "5s",
    "actin",
    "atp synthase", "ATP-synt", "ATP1", "ATP2", "ATPase", "atp6", "atp8", "atp9",
    "CBF5",
    "CDC48",
    "chaperone",
    "chloroplast",
    "cox1", "cox2", "cox3",
    "cytochrome",
    "DHH1",
    "dna polymerase",
    "dynein",
    "elongation factor", "translation elongation factor",
    "exonuclease",
    "FAL1",
    "GAPDH",
    "heat shock protein", "HSP70", "HSP82", "hsp",
    "histone",
    "HRR25",
    "kinesin",
    "ligase",
    "mitochondrial", "mitochondrion",
    "myosin",
    "nad", "nadh dehydrogenase",
    "NDH51",
    "NOG1",
    "plastid",
    "primase",
    "PRP8", "PRPF8", "prp10",
    "rbcl", "rubisco",
    "rdna", "rrna",
    "RNR1", "RNR2",
    "ribonucleoprotein",
    "RPL10", "RPL3",
    "RPS3","RPS2"
    "rna polymerase",
    "RPT1", "RPT6",
    "RVB2",
    "SNZ1",
    "Splicing",
    "SSA2",
    "topoisomerase",
    "tubulin",
    "ubiquitin",
    "VMA2",
]


def unique_path(path: Path) -> Path:
    if not path.exists():
        return path
    i = 2
    while True:
        p2 = path.with_name(f"{path.stem}_v{i}{path.suffix}")
        if not p2.exists():
            return p2
        i += 1


def write_default_keyword_file(outdir: Path) -> Path:
    out = unique_path(outdir / "default_housekeeping_keywords.csv")
    out.write_text(",".join(DEFAULT_KEYWORDS) + "\n", encoding="utf-8")
    print(f"[INFO] Default keyword file written: {out}")
    return out


def load_keywords_csv(path: Path) -> List[str]:
    txt = path.read_text(encoding="utf-8", errors="replace")
    parts = [p.strip() for p in txt.replace("\n", ",").split(",")]
    kws = [p for p in parts if p]
    # unique (case-insensitive), keep order
    seen = set()
    out = []
    for k in kws:
        kk = k.casefold()
        if kk not in seen:
            out.append(k)
            seen.add(kk)
    return out


def normalize_id(s: str) -> str:
    """
    Normalize both EggNOG and FASTA IDs:
    - keep before first '|'
    - remove trailing _<digits> if present
    """
    s = (s or "").strip()
    if "|" in s:
        s = s.split("|", 1)[0]
    if "_" in s:
        left, right = s.rsplit("_", 1)
        if right.isdigit():
            s = left
    return s


def compile_patterns(keywords: List[str]) -> List[Tuple[str, re.Pattern]]:
    pats = []
    for kw in keywords:
        esc = re.escape(kw.strip())
        # same boundary logic as before
        pat = re.compile(rf"(^|[^A-Za-z0-9_]){esc}($|[^A-Za-z0-9_])", re.I)
        pats.append((kw, pat))
    return pats


def parse_emapper_annotations(path: Path) -> Dict[str, Dict[str, str]]:
    """
    Parse EggNOG emapper.annotations without relying on the '#query' header line.
    Uses fixed column positions matching EggNOG header you showed:

    0 query
    1 seed_ortholog
    2 evalue
    3 score
    4 eggNOG_OGs
    5 max_annot_lvl
    6 COG_category
    7 Description
    8 Preferred_name
    ...
    19 PFAMs
    """
    info: Dict[str, Dict[str, str]] = {}

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 20:
                continue

            q_raw = cols[0].strip()
            q = normalize_id(q_raw)

            desc = cols[7].strip()
            pref = cols[8].strip()
            pfam = cols[19].strip() if len(cols) > 19 else ""

            # merge multiple ORFs (_0/_1/...) into same base id
            if q not in info:
                info[q] = {"description": desc, "preferred_name": pref, "pfams": pfam}
            else:
                if desc and desc not in info[q]["description"]:
                    info[q]["description"] += " ; " + desc
                if pref and pref not in info[q]["preferred_name"]:
                    info[q]["preferred_name"] += " ; " + pref
                if pfam and pfam not in info[q]["pfams"]:
                    info[q]["pfams"] += " ; " + pfam

    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", required=True, type=Path)
    ap.add_argument("--fasta-in", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--keyword-file", default=None, type=Path)
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.keyword_file:
        kw_path = args.keyword_file
        print(f"[INFO] Using USER keyword file ONLY: {kw_path}")
    else:
        kw_path = write_default_keyword_file(args.outdir)
        print(f"[INFO] Using DEFAULT keyword file: {kw_path}")

    keywords = load_keywords_csv(kw_path)
    patterns = compile_patterns(keywords)

    egg = parse_emapper_annotations(args.annotations)

    # find which base IDs are housekeeping
    housekeeping = {}
    for bid, meta in egg.items():
        hay = f"{meta.get('description','')} {meta.get('preferred_name','')} {meta.get('pfams','')}"
        for kw, pat in patterns:
            if pat.search(hay):
                housekeeping[bid] = kw
                break

    if args.debug:
        print(f"[DEBUG] EggNOG base IDs parsed: {len(egg)}")
        print(f"[DEBUG] Housekeeping matches found: {len(housekeeping)}")

    fasta_out = unique_path(args.outdir / "ht_clusters.housekeeping_filtered.fasta")
    summary_out = unique_path(args.outdir / "housekeeping_filter_summary.tsv")

    kept = removed = 0
    mapped = 0

    with fasta_out.open("w", encoding="utf-8") as fout, summary_out.open("w", encoding="utf-8") as s:
        s.write("id_raw\tid_norm\tstatus\tmatched_keyword\tdescription\tpreferred_name\tpfams\n")

        for rec in SeqIO.parse(str(args.fasta_in), "fasta"):
            rid_raw = rec.id
            rid = normalize_id(rid_raw)

            if rid in egg:
                mapped += 1

            meta = egg.get(rid, {})
            desc = meta.get("description", "")
            pref = meta.get("preferred_name", "")
            pfam = meta.get("pfams", "")

            if rid in housekeeping:
                removed += 1
                kw = housekeeping[rid]
                s.write(f"{rid_raw}\t{rid}\tremoved\t{kw}\t{desc}\t{pref}\t{pfam}\n")
            else:
                kept += 1
                SeqIO.write(rec, fout, "fasta")
                s.write(f"{rid_raw}\t{rid}\tkept\t-\t{desc}\t{pref}\t{pfam}\n")

    print("[SUCCESS] Step 8 complete.")
    print(f"[INFO] Output FASTA  : {fasta_out}")
    print(f"[INFO] Summary TSV   : {summary_out}")
    print(f"[INFO] Kept          : {kept}")
    print(f"[INFO] Removed       : {removed}")
    if args.debug:
        print(f"[DEBUG] FASTA normalized IDs found in EggNOG: {mapped} / {kept+removed}")


if __name__ == "__main__":
    main()

