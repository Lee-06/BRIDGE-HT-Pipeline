#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple


# TE signatures matched case-insensitively as substrings against the combined
# EggNOG Description + Preferred_name + PFAMs text. Covers free-text descriptions
# and common PFAM families EggNOG emits for mobile elements.
DEFAULT_TE_SIGNATURES = [
    # free-text descriptions
    "transposase", "transposon", "transposable element", "mobile element",
    "retrotransposon", "retroelement", "reverse transcriptase",
    "integrase", "gag-pol", "gag-pre", "polyprotein", "retrovirus",
    "helitron", "tyrosine recombinase",
    "ltr", "line-1", "non-ltr",
    "gypsy", "copia", "mariner", "piggybac", "harbinger", "cacta",
    "polinton", "maverick", "hat ", "mutator",
    # PFAM family tokens (EggNOG PFAMs column)
    "rve",              # Integrase core domain
    "rvt_",             # Reverse transcriptase (RVT_1/RVT_2/RVT_3)
    "dde_tnp",          # DDE transposase families
    "dde_",             # DDE superfamily
    "hth_tnp",          # transposase HTH
    "retrotrans_gag", "retrotransposon_gag",
    "transp_tc5", "dimer_tnp_hat", "zf-bed",
    "mule", "plant_tran",
    "rnase_h",
]


def run(cmd, cwd=None):
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Command failed ({e.returncode}): {' '.join(cmd)}")


def normalize_id(s: str) -> str:
    """Keep before first '|', strip trailing _<digits> (matches Script 8 behaviour)."""
    s = (s or "").strip()
    if "|" in s:
        s = s.split("|", 1)[0]
    if "_" in s:
        left, right = s.rsplit("_", 1)
        if right.isdigit():
            s = left
    return s


def parse_emapper_annotations(path: Path) -> Dict[str, Dict[str, str]]:
    """
    Parse EggNOG emapper.annotations using fixed column positions:
        0 query  7 Description  8 Preferred_name  19 PFAMs
    Multiple ORFs (_0/_1/...) collapse to the same normalized base id.
    """
    info: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 20:
                continue
            q = normalize_id(cols[0])
            desc = cols[7].strip()
            pref = cols[8].strip()
            pfam = cols[19].strip() if len(cols) > 19 else ""
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


def read_fasta_ids(path: Path) -> List[str]:
    ids: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def load_te_signatures(path: Path) -> List[str]:
    txt = path.read_text(encoding="utf-8", errors="replace")
    parts = [p.strip() for p in txt.replace("\n", ",").split(",")]
    sigs = [p for p in parts if p]
    seen, out = set(), []
    for s in sigs:
        k = s.casefold()
        if k not in seen:
            out.append(s)
            seen.add(k)
    return out


def match_te(meta: Dict[str, str], signatures: List[str]) -> Tuple[bool, str]:
    hay = " ".join([
        meta.get("description", ""),
        meta.get("preferred_name", ""),
        meta.get("pfams", ""),
    ]).casefold()
    for sig in signatures:
        if sig.casefold() in hay:
            return True, sig
    return False, ""


def label_tes(fasta: Path, annotations: Path, out_tsv: Path, signatures: List[str]) -> None:
    """Write a non-destructive TE/non_TE label TSV from EggNOG annotation."""
    egg = parse_emapper_annotations(annotations)
    print(f"[INFO] Parsed EggNOG annotations for {len(egg)} base IDs.")

    te_count = gene_count = unann_count = 0
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("candidate_id\tclassification\tmatched_signature\tdescription\tpreferred_name\tpfams\n")
        for raw_id in read_fasta_ids(fasta):
            bid = normalize_id(raw_id)
            meta = egg.get(bid)
            if meta is None:
                f.write(f"{raw_id}\tunannotated\t-\t-\t-\t-\n")
                unann_count += 1
                continue
            is_te, sig = match_te(meta, signatures)
            cls = "TE" if is_te else "non_TE"
            f.write(
                f"{raw_id}\t{cls}\t{sig or '-'}\t"
                f"{meta.get('description', '-') or '-'}\t"
                f"{meta.get('preferred_name', '-') or '-'}\t"
                f"{meta.get('pfams', '-') or '-'}\n"
            )
            if is_te:
                te_count += 1
            else:
                gene_count += 1

    print("[SUCCESS] TE labelling complete (no sequences removed).")
    print(f"  -> TE candidates labelled       : {te_count}")
    print(f"  -> non-TE candidates labelled   : {gene_count}")
    print(f"  -> unannotated (no EggNOG hit)  : {unann_count}")
    print(f"  -> Label table written to       : {out_tsv}")


def main():
    p = argparse.ArgumentParser(
        description=(
            "Script 7: Annotate clustered HT candidates (output of Script 6) with "
            "EggNOG-mapper, then label each candidate as TE or non_TE from the "
            "annotation (non-destructive: no sequence is removed)."
        )
    )

    # Script 6 output -> Script 7 input
    p.add_argument(
        "--clusters-fasta", required=True, type=Path,
        help="Clustered FASTA from Script 6 (e.g., Result_HT/ht_clusters.fasta)"
    )

    # EggNOG config
    p.add_argument(
        "--eggnog-data-dir", required=True, type=Path,
        help="EggNOG database directory (must contain eggnog_proteins.dmnd)"
    )
    p.add_argument(
        "--dmnd-db", default=None, type=Path,
        help="Optional path to eggnog_proteins.dmnd (default: <eggnog-data-dir>/eggnog_proteins.dmnd)"
    )

    # Output
    p.add_argument(
        "--outdir", default=".", type=Path,
        help="Output directory (default: current directory)"
    )
    p.add_argument(
        "--output-prefix", default="ht_annotations",
        help="EggNOG output prefix (default: ht_annotations)"
    )

    # emapper options
    p.add_argument("--cpu", type=int, default=8, help="CPUs for emapper (default: 8)")
    p.add_argument("--method", default="diamond", choices=["diamond", "hmmer"],
                   help="Search method (-m) (default: diamond)")
    p.add_argument("--itype", default="genome", help="Input type (--itype) (default: genome)")
    p.add_argument("--translate", action="store_true",
                   help="Translate nucleotide sequences (--translate). Recommended for genome DNA.")
    p.add_argument("--emapper", default="emapper.py", help="Path to emapper.py (default: emapper.py)")

    # TE labelling options
    p.add_argument(
        "--te-labels-out", default=None, type=Path,
        help="Path for the TE label TSV (default: <outdir>/te_labels.tsv)"
    )
    p.add_argument(
        "--te-keyword-file", default=None, type=Path,
        help="Optional CSV/newline file of TE signatures to use INSTEAD of the built-in defaults"
    )
    p.add_argument(
        "--skip-te-labelling", action="store_true",
        help="Skip the TE labelling step (annotation still runs normally)"
    )

    args = p.parse_args()

    # --- Checks ---
    if not args.clusters_fasta.exists():
        sys.exit(f"[ERROR] clusters-fasta not found: {args.clusters_fasta}")
    if not args.eggnog_data_dir.is_dir():
        sys.exit(f"[ERROR] eggnog-data-dir not found: {args.eggnog_data_dir}")

    dmnd_db = args.dmnd_db if args.dmnd_db else (args.eggnog_data_dir / "eggnog_proteins.dmnd")
    if not dmnd_db.exists():
        sys.exit(f"[ERROR] eggnog_proteins.dmnd not found: {dmnd_db}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Resolve to absolute paths (emapper runs with cwd=outdir)
    clusters_fasta_abs = args.clusters_fasta.resolve()
    eggnog_dir_abs = args.eggnog_data_dir.resolve()
    dmnd_db_abs = dmnd_db.resolve()

    # --- Run EggNOG-mapper ---
    cmd = [
        args.emapper,
        "-i", str(clusters_fasta_abs),
        "--output", args.output_prefix,
        "--cpu", str(args.cpu),
        "--dmnd_db", str(dmnd_db_abs),
        "-m", args.method,
        "--data_dir", str(eggnog_dir_abs),
        "--itype", args.itype,
    ]
    if args.translate:
        cmd.append("--translate")

    print("[INFO] Running EggNOG-mapper (Script 7) ...")
    print("       " + " ".join(cmd))
    run(cmd, cwd=str(args.outdir))
    print("[SUCCESS] EggNOG annotation complete.")
    print(f"[INFO] Output prefix : {args.output_prefix}")
    print(f"[INFO] Output directory : {args.outdir}")

    # --- TE labelling ---
    if args.skip_te_labelling:
        print("[INFO] TE labelling skipped (--skip-te-labelling).")
        return

    annotations_file = args.outdir / f"{args.output_prefix}.emapper.annotations"
    if not annotations_file.exists():
        print(f"[WARNING] Annotation file not found, skipping TE labelling: {annotations_file}")
        return

    if args.te_keyword_file:
        signatures = load_te_signatures(args.te_keyword_file)
        print(f"[INFO] Using USER TE signatures: {args.te_keyword_file} ({len(signatures)} terms)")
    else:
        signatures = DEFAULT_TE_SIGNATURES
        print(f"[INFO] Using DEFAULT TE signatures ({len(signatures)} terms)")

    te_out = args.te_labels_out if args.te_labels_out else (args.outdir / "te_labels.tsv")
    label_tes(clusters_fasta_abs, annotations_file, te_out, signatures)


if __name__ == "__main__":
    main()
