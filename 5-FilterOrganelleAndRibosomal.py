#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
from pathlib import Path
from typing import Set, Optional

from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def ensure_blast_db(fasta: Path, dbtype: str = "nucl") -> Path:
    """
    Ensure BLAST DB exists for a fasta file.
    If DB does not exist, create it using the same name/prefix as the fasta.
    """
    if dbtype == "nucl":
        required_ext = [".nhr", ".nin", ".nsq"]
    else:
        required_ext = [".phr", ".pin", ".psq"]

    missing = [ext for ext in required_ext if not Path(str(fasta) + ext).exists()]
    if not missing:
        return fasta

    print(
        f"[WARNING] BLAST database not found for {fasta}\n"
        f"          Missing files: {', '.join(missing)}\n"
        f"          -> Creating BLAST database now (this may take some time)."
    )

    try:
        subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", dbtype, "-out", str(fasta)],
            check=True,
        )
    except FileNotFoundError:
        raise SystemExit("[ERROR] makeblastdb not found (is BLAST+ installed and in PATH?)")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Failed to create BLAST DB for {fasta} (code {e.returncode})")

    return fasta


def run_blastn_hitlist(
    query_fasta: Path,
    db_fasta: Path,
    threads: int,
    evalue: str,
    description: str,
) -> Set[str]:
    """
    Run blastn and return set of query IDs that got at least one hit.
    Uses outfmt 6 qseqid and max_target_seqs 1 for speed.
    """
    print(f"[INFO] Running BLASTn vs {description} ...")

    db_prefix = ensure_blast_db(db_fasta, dbtype="nucl")

    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-outfmt", "6 qseqid",
        "-max_target_seqs", "1",
        "-evalue", str(evalue),
        "-num_threads", str(threads),
    ]

    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        raise SystemExit("[ERROR] blastn not found (is BLAST+ installed and in PATH?)")
    except subprocess.CalledProcessError as e:
        eprint(f"[ERROR] BLASTn failed vs {description}:")
        if e.stderr:
            eprint(e.stderr)
        raise SystemExit(1)

    hits = set()
    for line in res.stdout.splitlines():
        qid = line.strip()
        if qid:
            hits.add(qid)

    print(f"       -> Found {len(hits)} candidate hits in {description}.")
    return hits


def warn_about_missing_dbs(rDNA_db: Optional[Path], plast_db: Optional[Path]) -> None:
    if rDNA_db is None and plast_db is None:
        eprint(
            "\n[WARNING] No rDNA/rRNA database and no plastid/mitochondrial database were provided.\n"
            "          Script will run, BUT no filtering of rDNA/plastDNA will happen.\n"
            "          This can lead to a VERY HIGH number of candidates, mostly false positives,\n"
            "          and will considerably increase the runtime of the downstream pipeline.\n"
        )
    elif rDNA_db is None or plast_db is None:
        missing = "rDNA/rRNA" if rDNA_db is None else "plastid/mitochondrial"
        eprint(
            f"\n[WARNING] Only one database was provided. Missing: {missing} database.\n"
            "          Script will run, BUT filtering will be incomplete.\n"
            "          This can lead to elevated candidate counts (false positives) and increase\n"
            "          the runtime of the downstream pipeline.\n"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Script 5: Filter candidates by BLASTn against rDNA/rRNA and plastid/mitochondrial databases "
            "(TE filtering postponed)."
        )
    )

    parser.add_argument(
        "--fasta-in", required=True, type=Path,
        help="Candidate FASTA input (e.g., Result_HT/trf_clean/ht_candidates.cleaned.fasta)"
    )
    parser.add_argument(
        "--fasta-out", default="ht_candidates.filtered.fasta", type=Path,
        help="Filtered FASTA output (default: ht_candidates.filtered.fasta)"
    )
    parser.add_argument(
        "--summary", default="candidate_filter_summary.tsv", type=Path,
        help="TSV summary output (default: candidate_filter_summary.tsv)"
    )

    # Databases (optional)
    parser.add_argument(
        "--rDNA-db", default=None, type=Path,
        help="FASTA database for ribosomal DNA/RNA (rDNA/rRNA) (optional)"
    )
    parser.add_argument(
        "--plastDNA-db", default=None, type=Path,
        help="FASTA database for plastid and/or mitochondrial DNA (optional)"
    )

    # BLAST options
    parser.add_argument("--threads", type=int, default=8, help="BLAST threads (default: 8)")
    parser.add_argument("--evalue", default="1e-10", help="BLAST e-value threshold (default: 1e-10)")

    args = parser.parse_args()

    if not args.fasta_in.exists():
        sys.exit(f"[ERROR] fasta-in not found: {args.fasta_in}")

    # Validate DB FASTA paths: warn if provided but missing, then ignore
    rDNA_db: Optional[Path]
    plast_db: Optional[Path]

    if args.rDNA_db is not None and not args.rDNA_db.exists():
        eprint(
            f"[WARNING] rDNA database FASTA not found: {args.rDNA_db}\n"
            "          rDNA/rRNA filtering will be SKIPPED."
        )
        rDNA_db = None
    else:
        rDNA_db = args.rDNA_db

    if args.plastDNA_db is not None and not args.plastDNA_db.exists():
        eprint(
            f"[WARNING] plastid/mitochondrial database FASTA not found: {args.plastDNA_db}\n"
            "          plastid/mitochondrial filtering will be SKIPPED."
        )
        plast_db = None
    else:
        plast_db = args.plastDNA_db

    warn_about_missing_dbs(rDNA_db, plast_db)

    # Run BLAST checks (only if DB provided)
    rDNA_hits: Set[str] = set()
    plast_hits: Set[str] = set()

    if rDNA_db is not None:
        rDNA_hits = run_blastn_hitlist(
            query_fasta=args.fasta_in,
            db_fasta=rDNA_db,
            threads=args.threads,
            evalue=args.evalue,
            description="rDNA/rRNA database"
        )

    if plast_db is not None:
        plast_hits = run_blastn_hitlist(
            query_fasta=args.fasta_in,
            db_fasta=plast_db,
            threads=args.threads,
            evalue=args.evalue,
            description="plastid/mitochondrial database"
        )

    ids_to_remove = rDNA_hits.union(plast_hits)

    # Write outputs
    kept = 0
    removed = 0

    args.fasta_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    with args.fasta_out.open("w", encoding="utf-8") as fout, \
         args.summary.open("w", encoding="utf-8") as s:
        s.write("id\tstatus\treasons\n")

        for rec in SeqIO.parse(str(args.fasta_in), "fasta"):
            rid = rec.id
            if rid in ids_to_remove:
                reasons = []
                if rid in rDNA_hits:
                    reasons.append("rDNA_rRNA_hit")
                if rid in plast_hits:
                    reasons.append("plastid_mito_hit")
                s.write(f"{rid}\tremoved\t{'+'.join(reasons)}\n")
                removed += 1
            else:
                SeqIO.write(rec, fout, "fasta")
                s.write(f"{rid}\tkept\t-\n")
                kept += 1

    print(f"[SUCCESS] Filtered FASTA saved to: {args.fasta_out}")
    print(f"[INFO] Kept: {kept}")
    print(f"[INFO] Removed: {removed} (rDNA/plastDNA hits)")
    print(f"[INFO] Summary written to: {args.summary}")


if __name__ == "__main__":
    main()
