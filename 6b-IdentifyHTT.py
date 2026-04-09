#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
from pathlib import Path
from typing import Set
from Bio import SeqIO

def ensure_blast_db(fasta: Path) -> Path:
    """Ensure BLAST DB exists for a fasta file."""
    required_ext = [".nhr", ".nin", ".nsq"]
    missing = [ext for ext in required_ext if not Path(str(fasta) + ext).exists()]
    if not missing:
        return fasta

    print(f"[INFO] Creating BLAST database for {fasta} ...")
    try:
        subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(fasta)],
            check=True, stdout=subprocess.DEVNULL
        )
    except Exception as e:
        sys.exit(f"[ERROR] Failed to create BLAST DB for {fasta}: {e}")
    return fasta

def main():
    parser = argparse.ArgumentParser(
        description="Script 6b: Identify Horizontal Transposon Transfers (HTT) by cross-referencing Repbase and flagging them in a TSV."
    )
    
    parser.add_argument("--fasta-in", required=True, type=Path,
                        help="Clustered candidates FASTA (output of Script 6)")
    parser.add_argument("--repbase-db", required=True, type=Path,
                        help="Repbase nucleotide FASTA database")
    parser.add_argument("--summary", default="htt_identification_summary.tsv", type=Path,
                        help="Summary TSV logging TE classification")
    
    # BLAST params
    parser.add_argument("--threads", type=int, default=8, help="BLAST threads (default: 8)")
    parser.add_argument("--evalue", default="1e-10", help="BLAST e-value threshold (default: 1e-10)")
    parser.add_argument("--min-pident", type=float, default=70.0, help="Minimum identity % for TE hit (default: 70)")
    
    args = parser.parse_args()

    if not args.fasta_in.exists():
        sys.exit(f"[ERROR] Input FASTA not found: {args.fasta_in}")
    if not args.repbase_db.exists():
        sys.exit(f"[ERROR] Repbase database not found: {args.repbase_db}")

    ensure_blast_db(args.repbase_db)

    print("[INFO] Running BLASTn against Repbase...")
    cmd = [
        "blastn",
        "-query", str(args.fasta_in),
        "-db", str(args.repbase_db),
        "-outfmt", "6 qseqid pident",
        "-max_target_seqs", "1",
        "-evalue", str(args.evalue),
        "-num_threads", str(args.threads),
    ]

    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] BLASTn failed: {e.stderr}")

    # Parse hits
    te_hits: Set[str] = set()
    for line in res.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            qseqid = parts[0]
            pident = float(parts[1])
            if pident >= args.min_pident:
                te_hits.add(qseqid)

    print(f"[INFO] Identified {len(te_hits)} sequences as Transposable Elements (HTT).")

    # Write Output TSV
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    
    genes_count, tes_count = 0, 0
    with args.summary.open("w", encoding="utf-8") as f_sum:
        f_sum.write("candidate_id\tclassification\n")
        
        # We iterate through the original FASTA to log EVERY candidate
        for rec in SeqIO.parse(str(args.fasta_in), "fasta"):
            if rec.id in te_hits:
                f_sum.write(f"{rec.id}\tHTT_Transposon\n")
                tes_count += 1
            else:
                f_sum.write(f"{rec.id}\tHGT_Gene\n")
                genes_count += 1

    print("[SUCCESS] HTT Flagging Complete.")
    print(f"  -> HGT standard genes flagged: {genes_count}")
    print(f"  -> HTT transposable elements flagged: {tes_count}")
    print(f"  -> Summary written to: {args.summary}")
    print(f"[NOTE] Please pass the ORIGINAL {args.fasta_in.name} to Step 7 (EggNOG).")

if __name__ == "__main__":
    main()
