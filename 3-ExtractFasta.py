#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


parser = argparse.ArgumentParser(
    description="Extract FUNGI (query) regions from filtered BLAST results."
)

parser.add_argument(
    "-i", "--input-tsv",
    required=True,
    help="Filtered BLAST TSV file (output of script 2)"
)

parser.add_argument(
    "--fungi-dir",
    required=True,
    help="Directory containing FUNGI genome FASTA files"
)

parser.add_argument(
    "-o", "--outdir",
    default="selected_sequences",
    help="Output directory (default: selected_sequences)"
)

args = parser.parse_args()

input_tsv = Path(args.input_tsv)
fungi_dir = Path(args.fungi_dir)
outdir = Path(args.outdir)

# =========================
# Checks
# =========================
if not input_tsv.exists():
    sys.exit(f"[ERROR] Input TSV not found: {input_tsv}")

if not fungi_dir.is_dir():
    sys.exit(f"[ERROR] Fungi directory not found: {fungi_dir}")

outdir.mkdir(parents=True, exist_ok=True)

# =========================
# Read TSV
# =========================
print(f"[INFO] Reading {input_tsv} ...")
try:
    df = pd.read_csv(input_tsv, sep="\t")
except Exception as e:
    sys.exit(f"[ERROR] Could not read TSV: {e}")

required_cols = {"fungi_genome", "qseqid", "qstart", "qend"}
missing = required_cols - set(df.columns)
if missing:
    sys.exit(f"[ERROR] Missing columns in TSV: {', '.join(sorted(missing))}")

# =========================
# Group by fungi genome
# =========================
grouped = df.groupby("fungi_genome")
total_extracted = 0

for fungi_genome, group in grouped:
    # Locate genome FASTA
    genome_path = fungi_dir / str(fungi_genome)
    if not genome_path.exists():
        genome_path = fungi_dir / f"{fungi_genome}.fasta"

    if not genome_path.exists():
        print(f"[WARNING] Genome FASTA not found for {fungi_genome}. Skipping.")
        continue

    print(f"[INFO] Processing {fungi_genome}")

    try:
        seq_dict = SeqIO.to_dict(SeqIO.parse(str(genome_path), "fasta"))
    except Exception as e:
        print(f"[ERROR] Failed to parse {genome_path}: {e}")
        continue

    out_fasta = outdir / f"selected_{fungi_genome}.fasta"
    extracted = 0

    with out_fasta.open("w") as out_f:
        for _, row in group.iterrows():
            scaffold = str(row["qseqid"])
            start = int(row["qstart"])
            end = int(row["qend"])

            if scaffold not in seq_dict:
                continue

            if start > end:
                start, end = end, start

            fragment = seq_dict[scaffold].seq[start - 1:end]
            header = f"{fungi_genome}__{scaffold}_{start}-{end}"

            out_f.write(f">{header}\n{fragment}\n")
            extracted += 1
            total_extracted += 1

    print(f"[INFO]  → {extracted} fragments written to {out_fasta}")

# =========================
# Summary
# =========================
print(f"[SUCCESS] Extraction finished")
print(f"[INFO] Total fragments extracted: {total_extracted}")
print(f"[INFO] Output directory: {outdir}")
