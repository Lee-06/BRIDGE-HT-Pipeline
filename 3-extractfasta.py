#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract SPECIFIC FRAGMENTS from genomes.")
parser.add_argument("-i", "--input_tsv", default="filtered_blast_results_with_fungi.tsv", help="Input filtered BLAST results")
parser.add_argument("-g", "--genomes_dir", required=True, help="Directory containing Genome FASTA files")
parser.add_argument("-o", "--outdir", default="selected_sequences", help="Output directory")

args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

print(f"[INFO] Reading {args.input_tsv}...")
try:
    df = pd.read_csv(args.input_tsv, sep="\t")
    if "plant_genome" not in df.columns:
        sys.exit("[ERROR] Input TSV missing 'plant_genome' column. Please run corrected Script 2.")
except Exception as e:
    exit(f"[ERROR] Could not read input TSV: {e}")

print(f"[INFO] Extracting fragments...")

grouped = df.groupby("plant_genome")
count = 0

for genome_name, group in grouped:
    genome_path = os.path.join(args.genomes_dir, str(genome_name) + ".fasta")
    if not os.path.exists(genome_path):
         genome_path = os.path.join(args.genomes_dir, str(genome_name))
    if not os.path.exists(genome_path):
        print(f"[WARNING] Genome file for {genome_name} not found. Skipping {len(group)} hits.")
        continue

    try:
        seq_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    except Exception as e:
        print(f"[ERROR] Failed to parse {genome_name}: {e}")
        continue

    out_fasta = os.path.join(args.outdir, f"selected_{genome_name}.fasta")
    
    with open(out_fasta, "w") as out_f:
        for _, row in group.iterrows():
            scaffold = str(row["sseqid"])
            start = int(row["sstart"])
            end = int(row["send"])
            
            if start > end:
                start, end = end, start
            
            if scaffold in seq_dict:
                full_seq = seq_dict[scaffold].seq
                fragment = full_seq[start-1 : end]
                safe_id = f"{genome_name}__{scaffold}_{start}-{end}"
                out_f.write(f">{safe_id}\n{fragment}\n")
                count += 1
            else:
                pass

print(f"[SUCCESS] Extracted {count} fragments to {args.outdir}/")
