#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract sequences from genomes (Bidirectional/Generic).")
parser.add_argument("-i", "--input_tsv", default="filtered_blast_results_with_fungi.tsv", help="Input filtered BLAST results")
parser.add_argument("-g", "--genomes_dir", required=True, help="Directory containing Genome FASTA files (Subject Genomes)")
parser.add_argument("-o", "--outdir", default="selected_sequences", help="Output directory")

args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

print(f"[INFO] Reading {args.input_tsv}...")
try:
    filtered_results = pd.read_csv(args.input_tsv, sep="\t")
    selected_ids = set(filtered_results["sseqid"].astype(str))
except Exception as e:
    exit(f"[ERROR] Could not read input TSV: {e}")

print(f"[INFO] Scanning genomes in {args.genomes_dir}...")

found_count = 0

for genome_file in os.listdir(args.genomes_dir):
    if genome_file.lower().endswith((".fasta", ".fa", ".fna")):
        fasta_path = os.path.join(args.genomes_dir, genome_file)
        species_name = os.path.splitext(genome_file)[0].replace(" ", "_").replace(".", "_")
        output_fasta = os.path.join(args.outdir, f"selected_{species_name}.fasta")
        hits_in_genome = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in selected_ids:
                original_id = record.id
                safe_id = original_id.replace("|", "_").replace(":", "_")
                new_id = f"{species_name}__{safe_id}"
                record.id = new_id
                record.description = new_id
                hits_in_genome.append(record)
                found_count += 1
        if hits_in_genome:
            with open(output_fasta, "w") as output_handle:
                SeqIO.write(hits_in_genome, output_handle, "fasta")

print(f"[SUCCESS] Extracted and renamed {found_count} sequences.")
