#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description="Extract final candidates and anonymize headers with mapping.")
parser.add_argument("-i", "--input_candidates", required=True, help="Candidate TSV file")
parser.add_argument("-g", "--genomes_dir", required=True, help="Directory containing Genome FASTA files (Source of candidates)")
parser.add_argument("-o", "--output", default="hgt_candidates.fasta", help="Output Multi-FASTA file")
parser.add_argument("--mapping_out", default="hgt_id_mapping.tsv", help="Output file for ID mapping")

args = parser.parse_args()

df = pd.read_csv(args.input_candidates, sep="\t")
required_cols = ["sseqid_fungi", "sstart_fungi", "send_fungi", "fungi_genome", "qseqid"]
if not all(col in df.columns for col in required_cols):
    exit(f"[ERROR] Input TSV missing columns. Required: {required_cols}")

print(f"[INFO] Extracting {len(df)} sequences...")

grouped = df.groupby("fungi_genome")
mapping_data = []
counter = 1

with open(args.output, "w") as out_f:
    for genome_name, group in grouped:
        genome_path = os.path.join(args.genomes_dir, str(genome_name) + ".fasta")
        if not os.path.exists(genome_path):
             genome_path = os.path.join(args.genomes_dir, str(genome_name)) 
        
        if not os.path.exists(genome_path):
            print(f"[WARNING] Genome {genome_name} not found. Skipping {len(group)} candidates.")
            continue
            
        print(f"  -> Processing {genome_name}...")
        seq_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        
        for _, row in group.iterrows():
            scaffold = str(row["sseqid_fungi"]).strip()
            start = int(row["sstart_fungi"])
            end = int(row["send_fungi"])
            query_ref = str(row["qseqid"])
            
            if start > end: start, end = end, start
            
            if scaffold in seq_dict:
                seq_record = seq_dict[scaffold]
                fragment = seq_record.seq[start-1:end]
                safe_id = f"CAND_{counter:05d}"
                clean_genome = str(genome_name).replace(" ", "_")
                original_info = f"{clean_genome}|{scaffold}|{start}-{end}|vs|{query_ref}"
                out_f.write(f">{safe_id}\n{fragment}\n")
                mapping_data.append({"safe_id": safe_id, "original_header": original_info})
                counter += 1
            else:
                print(f"[WARNING] Scaffold {scaffold} not found in {genome_name}")

map_df = pd.DataFrame(mapping_data)
map_df.to_csv(args.mapping_out, sep="\t", index=False)
print(f"[SUCCESS] Extracted sequences to {args.output}")
print(f"[SUCCESS] ID Mapping saved to {args.mapping_out}")
