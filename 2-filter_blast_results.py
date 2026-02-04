#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Filter BLAST results (Identity, Alignment, Bilateral Scaffold Limit).")
parser.add_argument("--blast_dir", required=True, help="Directory containing .blast output files")
parser.add_argument("--fungi_fai", required=True, help="Directory containing Fungi .fasta.fai index files")
parser.add_argument("--plant_fai", required=True, help="Directory containing Plant .fasta.fai index files")
parser.add_argument("--output", default="filtered_blast_results_with_fungi.tsv", help="Output TSV filename")

args = parser.parse_args()

IDENTITY_THRESHOLD = 80
ALIGNMENT_LENGTH_THRESHOLD = 500
SCAFFOLD_LENGTH_THRESHOLD = 20000 

def load_fai(fai_path):
    lengths = {}
    if os.path.exists(fai_path):
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
    return lengths

filtered_frames = []

for filename in os.listdir(args.blast_dir):
    if not filename.endswith(".blast"): continue
    base = filename[:-6]
    if "_VS_" in base:
        subject_name, query_name = base.split("_VS_", 1)
    else:
        continue

    q_fai = os.path.join(args.fungi_fai, query_name + ".fasta.fai")
    if not os.path.exists(q_fai): q_fai = os.path.join(args.fungi_fai, query_name + ".fai")
    q_len = load_fai(q_fai)

    s_fai = os.path.join(args.plant_fai, subject_name + ".fasta.fai")
    if not os.path.exists(s_fai): s_fai = os.path.join(args.plant_fai, subject_name + ".fai")
    s_len = load_fai(s_fai)

    try:
        if os.path.getsize(os.path.join(args.blast_dir, filename)) == 0: continue
        
        df = pd.read_csv(os.path.join(args.blast_dir, filename), sep="\t", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
        
        df["fungi_genome"] = query_name
        df["plant_genome"] = subject_name
        df = df[(df["pident"] >= IDENTITY_THRESHOLD) & (df["length"] >= ALIGNMENT_LENGTH_THRESHOLD)]
        
        if df.empty: continue

        valid_q = df["qseqid"].map(lambda x: q_len.get(str(x), 0)) >= SCAFFOLD_LENGTH_THRESHOLD
        valid_s = df["sseqid"].map(lambda x: s_len.get(str(x), 0)) >= SCAFFOLD_LENGTH_THRESHOLD
        
        df = df[valid_q & valid_s]

        if not df.empty:
            filtered_frames.append(df)
            
    except Exception as e:
        print(f"[ERROR] {filename}: {e}", file=sys.stderr)

if filtered_frames:
    final = pd.concat(filtered_frames, ignore_index=True)
    final.to_csv(args.output, sep="\t", index=False)
    print(f"[SUCCESS] Saved {len(final)} hits to {args.output}")
else:
    print("[WARNING] No hits passed filtering.")
