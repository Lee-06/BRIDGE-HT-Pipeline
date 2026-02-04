#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
Description: Filters raw BLAST results based on identity, alignment length, 
and scaffold length (for BOTH Query and Subject) to identify potential HT events.
"""
import os
import sys
import argparse
import pandas as pd

# ─── ARGUMENT PARSING ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Filter BLAST results for HT candidates.")
parser.add_argument("--blast_dir", required=True, help="Directory containing .blast output files")
parser.add_argument("--fungi_fai", required=True, help="Directory containing Fungi .fasta.fai index files (Query genomes)")
parser.add_argument("--plant_fai", required=True, help="Directory containing Plant .fasta.fai index files (Subject/DB genomes)")
parser.add_argument("--output", default="filtered_blast_results_with_fungi.tsv", help="Output TSV filename")

args = parser.parse_args()

BLAST_DIR = args.blast_dir
FUNGI_FAI_DIR = args.fungi_fai
PLANT_FAI_DIR = args.plant_fai
OUTPUT_FILE = args.output

# ─── CONFIGURATION ──────────────────────────────────────────────────────────────
# Thresholds defined in the study context, you might want different ones
IDENTITY_THRESHOLD = 80
ALIGNMENT_LENGTH_THRESHOLD = 500
SCAFFOLD_LENGTH_THRESHOLD = 20000  # 20 kb filter applied to both Query and Subject scaffolds

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# ─── HELPER FUNCTIONS ───────────────────────────────────────────────────────────
def load_fai(fai_path):
    """Parses a .fai file and returns a dictionary of {seq_id: length}."""
    lengths = {}
    if os.path.exists(fai_path):
        try:
            with open(fai_path, 'r') as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        lengths[parts[0]] = int(parts[1])
        except Exception as e:
            print(f"[WARNING] Could not read FAI file {fai_path}: {e}", file=sys.stderr)
    return lengths

# ─── MAIN EXECUTION ─────────────────────────────────────────────────────────────
def main():
    if not os.path.isdir(BLAST_DIR):
        sys.exit(f"[ERROR] BLAST directory not found: {BLAST_DIR}")

    print(f"[INFO] Starting filtering process...")
    print(f"[INFO] Thresholds: Identity>={IDENTITY_THRESHOLD}%, Len>={ALIGNMENT_LENGTH_THRESHOLD}bp")
    print(f"[INFO] Scaffold Length Threshold: >={SCAFFOLD_LENGTH_THRESHOLD}bp (Applied to BOTH Query and Subject)")

    filtered_frames = []
    file_count = 0

    for filename in os.listdir(BLAST_DIR):
        if not filename.endswith(".blast"):
            continue
            
        file_path = os.path.join(BLAST_DIR, filename)
        base = filename[:-6]
        if "_VS_" in base:
            subject_name, query_name = base.split("_VS_", 1)
        else:
            print(f"[WARNING] Filename {filename} does not match format DB_VS_QUERY. Skipping.")
            continue

        # 1. Load Query Index (Fungi)
        query_fai_path = os.path.join(FUNGI_FAI_DIR, query_name + ".fasta.fai")
        if not os.path.exists(query_fai_path):
             query_fai_path = os.path.join(FUNGI_FAI_DIR, query_name + ".fai") # try alt extension
        
        query_lengths = load_fai(query_fai_path)

        # 2. Load Subject Index (Plant)
        subject_fai_path = os.path.join(PLANT_FAI_DIR, subject_name + ".fasta.fai")
        if not os.path.exists(subject_fai_path):
             subject_fai_path = os.path.join(PLANT_FAI_DIR, subject_name + ".fai")

        subject_lengths = load_fai(subject_fai_path)
        
        if not query_lengths:
            # print(f"[DEBUG] Missing index for Query {query_name}. Scaffold filter might fail.", file=sys.stderr)
            pass
        if not subject_lengths:
            # print(f"[DEBUG] Missing index for Subject {subject_name}. Scaffold filter might fail.", file=sys.stderr)
            pass

        try:
            if os.path.getsize(file_path) == 0:
                continue

            blast_results = pd.read_csv(file_path, sep="\t", names=COLUMNS)
            
            # Tag the genome source (useful for later steps)
            blast_results["fungi_genome"] = query_name
            
            # Apply Filters
            # ---------------------------------------------------------
            # Filter 1: Basic BLAST metrics
            mask_basic = (blast_results["pident"] >= IDENTITY_THRESHOLD) & \
                         (blast_results["length"] >= ALIGNMENT_LENGTH_THRESHOLD)
            
            temp_filtered = blast_results[mask_basic]

            if temp_filtered.empty:
                continue

            # Filter 2: Scaffold Lengths (CRITICAL FIX)
            # Check QSEQID (Query/Fungi) length
            valid_query = temp_filtered["qseqid"].map(lambda x: query_lengths.get(str(x), 0)) >= SCAFFOLD_LENGTH_THRESHOLD
            
            # Check SSEQID (Subject/Plant) length
            valid_subject = temp_filtered["sseqid"].map(lambda x: subject_lengths.get(str(x), 0)) >= SCAFFOLD_LENGTH_THRESHOLD
            
            # Keep only rows satisfying BOTH
            temp_filtered = temp_filtered[valid_query & valid_subject]

            if not temp_filtered.empty:
                filtered_frames.append(temp_filtered)
        
        except Exception as e:
            print(f"[ERROR] Processing {filename}: {e}", file=sys.stderr)
        
        file_count += 1
        if file_count % 1000 == 0:
            print(f"[INFO] Processed {file_count} files...")

    # Concatenate all results
    if filtered_frames:
        final_df = pd.concat(filtered_frames, ignore_index=True)
        final_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
        print(f"[SUCCESS] Filtered results saved to {OUTPUT_FILE}")
        print(f"[INFO] Total hits kept: {len(final_df)}")
    else:
        print("[WARNING] No hits passed the filters. Output file not created.")

if __name__ == "__main__":
    main()
