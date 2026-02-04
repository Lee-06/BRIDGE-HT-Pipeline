#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter out housekeeping genes (EggNOG) and rRNA (SILVA), while checking for TEs (Repbase).")
parser.add_argument("--annotations", required=True, help="EggNOG annotation file (.annotations)")
parser.add_argument("--fasta_in", required=True, help="Clustered FASTA input (hgt_clusters.fasta)")
parser.add_argument("--fasta_out", default="hgt_filtered.fasta", help="Final Filtered FASTA output")
parser.add_argument("--silva", required=True, help="Path to local SILVA BLAST database")
parser.add_argument("--repbase", required=True, help="Path to local Repbase BLAST database")
parser.add_argument("--threads", default=4, help="Number of threads for BLAST")

args = parser.parse_args()

# 1. Keywords to exclude (EggNOG based)
HOUSEKEEPING_KEYWORDS = [
    "ribosomal", "18s", "28s", "5s", "rrna", "rdna", "ribonucleoprotein",
    "translation elongation factor", "mitochondrion", "mitochondrial",
    "cytochrome", "cox1", "nad", "atp6", "atp9", "chloroplast", "plastid",
    "glycolysis", "atp synthase", "nadh dehydrogenase", "oxidoreductase",
    "succinate dehydrogenase", "malate dehydrogenase", "dna polymerase",
    "rna polymerase", "helicase", "topoisomerase", "exonuclease", "ligase",
    "primase", "actin", "tubulin", "kinesin", "dynein", "myosin", "chaperone",
    "heat shock protein", "ubiquitin", "kinase", "phosphatase"
]
keywords_lower = [kw.lower() for kw in HOUSEKEEPING_KEYWORDS]

# --- Helper Function: Run BLASTn ---
def run_blast_check(fasta_input, db_path, threads, description):
    """Runs blastn and returns a set of hit IDs."""
    print(f"[INFO] Running BLASTn against {description} ({db_path})...")
    
    # Check if DB exists
    if not os.path.exists(db_path + ".nhr") and not os.path.exists(db_path + ".nsq"):
         # Try without extension or warn
         pass # BLAST might handle the base name correctly if configured
         
    cmd = [
        "blastn", 
        "-query", fasta_input,
        "-db", db_path,
        "-outfmt", "6 qseqid",
        "-max_target_seqs", "1",
        "-evalue", "1e-10",
        "-num_threads", str(threads)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        hits = set()
        for line in result.stdout.strip().split("\n"):
            if line:
                hits.add(line.split("\t")[0])
        print(f"       -> Found {len(hits)} hits in {description}.")
        return hits
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] BLAST against {description} failed: {e.stderr}")
        sys.exit(1)

# --- MAIN ---

# Step 1: Parse EggNOG Annotations
print("[INFO] Parsing EggNOG annotations...")
try:
    df = pd.read_csv(args.annotations, sep="\t", comment="#", header=None)
    eggnog_ids = set()
    for index, row in df.iterrows():
        row_text = " ".join(row.astype(str)).lower()
        if any(kw in row_text for kw in keywords_lower):
            eggnog_ids.add(str(row[0]))
    print(f"       -> Found {len(eggnog_ids)} housekeeping candidates via EggNOG.")
except Exception as e:
    sys.exit(f"[ERROR] Could not read annotations: {e}")

# Step 2: Run BLAST vs SILVA (rRNA removal)
silva_ids = run_blast_check(args.fasta_in, args.silva, args.threads, "SILVA (rRNA)")

# Step 3: Run BLAST vs Repbase
repbase_ids = run_blast_check(args.fasta_in, args.repbase, args.threads, "Repbase (Transposons)")

# Step 4: Combine Removal Lists
ids_to_remove = eggnog_ids.union(silva_ids)

te_kept = repbase_ids - ids_to_remove
print(f"[INFO] Verified {len(te_kept)} Transposable Elements via Repbase (these are RETAINED).")

# Step 5: Write Filtered FASTA
kept_count = 0
excluded_count = 0

with open(args.fasta_out, "w") as out_handle:
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        if record.id not in ids_to_remove:
            SeqIO.write(record, out_handle, "fasta")
            kept_count += 1
        else:
            excluded_count += 1

print(f"[SUCCESS] Filtered FASTA saved to {args.fasta_out}")
print(f"         Total Kept: {kept_count}")
print(f"         Removed: {excluded_count} (EggNOG + SILVA)")
