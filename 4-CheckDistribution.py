#!/usr/bin/env python3
import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Check distribution: Reject if conserved in BOTH kingdoms (>20%).")
parser.add_argument("-s", "--selected_fasta", required=True, help="Sequences to check (FASTA)")
parser.add_argument("-p", "--plant_genomes", required=True, help="Directory of Plant Genomes")
parser.add_argument("-f", "--fungi_genomes", required=True, help="Directory of Fungi Genomes")
parser.add_argument("-o", "--output_summary", default="distribution_summary.tsv", help="Output summary")
parser.add_argument("-k", "--keep_list", default="candidates_patchy.txt", help="List of IDs to keep")
parser.add_argument("-t", "--threads", default=8, help="Threads")

args = parser.parse_args()

# Threshold
CONSERVATION_THRESHOLD = 0.10

def make_blast_db(fasta_dir, db_name):
    """Concatenates all genomes in a dir and makes a BLAST DB."""
    mega_fasta = f"{db_name}_all.fasta"
    if not os.path.exists(mega_fasta):
        print(f"[INFO] building {db_name} database (concatenating genomes)...")
        # This is heavy but necessary for global stats. 
        # Better approach: makeblastdb on list of files if blast+ supports it, 
        # or cat them.
        os.system(f"cat {fasta_dir}/*.fasta > {mega_fasta}")
        cmd = ["makeblastdb", "-in", mega_fasta, "-dbtype", "nucl", "-out", db_name]
        subprocess.run(cmd, check=True)
    return db_name

def run_blast_count(query, db, threads):
    """Runs blast and returns count of distinct subjects hit."""
    cmd = [
        "blastn", "-query", query, "-db", db, 
        "-outfmt", "6 qseqid sseqid", 
        "-num_threads", str(threads),
        "-evalue", "1e-5", "-perc_identity", "70",
        "-max_target_seqs", "5000" 
    ]
    process = subprocess.run(cmd, capture_output=True, text=True)

    hits = {}
    for line in process.stdout.strip().split("\n"):
        if not line: continue
        q, s = line.split("\t")
        if q not in hits: hits[q] = set()
        hits[q].add(s)
    
    return {k: len(v) for k, v in hits.items()}

# 1. Count Total Genomes
num_plants = len([f for f in os.listdir(args.plant_genomes) if f.endswith(".fasta")])
num_fungi = len([f for f in os.listdir(args.fungi_genomes) if f.endswith(".fasta")])

print(f"[INFO] Total Plants: {num_plants}, Total Fungi: {num_fungi}")

# 2. Make/Get DBs
plant_db = make_blast_db(args.plant_genomes, "temp_plant_db")
fungi_db = make_blast_db(args.fungi_genomes, "temp_fungi_db")

# 3. Blast
print("[INFO] Blasting against Plants...")
plant_counts = run_blast_count(args.selected_fasta, plant_db, args.threads)

print("[INFO] Blasting against Fungi...")
fungi_counts = run_blast_count(args.selected_fasta, fungi_db, args.threads)

# 4. Filter
kept_ids = []
with open(args.output_summary, "w") as f:
    f.write("qseqid\tplant_hits\tplant_pct\tfungi_hits\tfungi_pct\tstatus\n")
    
    for record in SeqIO.parse(args.selected_fasta, "fasta"):
        q = record.id
        p_hits = plant_counts.get(q, 0)
        f_hits = fungi_counts.get(q, 0)
        
        p_pct = p_hits / num_plants if num_plants > 0 else 0
        f_pct = f_hits / num_fungi if num_fungi > 0 else 0
        
        if p_pct > CONSERVATION_THRESHOLD and f_pct > CONSERVATION_THRESHOLD:
            status = "REJECT_CONSERVED"
        elif p_hits == 0 and f_hits == 0:
             status = "REJECT_NO_HIT"
        else:
            status = "KEEP_PATCHY"
            kept_ids.append(q)
            
        f.write(f"{q}\t{p_hits}\t{p_pct:.2f}\t{f_hits}\t{f_pct:.2f}\t{status}\n")

with open(args.keep_list, "w") as f:
    for k in kept_ids: f.write(k + "\n")

print(f"[SUCCESS] Kept {len(kept_ids)} patchy candidates.")
