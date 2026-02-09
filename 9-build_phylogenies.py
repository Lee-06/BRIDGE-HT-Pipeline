#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Build Trees and Restore Headers.")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file (output of Script 6)")
parser.add_argument("-m", "--mapping", required=True, help="Mapping TSV file (output of Script 6)")
parser.add_argument("-db", "--database", required=True, help="Path to local BLAST database (nt)")
parser.add_argument("-o", "--outdir", default="phylogenies", help="Output directory")
parser.add_argument("-t", "--threads", default=4, help="Number of threads")
parser.add_argument("--max_hits", default=50, type=int, help="Max homologs")

args = parser.parse_args()

def check_tool(name):
    from shutil import which
    if which(name) is None:
        sys.exit(f"Error: {name} is not installed or not in your PATH.")

check_tool("mafft")
check_tool("iqtree")
check_tool("blastn")
check_tool("trimal")

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# Load Mapping: SafeID -> Full Name
print(f"[INFO] Loading ID mapping from {args.mapping}...")
map_df = pd.read_csv(args.mapping, sep="\t")
id_map = pd.Series(map_df.original_header.values, index=map_df.safe_id).to_dict()

def run_blast_local(query_seq, db_path, threads):
    """Fetches homologs and ensures headers are clean/readable."""
    cmd = [
        "blastn", "-db", db_path, "-query", "-", 
        "-outfmt", "6 sseqid sseq sscinames", 
        "-max_target_seqs", str(args.max_hits), 
        "-num_threads", str(threads)
    ]
    
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=query_seq)
    
    hits = []
    if process.returncode != 0: return hits
    
    for line in stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                acc, seq, sciname = parts[0], parts[1], parts[2]
                
                clean_name = sciname.replace(" ", "_").replace(":", "").replace("(", "").replace(")", "")
                clean_acc = acc.replace("|", "_")
                
                # Readable Homolog ID
                full_header = f"{clean_acc}_{clean_name}"
                hits.append((full_header, seq))
    return hits

def restore_tree_headers(tree_path, mapping_dict):
    """Replaces SafeIDs with Original Biological Headers in the tree file."""
    try:
        with open(tree_path, 'r') as f:
            tree_data = f.read()
        
        for safe_id, original_header in mapping_dict.items():
            if safe_id in tree_data:
                # Sanitize original header for Newick format (no colons/parentheses)
                clean_original = original_header.replace(":", "_").replace("(", "_").replace(")", "_")
                tree_data = tree_data.replace(safe_id, clean_original)
        
        with open(tree_path, 'w') as f:
            f.write(tree_data)
        return True
    except Exception as e:
        print(f"    [ERROR] Could not restore headers: {e}")
        return False

# Main Execution
for record in SeqIO.parse(args.input, "fasta"):
    safe_id = record.id
    print(f"--> Processing: {safe_id}")
    
    # 1. Fetch Homologs
    homologs = run_blast_local(str(record.seq), args.database, args.threads)
    
    if not homologs:
        print(f"    No homologs. Skipping.")
        continue

    # 2. Write Unaligned Fasta (SafeID + Homologs)
    fasta_path = os.path.join(args.outdir, f"{safe_id}_homologs.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{safe_id}\n{str(record.seq)}\n")
        seen = set()
        for h_head, h_seq in homologs:
            if h_head not in seen:
                f.write(f">{h_head}\n{h_seq}\n")
                seen.add(h_head)
    
    # 3. Align (MAFFT)
    aln_path = os.path.join(args.outdir, f"{safe_id}.aln")
    subprocess.run(["mafft", "--thread", str(args.threads), "--auto", fasta_path], stdout=open(aln_path, "w"), stderr=subprocess.DEVNULL)
    
    # 4. Trim (TrimAl)
    trimmed_path = os.path.join(args.outdir, f"{safe_id}.trimmed.aln")
    subprocess.run(["trimal", "-in", aln_path, "-out", trimmed_path, "-automated1"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # 5. Tree (IQ-TREE)
    subprocess.run(["iqtree2", "-s", trimmed_path, "-bb", "1000", "-nt", str(args.threads), "-quiet"])
    
    # 6. RESTORE HEADERS
    tree_file = trimmed_path + ".treefile"
    if os.path.exists(tree_file):
        restore_tree_headers(tree_file, id_map)
        print(f"    [SUCCESS] Tree built & headers restored: {tree_file}")
    else:
        print("    [FAIL] Tree not generated.")

print("Done.")
