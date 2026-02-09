#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description="Extract candidates, filter Ns, and anonymize.")
parser.add_argument("-i", "--input_candidates", required=True, help="Candidate TSV file")
parser.add_argument("-g", "--genomes_dir", required=True, help="Directory containing Source Genome FASTAs")
parser.add_argument("-o", "--output", default="ht_candidates.fasta", help="Output Multi-FASTA file")
parser.add_argument("--mapping_out", default="ht_id_mapping.tsv", help="Output file for ID mapping")

args = parser.parse_args()

def check_tandem_repeats(sequence, temp_id="temp_seq"):
    """
    Runs TRF (Tandem Repeats Finder) on a sequence string.
    Returns True if significant repeats are found (>50% of seq covered).
    """
    temp_fasta = f"{temp_id}.fasta"
    with open(temp_fasta, "w") as f:
        f.write(f">{temp_id}\n{sequence}\n")

    cmd = ["trf", temp_fasta, "2", "7", "7", "80", "10", "50", "500", "-h", "-ngs"]
    
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        dat_file = f"{temp_fasta}.2.7.7.80.10.50.500.dat"
        if os.path.exists(dat_file):
            with open(dat_file, "r") as df:
                lines = df.readlines()
                data_lines = [l for l in lines if len(l.split()) > 10]
                
            os.remove(dat_file)
            os.remove(temp_fasta)
            
            if len(data_lines) > 0:
                return True
            else:
                return False
        else:
            if os.path.exists(temp_fasta): os.remove(temp_fasta)
            return False

    except Exception:
        if os.path.exists(temp_fasta): os.remove(temp_fasta)
        return False

df = pd.read_csv(args.input_candidates, sep="\t")

if "sseqid_fungi" in df.columns:
    col_id, col_start, col_end, col_genome = "sseqid_fungi", "sstart_fungi", "send_fungi", "fungi_genome"
else:
    col_id, col_start, col_end, col_genome = "sseqid", "sstart", "send", "plant_genome"

print(f"[INFO] Processing {len(df)} candidates...")

grouped = df.groupby(col_genome)
mapping_data = []
counter = 1
kept_count = 0
rejected_count = 0

with open(args.output, "w") as out_f:
    for genome_name, group in grouped:
        genome_path = os.path.join(args.genomes_dir, str(genome_name) + ".fasta")
        if not os.path.exists(genome_path): genome_path = os.path.join(args.genomes_dir, str(genome_name))
        
        if not os.path.exists(genome_path):
            print(f"[WARNING] Genome {genome_name} not found.")
            continue
            
        seq_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        
        for _, row in group.iterrows():
            scaffold = str(row[col_id]).strip()
            start = int(row[col_start])
            end = int(row[col_end])
            if start > end: start, end = end, start
            
            if scaffold in seq_dict:
                full_scaffold_seq = seq_dict[scaffold].seq
                fragment = full_scaffold_seq[start-1:end]
                flank_upstream_start = max(0, start-1 - 5000)
                flank_upstream = full_scaffold_seq[flank_upstream_start : start-1]
                flank_downstream_end = min(len(full_scaffold_seq), end + 5000)
                flank_downstream = full_scaffold_seq[end : flank_downstream_end]

                if 'n' in fragment.lower() or 'n' in flank_upstream.lower() or 'n' in flank_downstream.lower():
                    rejected_count += 1
                    continue

                if check_tandem_repeats(str(fragment), f"cand_{counter}"):
                    print(f"    [REJECT] Candidate {counter} is repetitive.")
                    rejected_count += 1
                    continue

                safe_id = f"CAND_{counter:05d}"
                clean_genome = str(genome_name).replace(" ", "_")
                original_info = f"{clean_genome}|{scaffold}|{start}-{end}"
                
                out_f.write(f">{safe_id}\n{fragment}\n")
                mapping_data.append({"safe_id": safe_id, "original_header": original_info})
                counter += 1
                kept_count += 1

map_df = pd.DataFrame(mapping_data)
map_df.to_csv(args.mapping_out, sep="\t", index=False)
print(f"[SUCCESS] Kept {kept_count} candidates. Rejected {rejected_count} due to 'N's.")
print(f"[SUCCESS] Mapping saved to {args.mapping_out}")
