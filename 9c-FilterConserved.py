#!/usr/bin/env python3
import os
import argparse
import shutil
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description="Move conserved candidates to avoid unnecessary tree building.")
parser.add_argument("--summary", required=True, help="Summary file from Script 9b")
parser.add_argument("--homologs_dir", required=True, help="Directory containing the fasta files from Script 9b")
parser.add_argument("--rejected_dir", default="Result_HT/homologs_conserved_skipped", help="Where to move rejected files")
parser.add_argument("--max_plants", type=int, default=200, help="Max plant species allowed (approx 50% of 400)")
parser.add_argument("--max_fungi", type=int, default=500, help="Max fungi species allowed (approx 50% of 1000)")

args = parser.parse_args()

rejected_path = Path(args.rejected_dir)
rejected_path.mkdir(parents=True, exist_ok=True)
homologs_path = Path(args.homologs_dir)

print(f"[INFO] Reading summary: {args.summary}")
cols = ["filename", "qseqid", "status", "n_base", "n_sel", "n_added", "n_plants", "n_fungi", "n_other", "reason"]

moved_count = 0

try:
    with open(args.summary, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8: continue
            
            filename = parts[0]
            try:
                n_plants = int(parts[6])
                n_fungi = int(parts[7])
            except ValueError:
                continue

            if n_plants > args.max_plants and n_fungi > args.max_fungi:
                src_file = homologs_path / filename
                dst_file = rejected_path / filename
                
                if src_file.exists():
                    shutil.move(str(src_file), str(dst_file))
                    moved_count += 1

    print(f"[SUCCESS] Moved {moved_count} conserved files to {rejected_path}")
    print(f"[INFO] Only remaining files in {homologs_path} will be processed by Script 10.")

except Exception as e:
    print(f"[ERROR] {e}")
