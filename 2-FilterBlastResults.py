#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2-FilterBlastResults.py

Filters BLAST tabular (outfmt 6) results based on:
- percent identity
- alignment length (bp)
- scaffold/chr length threshold for BOTH Query (qseqid) and Subject (sseqid)

Also:
- Creates missing .fai indexes automatically using `samtools faidx` if --build-fai is set
- Adds provenance columns: fungi_genome (query file), plant_genome (subject file)

Expected BLAST filenames by default:
  <SUBJECT>_VS_<QUERY>.blast
(where SUBJECT = plant genome DB, QUERY = fungi genome query)

Output: one TSV combining all filtered hits.
"""

import argparse
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


DEFAULT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run(cmd: List[str]) -> None:
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Command failed ({e.returncode}): {' '.join(cmd)}")


def ensure_fai(fasta_path: Path, samtools: str = "samtools") -> Path:
    """
    Ensure <fasta>.fai exists. Create it with samtools faidx if missing.
    Returns the .fai path.
    """
    fai_path = Path(str(fasta_path) + ".fai")
    if fai_path.exists():
        return fai_path

    # Create .fai
    run([samtools, "faidx", str(fasta_path)])

    if not fai_path.exists():
        raise SystemExit(f"[ERROR] Failed to create FAI for {fasta_path}")
    return fai_path


def load_fai(fai_path: Path) -> Dict[str, int]:
    """
    Load FAI lengths: dict[seqname] = length
    """
    lengths: Dict[str, int] = {}
    try:
        with fai_path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        lengths[parts[0]] = int(parts[1])
                    except ValueError:
                        continue
    except Exception as e:
        eprint(f"[WARNING] Could not read FAI file {fai_path}: {e}")
    return lengths


def parse_subject_query(fname: str, blast_suffix: str) -> Optional[tuple[str, str]]:
    """
    From <SUBJECT>_VS_<QUERY><suffix> return (SUBJECT, QUERY).
    """
    if not fname.endswith(blast_suffix):
        return None
    base = fname[: -len(blast_suffix)]
    if "_VS_" not in base:
        return None
    subject, query = base.split("_VS_", 1)
    return subject, query


def find_fasta(genome_dir: Path, genome_name: str, exts: List[str]) -> Optional[Path]:
    """
    Find a FASTA file in genome_dir matching genome_name.
    genome_name may already include extension (e.g., X.fasta).
    We try:
      - exact match genome_dir/genome_name
      - if genome_name has no recognized ext: try adding each ext
    """
    exact = genome_dir / genome_name
    if exact.exists():
        return exact

    # If no extension in name, try adding known extensions
    if not any(genome_name.endswith(ext) for ext in exts):
        for ext in exts:
            candidate = genome_dir / f"{genome_name}{ext}"
            if candidate.exists():
                return candidate

    return None


def main():
    p = argparse.ArgumentParser(description="Filter BLAST results for HT candidates + auto-create .fai with samtools faidx.")

    # Required I/O
    p.add_argument("--blast-dir", required=True, type=Path,
                   help="Directory containing BLAST output files (e.g., Result_HT/blast)")
    p.add_argument("--fungi-dir", required=True, type=Path,
                   help="Directory containing fungi FASTA files (queries)")
    p.add_argument("--plant-dir", required=True, type=Path,
                   help="Directory containing plant FASTA files (subjects / DB)")
    p.add_argument("--out", default="filtered_blast_results.tsv", type=Path,
                   help="Output TSV file path")

    # Filters as CLI options
    p.add_argument("--identity", type=float, default=70.0,
                   help="Minimum percent identity (default: 70)")
    p.add_argument("--min-align-len", type=int, default=500,
                   help="Minimum alignment length (bp) (default: 500)")
    p.add_argument("--min-scaffold-len", type=int, default=20000,
                   help="Minimum scaffold/chr length (bp) for BOTH qseqid and sseqid (default: 20000)")

    # BLAST parsing options
    p.add_argument("--blast-suffix", default=".blast",
                   help="BLAST filename suffix to include (default: .blast)")
    p.add_argument("--columns", default=",".join(DEFAULT_COLUMNS),
                   help="Comma-separated columns for BLAST tabular (default matches outfmt 6 used in pipeline)")

    # FASTA extension handling
    p.add_argument("--extensions", default=".fasta,.fa,.fna",
                   help="Comma-separated FASTA extensions to try (default: .fasta,.fa,.fna)")

    # FAI building options
    p.add_argument("--build-fai", action="store_true",
                   help="If set, create missing .fai using samtools faidx")
    p.add_argument("--samtools", default="samtools",
                   help="samtools executable (default: samtools)")
    p.add_argument("--require-fai", action="store_true",
                   help="Fail if .fai cannot be found/created (otherwise missing => lengths=0 => filtered out)")

    # Output behavior
    p.add_argument("--keep-empty", action="store_true",
                   help="Write output TSV with header even if 0 hits pass filters")

    args = p.parse_args()

    # Validate dirs
    if not args.blast_dir.is_dir():
        sys.exit(f"[ERROR] blast-dir not found: {args.blast_dir}")
    if not args.fungi_dir.is_dir():
        sys.exit(f"[ERROR] fungi-dir not found: {args.fungi_dir}")
    if not args.plant_dir.is_dir():
        sys.exit(f"[ERROR] plant-dir not found: {args.plant_dir}")

    exts = [x.strip() for x in args.extensions.split(",") if x.strip()]
    columns = [c.strip() for c in args.columns.split(",") if c.strip()]

    required_cols = {"qseqid", "sseqid", "pident", "length"}
    if not required_cols.issubset(set(columns)):
        sys.exit(f"[ERROR] --columns must include: {', '.join(sorted(required_cols))}")

    filtered_frames = []
    parsed_files = 0
    used_files = 0

    for blast_file in sorted(args.blast_dir.iterdir()):
        if not blast_file.is_file():
            continue

        parsed = parse_subject_query(blast_file.name, args.blast_suffix)
        if parsed is None:
            continue

        subject_name, query_name = parsed
        parsed_files += 1

        # Locate FASTA files
        fungi_fasta = find_fasta(args.fungi_dir, query_name, exts)
        plant_fasta = find_fasta(args.plant_dir, subject_name, exts)

        if fungi_fasta is None or plant_fasta is None:
            msg = f"[WARNING] Missing FASTA for "
            if fungi_fasta is None:
                msg += f"query={query_name} in {args.fungi_dir} "
            if plant_fasta is None:
                msg += f"subject={subject_name} in {args.plant_dir} "
            eprint(msg.strip())
            if args.require_fai:
                sys.exit("[ERROR] require-fai set and FASTA missing (cannot build/read .fai).")
            continue

        # Ensure / load FAI lengths
        fungi_fai = Path(str(fungi_fasta) + ".fai")
        plant_fai = Path(str(plant_fasta) + ".fai")

        if args.build_fai:
            # Create if missing
            if not fungi_fai.exists():
                ensure_fai(fungi_fasta, args.samtools)
            if not plant_fai.exists():
                ensure_fai(plant_fasta, args.samtools)

        if args.require_fai:
            if not fungi_fai.exists():
                sys.exit(f"[ERROR] Missing FAI for fungi genome: {fungi_fasta}")
            if not plant_fai.exists():
                sys.exit(f"[ERROR] Missing FAI for plant genome: {plant_fasta}")

        fungi_lengths = load_fai(fungi_fai) if fungi_fai.exists() else {}
        plant_lengths = load_fai(plant_fai) if plant_fai.exists() else {}

        # Read BLAST
        try:
            if blast_file.stat().st_size == 0:
                continue
            df = pd.read_csv(blast_file, sep="\t", names=columns, header=None)
        except pd.errors.EmptyDataError:
            continue
        except Exception as e:
            eprint(f"[ERROR] Reading {blast_file.name}: {e}")
            continue

        if df.empty:
            continue

        used_files += 1

        # Provenance
        df["fungi_genome"] = query_name
        df["plant_genome"] = subject_name

        # Filter identity + alignment length
        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df["length"] = pd.to_numeric(df["length"], errors="coerce")

        df = df.dropna(subset=["pident", "length"])
        df = df[(df["pident"] >= args.identity) & (df["length"] >= args.min_align_len)]
        if df.empty:
            continue

        # Filter scaffold/chr length on BOTH sides
        min_scaf = args.min_scaffold_len
        df = df[
            df["qseqid"].map(lambda x: fungi_lengths.get(str(x), 0) >= min_scaf) &
            df["sseqid"].map(lambda x: plant_lengths.get(str(x), 0) >= min_scaf)
        ]

        if not df.empty:
            filtered_frames.append(df)

    # Write output
    if filtered_frames:
        out_df = pd.concat(filtered_frames, ignore_index=True)
        out_df.to_csv(args.out, sep="\t", index=False)
        print(f"[SUCCESS] Saved: {args.out}")
        print(f"[INFO] Parsed BLAST files: {parsed_files}, used: {used_files}")
        print(f"[INFO] Hits kept: {len(out_df)}")
    else:
        if args.keep_empty:
            # write header-only TSV
            header_cols = columns + ["fungi_genome", "plant_genome"]
            pd.DataFrame(columns=header_cols).to_csv(args.out, sep="\t", index=False)
            print(f"[WARNING] No hits passed filters. Wrote empty TSV with header: {args.out}")
        else:
            print("[WARNING] No hits passed filters. Output file not created.")


if __name__ == "__main__":
    main()
