#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd, cwd=None):
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Command failed ({e.returncode}): {' '.join(cmd)}")


def main():
    p = argparse.ArgumentParser(
        description="Step 7: Annotate clustered HT candidates (output of Step 6) with EggNOG-mapper."
    )

    # Step 6 output -> Step 7 input
    p.add_argument(
        "--clusters-fasta", required=True, type=Path,
        help="Clustered FASTA from Step 6 (e.g., Result_HT/ht_clusters.fasta)"
    )

    # EggNOG config
    p.add_argument(
        "--eggnog-data-dir", required=True, type=Path,
        help="EggNOG database directory (must contain eggnog_proteins.dmnd)"
    )
    p.add_argument(
        "--dmnd-db", default=None, type=Path,
        help="Optional path to eggnog_proteins.dmnd (default: <eggnog-data-dir>/eggnog_proteins.dmnd)"
    )

    # Output
    p.add_argument(
        "--outdir", default=".", type=Path,
        help="Output directory (default: current directory)"
    )
    p.add_argument(
        "--output-prefix", default="ht_annotations",
        help="EggNOG output prefix (default: ht_annotations)"
    )

    # emapper options
    p.add_argument("--cpu", type=int, default=8, help="CPUs for emapper (default: 8)")
    p.add_argument("--method", default="diamond", choices=["diamond", "hmmer"],
                   help="Search method (-m) (default: diamond)")
    p.add_argument("--itype", default="genome", help="Input type (--itype) (default: genome)")
    p.add_argument("--translate", action="store_true",
                   help="Translate nucleotide sequences (--translate). Recommended for genome DNA.")
    p.add_argument("--emapper", default="emapper.py", help="Path to emapper.py (default: emapper.py)")

    args = p.parse_args()

    # Checks
    if not args.clusters_fasta.exists():
        sys.exit(f"[ERROR] clusters-fasta not found: {args.clusters_fasta}")

    if not args.eggnog_data_dir.is_dir():
        sys.exit(f"[ERROR] eggnog-data-dir not found: {args.eggnog_data_dir}")

    dmnd_db = args.dmnd_db if args.dmnd_db else (args.eggnog_data_dir / "eggnog_proteins.dmnd")
    if not dmnd_db.exists():
        sys.exit(f"[ERROR] eggnog_proteins.dmnd not found: {dmnd_db}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # IMPORTANT: resolve paths (absolute) because we run emapper with cwd=outdir
    clusters_fasta_abs = args.clusters_fasta.resolve()
    eggnog_dir_abs = args.eggnog_data_dir.resolve()
    dmnd_db_abs = dmnd_db.resolve()

    # Run emapper
    cmd = [
        args.emapper,
        "-i", str(clusters_fasta_abs),
        "--output", args.output_prefix,
        "--cpu", str(args.cpu),
        "--dmnd_db", str(dmnd_db_abs),
        "-m", args.method,
        "--data_dir", str(eggnog_dir_abs),
        "--itype", args.itype,
    ]
    if args.translate:
        cmd.append("--translate")

    print("[INFO] Running EggNOG-mapper (Step 7) ...")
    print("       " + " ".join(cmd))
    run(cmd, cwd=str(args.outdir))

    print("[SUCCESS] Step 7 annotation complete.")
    print(f"[INFO] Output prefix: {args.output_prefix}")
    print(f"[INFO] Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
