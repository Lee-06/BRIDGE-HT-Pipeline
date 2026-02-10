#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 10: Build phylogenies from homolog multi-FASTAs (Step 9 output).

- Input : directory with one multifasta per candidate
- Output: ONE directory with all alignments + trees (no subfolders)
- MAFFT -> TrimAl (-keepheader) -> IQ-TREE2
"""

import argparse
import subprocess
from pathlib import Path
from shutil import which
from Bio import SeqIO
import sys


# -------------------------
# Utils
# -------------------------

def check_tool(name):
    if which(name) is None:
        sys.exit(f"[ERROR] Required tool not found in PATH: {name}")


def run(cmd, quiet=False):
    if not quiet:
        print("       " + " ".join(cmd))
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL if quiet else None,
        stderr=subprocess.DEVNULL if quiet else None,
        check=False
    )


# -------------------------
# MAIN
# -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Step 10: Build phylogenies from homolog multi-FASTAs (no subfolders, TrimAl keepheader)."
    )
    ap.add_argument("--homologs-dir", required=True, type=Path,
                    help="Directory with per-candidate homolog multifastas (Step 9)")
    ap.add_argument("--outdir", default="phylogenie", type=Path,
                    help="Output directory for ALL phylogenies")
    ap.add_argument("--pattern", default="*.fasta",
                    help="Pattern for homolog multifastas (default: *.fasta)")

    ap.add_argument("--min-seqs", type=int, default=4,
                    help="Minimum sequences required to build a tree (default: 4)")
    ap.add_argument("--drop-query", action="store_true",
                    help="Remove the candidate sequence before alignment")

    ap.add_argument("--mafft-threads", type=int, default=8)
    ap.add_argument("--iqtree-threads", type=int, default=8)
    ap.add_argument("--model", default=None,
                    help="IQ-TREE2 model (e.g. MFP, GTR+G)")
    ap.add_argument("--bb", type=int, default=1000,
                    help="UFBoot replicates (default: 1000)")
    ap.add_argument("--resume", action="store_true",
                    help="Skip if treefile already exists")
    ap.add_argument("--quiet", action="store_true")

    args = ap.parse_args()

    # Checks
    if not args.homologs_dir.is_dir():
        sys.exit(f"[ERROR] homologs-dir not found: {args.homologs_dir}")

    check_tool("mafft")
    check_tool("trimal")
    check_tool("iqtree2")

    args.outdir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(args.homologs_dir.glob(args.pattern))
    if not fasta_files:
        sys.exit("[ERROR] No homolog multifastas found")

    print(f"[INFO] Found {len(fasta_files)} candidates")

    built = skipped = failed = 0

    for mf in fasta_files:
        cid = mf.stem
        print(f"\n[INFO] Processing {cid}")

        aln = args.outdir / f"{cid}.aln"
        trimmed = args.outdir / f"{cid}.trimmed.aln"
        treefile = Path(str(trimmed) + ".treefile")

        if args.resume and treefile.exists():
            print("    [SKIP] Tree already exists")
            skipped += 1
            continue

        records = list(SeqIO.parse(mf, "fasta"))
        if args.drop_query:
            records = [r for r in records if r.id != cid]

        if len(records) < args.min_seqs:
            print(f"    [SKIP] Not enough sequences ({len(records)})")
            skipped += 1
            continue

        tmp_fasta = args.outdir / f"{cid}.input.fasta"
        SeqIO.write(records, tmp_fasta, "fasta")

        # MAFFT
        print("    [STEP] MAFFT")
        with open(aln, "w") as out:
            subprocess.run(
                ["mafft", "--auto", "--thread", str(args.mafft_threads), str(tmp_fasta)],
                stdout=out,
                stderr=subprocess.DEVNULL if args.quiet else None
            )

        if not aln.exists() or aln.stat().st_size == 0:
            print("    [FAIL] MAFFT failed")
            failed += 1
            continue

        # TrimAl (KEEP HEADERS!)
        print("    [STEP] TrimAl (-keepheader)")
        run([
            "trimal",
            "-in", str(aln),
            "-out", str(trimmed),
            "-automated1",
            "-keepheader"
        ], quiet=args.quiet)

        if not trimmed.exists() or trimmed.stat().st_size == 0:
            print("    [FAIL] TrimAl failed")
            failed += 1
            continue

        # IQ-TREE2
        print("    [STEP] IQ-TREE2")
        iq_cmd = [
            "iqtree2",
            "-s", str(trimmed),
            "-nt", str(args.iqtree_threads),
            "-bb", str(args.bb)
        ]
        if args.model:
            iq_cmd += ["-m", args.model]

        run(iq_cmd, quiet=args.quiet)

        if not treefile.exists():
            print("    [FAIL] IQ-TREE2 failed")
            failed += 1
            continue

        print(f"    [SUCCESS] Tree built: {treefile.name}")
        built += 1

    print("\n[SUMMARY]")
    print(f"  Built   : {built}")
    print(f"  Skipped : {skipped}")
    print(f"  Failed  : {failed}")


if __name__ == "__main__":
    main()
