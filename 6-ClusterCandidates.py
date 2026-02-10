#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd):
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] Command not found: {cmd[0]} (is CD-HIT installed and in PATH?)")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] Command failed ({e.returncode}): {' '.join(cmd)}")


def fasta_sort_by_length_desc(in_fasta: Path, out_fasta: Path) -> None:
    """
    Sort FASTA records by sequence length (descending).
    This helps CD-HIT keep the longest sequence as cluster representative.
    """
    records = []
    cur_header = None
    cur_seq = []

    def flush():
        nonlocal cur_header, cur_seq
        if cur_header is None:
            return
        seq = "".join(cur_seq).replace(" ", "").replace("\r", "")
        records.append((cur_header, seq))
        cur_header = None
        cur_seq = []

    with in_fasta.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                cur_header = line
                cur_seq = []
            else:
                cur_seq.append(line)

    flush()

    records.sort(key=lambda x: len(x[1]), reverse=True)

    with out_fasta.open("w", encoding="utf-8") as out:
        for h, s in records:
            out.write(h + "\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")


def main():
    p = argparse.ArgumentParser(
        description="Cluster filtered HT candidates with cd-hit-est to remove redundancy (supports keeping longest representatives)."
    )

    p.add_argument("--input", "-i", required=True, type=Path,
                   help="Input FASTA (filtered candidates, e.g. ht_candidates.filtered.fasta)")
    p.add_argument("--output", "-o", default="ht_clusters.fasta", type=Path,
                   help="Output clustered FASTA (default: ht_clusters.fasta)")

    # cd-hit-est binary
    p.add_argument("--cdhit", default="cd-hit-est",
                   help="Path to cd-hit-est (default: cd-hit-est)")

    # Parameters (DEFAULTS updated for 'keep longest / containment-like' behavior)
    p.add_argument("--c", type=float, default=0.8,
                   help="Sequence identity threshold (-c). NOTE: cd-hit-est requires >= 0.8 (default: 0.8)")
    p.add_argument("--d", type=int, default=0,
                   help="Description length in .clstr (-d) (default: 0)")
    p.add_argument("--T", type=int, default=0,
                   help="Threads (-T). 0 = all CPUs (default: 0)")
    p.add_argument("--M", type=int, default=0,
                   help="Memory limit in MB (-M). 0 = unlimited (default: 0)")

    # Key options for "small included in large"
    p.add_argument("--G", type=int, default=0,
                   help="Use global (1) or local (0) identity (-G) (default: 0)")
    p.add_argument("--aS", type=float, default=0.4,
                   help="Alignment coverage for shorter sequence (-aS) (default: 0.4)")
    p.add_argument("--aL", type=float, default=0.1,
                   help="Alignment coverage for longer sequence (-aL) (default: 0.1)")

    # Sorting behavior (important!)
    p.add_argument("--no-sort", action="store_true",
                   help="Disable sorting by sequence length (NOT recommended if you want to keep longest representatives).")
    p.add_argument("--keep-sorted", action="store_true",
                   help="Keep the temporary length-sorted FASTA (default: delete it after clustering).")

    # Optional: output clstr file path
    p.add_argument("--clstr", default=None, type=Path,
                   help="Optional: copy the .clstr file to this path")

    args = p.parse_args()

    if not args.input.exists():
        sys.exit(f"[ERROR] Input FASTA not found: {args.input}")

    if args.c < 0.8:
        sys.exit("[ERROR] cd-hit-est does not accept -c < 0.8. Use >= 0.8 (or use mmseqs2/vsearch for lower identity).")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.clstr:
        args.clstr.parent.mkdir(parents=True, exist_ok=True)

    # Sort input by length (desc) unless disabled
    temp_sorted = None
    if args.no_sort:
        input_for_cdhit = args.input
        print("[INFO] Sorting disabled (--no-sort). Representatives may NOT be the longest if input is not ordered.")
    else:
        temp_sorted = args.output.parent / (args.input.stem + ".sorted_by_len.fasta")
        print(f"[INFO] Sorting input FASTA by length (desc): {args.input} -> {temp_sorted}")
        fasta_sort_by_length_desc(args.input, temp_sorted)
        input_for_cdhit = temp_sorted

    cmd = [
        args.cdhit,
        "-i", str(input_for_cdhit),
        "-o", str(args.output),
        "-c", str(args.c),
        "-d", str(args.d),
        "-T", str(args.T),
        "-M", str(args.M),
        "-G", str(args.G),
        "-aS", str(args.aS),
        "-aL", str(args.aL),
    ]

    print("[INFO] Running cd-hit-est:")
    print("       " + " ".join(cmd))
    run(cmd)

    # cd-hit-est writes a .clstr file next to output
    clstr_src = Path(str(args.output) + ".clstr")
    if args.clstr and clstr_src.exists():
        args.clstr.write_text(
            clstr_src.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8"
        )
        print(f"[INFO] Copied cluster file to: {args.clstr}")

    print(f"[SUCCESS] Clustered FASTA saved to: {args.output}")
    if clstr_src.exists():
        print(f"[INFO] Cluster report: {clstr_src}")

    # Cleanup sorted temp file
    if temp_sorted and temp_sorted.exists() and (not args.keep_sorted):
        try:
            temp_sorted.unlink()
            print(f"[INFO] Removed temporary sorted FASTA: {temp_sorted}")
        except Exception:
            print(f"[WARNING] Could not remove temporary file: {temp_sorted}", file=sys.stderr)
    elif temp_sorted and temp_sorted.exists():
        print(f"[INFO] Kept temporary sorted FASTA: {temp_sorted} (--keep-sorted)")


if __name__ == "__main__":
    main()
