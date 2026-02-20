#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import sys
from typing import Iterator, Tuple, List, Set


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, seq) from a fasta file. header excludes leading '>'."""
    header = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def write_fasta(records: List[Tuple[str, str]], out_path: Path, wrap: int = 60) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i + wrap] + "\n")


def classify_header(header: str) -> str:
    """Return 'fungi', 'plant', or 'unknown' using the pipeline tokens."""
    if "__fungi__" in header:
        return "fungi"
    if "__plant__" in header:
        return "plant"
    return "unknown"


def main():
    ap = argparse.ArgumentParser(
        description="Clean homolog FASTA: remove contaminant species and keep only files with >=1 fungi and >=1 plant."
    )

    ap.add_argument("--in-dir", required=True, type=Path,
                    help="Input directory containing *.fasta")

    # 👉 ajout ici (seul changement)
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="Output directory")

    ap.add_argument("--remove-species", type=str, default="",
                    help="Comma-separated substrings to remove if found in header "
                         "(e.g. 'Betula.nana,Pseudotsuga.menziesii,Quercus.suber')")

    ap.add_argument("--report", type=Path, default=Path("homologs_cleaned.report.tsv"),
                    help="TSV report file (default: homologs_cleaned.report.tsv)")

    args = ap.parse_args()

    in_dir = args.in_dir.resolve()
    out_dir = args.out_dir.resolve()
    report_path = args.report.resolve()

    if not in_dir.is_dir():
        print(f"[ERROR] Input directory not found: {in_dir}", file=sys.stderr)
        sys.exit(1)

    fasta_files = sorted(in_dir.glob("*.fasta"))
    if not fasta_files:
        print(f"[ERROR] No *.fasta found in {in_dir}", file=sys.stderr)
        sys.exit(1)

    # parse species list
    remove_list = []
    if args.remove_species.strip():
        remove_list = [s.strip() for s in args.remove_species.split(",") if s.strip()]

    out_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with report_path.open("w", encoding="utf-8") as rep:
        rep.write("file\tkept_total\tfungi\tplant\tstatus\n")

        total_files = 0
        written_files = 0
        skipped_files = 0
        removed_contam = 0
        removed_dups = 0

        for fp in fasta_files:
            total_files += 1

            kept: List[Tuple[str, str]] = []
            seen_pairs: Set[Tuple[str, str]] = set()

            fungi = 0
            plant = 0

            for header, seq in parse_fasta(fp):

                # remove unwanted species
                if remove_list and any(sp in header for sp in remove_list):
                    removed_contam += 1
                    continue

                # remove only exact duplicates (id + sequence)
                fid = header.split()[0] if header else ""
                key = (fid, seq)
                if key in seen_pairs:
                    removed_dups += 1
                    continue
                seen_pairs.add(key)

                kept.append((header, seq))

                # count types
                k = classify_header(header)
                if k == "fungi":
                    fungi += 1
                elif k == "plant":
                    plant += 1

            # keep only if both present
            if fungi > 0 and plant > 0:
                status = "OK"
                out_path = out_dir / fp.name
                write_fasta(kept, out_path)
                written_files += 1
            else:
                if fungi == 0 and plant == 0:
                    status = "FAIL_NO_BOTH"
                elif fungi == 0:
                    status = "FAIL_NO_FUNGI"
                else:
                    status = "FAIL_NO_PLANT"
                skipped_files += 1

            rep.write(f"{fp.name}\t{len(kept)}\t{fungi}\t{plant}\t{status}\n")

    print("[DONE]")
    print(f"Input dir : {in_dir}")
    print(f"Output dir: {out_dir}")
    print(f"Files processed: {total_files}")
    print(f"Files kept (OK): {written_files}")
    print(f"Files removed (fail): {skipped_files}")
    print(f"Removed contaminants: {removed_contam}")
    print(f"Removed exact duplicates: {removed_dups}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
