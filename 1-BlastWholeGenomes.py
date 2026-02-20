#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

DEFAULT_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore"
)

def run(cmd):
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        raise SystemExit(f"Error: command not found: {cmd[0]}")
    except subprocess.CalledProcessError:
        raise SystemExit(f"Error running: {' '.join(cmd)}")

def blast_db_exists(db_prefix: Path) -> bool:
    """
    With -out <db_prefix>, makeblastdb typically produces:
      <db_prefix>.nin/.nsq/.nhr (+ .ndb/.not/.ntf/.nto depending on version)
    We consider the DB existing if the core files exist.
    """
    core_exts = (".nin", ".nsq", ".nhr")
    return all(db_prefix.with_name(db_prefix.name + ext).exists() for ext in core_exts)

def build_blast_db(fasta: Path, db_prefix: Path):
    if blast_db_exists(db_prefix):
        return

    print(f"Building BLAST database for {fasta.name}")
    run([
        "makeblastdb",
        "-in", str(fasta),
        "-dbtype", "nucl",
        "-out", str(db_prefix)
    ])

def blast_done(final_out: Path) -> bool:
    # "Done" means file exists and is non-empty
    return final_out.exists() and final_out.stat().st_size > 0

def filter_by_alignment_length(raw_file: Path, final_file: Path, min_len: int):
    with raw_file.open() as fin, final_file.open("w") as fout:
        for line in fin:
            cols = line.rstrip().split("\t")
            if len(cols) >= 4 and int(cols[3]) >= min_len:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Whole genome BLAST (fungi vs plants) with BLAST+"
    )

    # ===== Input / Output options =====
    parser.add_argument(
        "--plants-dir",
        required=True,
        type=Path,
        help="Directory containing plant genomes (.fasta)"
    )
    parser.add_argument(
        "--fungi-dir",
        required=True,
        type=Path,
        help="Directory containing fungi genomes (.fasta)"
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output directory name (e.g. Result_HT)"
    )

    # ===== BLAST options =====
    parser.add_argument("--evalue", default="1e-50", help="E-value threshold")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--identity", type=float, default=70.0, help="Min percent identity")
    parser.add_argument("--min-align-len", type=int, default=500, help="Min alignment length (bp)")
    parser.add_argument("--outfmt", default=DEFAULT_OUTFMT, help="BLAST outfmt")

    # ===== Force rerun =====
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rerun (rebuild DB and rerun BLAST) even if outputs exist"
    )

    args = parser.parse_args()

    # ===== Directories =====
    blast_dir = args.output / "blast"
    db_dir = args.output / "blastdb"
    tmp_dir = args.output / "tmp"

    for d in (blast_dir, db_dir, tmp_dir):
        d.mkdir(parents=True, exist_ok=True)

    plant_fastas = sorted(args.plants_dir.glob("*.fasta"))
    fungi_fastas = sorted(args.fungi_dir.glob("*.fasta"))

    if not plant_fastas:
        raise SystemExit("No plant FASTA found")
    if not fungi_fastas:
        raise SystemExit("No fungi FASTA found")

    # ===== Main loop =====
    for plant in plant_fastas:
        # Keep the existing naming style so it matches your current DB files:
        # e.g. Arabidopsis.thaliana.fasta.nin
        db_prefix = db_dir / plant.name

        if args.force:
            # If forcing, remove existing DB core files so makeblastdb rebuilds cleanly
            for ext in (".nin", ".nsq", ".nhr", ".ndb", ".not", ".ntf", ".nto"):
                p = db_prefix.with_name(db_prefix.name + ext)
                if p.exists():
                    p.unlink()
        build_blast_db(plant, db_prefix)

        for fungi in fungi_fastas:
            final_out = blast_dir / f"{plant.name}_VS_{fungi.name}.blast"
            raw_out = tmp_dir / f"{plant.name}_VS_{fungi.name}.raw.blast"

            # Skip if already done (unless --force)
            if blast_done(final_out) and not args.force:
                print(f"SKIP (already done): {fungi.name} vs {plant.name}")
                continue

            print(f"BLAST {fungi.name} vs {plant.name}")

            run([
                "blastn",
                "-query", str(fungi),
                "-db", str(db_prefix),
                "-evalue", str(args.evalue),
                "-num_threads", str(args.threads),
                "-perc_identity", str(args.identity),
                "-outfmt", str(args.outfmt),
                "-out", str(raw_out)
            ])

            filter_by_alignment_length(raw_out, final_out, args.min_align_len)
            raw_out.unlink(missing_ok=True)

if __name__ == "__main__":
    main()
