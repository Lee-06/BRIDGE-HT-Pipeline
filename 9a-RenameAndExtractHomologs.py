#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

from Bio import SeqIO


# ----------------------------
# Utils: safe output (no overwrite)
# ----------------------------
def unique_path(path: Path) -> Path:
    """If file exists, return path_v2, path_v3, ..."""
    if not path.exists():
        return path
    parent = path.parent
    stem = path.stem
    suffix = path.suffix
    i = 2
    while True:
        p = parent / f"{stem}_v{i}{suffix}"
        if not p.exists():
            return p
        i += 1


# ----------------------------
# String sanitization for IDs
# ----------------------------
_bad = re.compile(r"[^A-Za-z0-9._-]+")
_multi_us = re.compile(r"_+")

def sanitize_token(s: str, maxlen: int = 40) -> str:
    s = (s or "").strip()
    if not s or s == "-" or s.lower() == "nan":
        return "NA"
    s = s.replace("/", "_").replace("\\", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace("[", "").replace("]", "")
    s = s.replace("{", "").replace("}", "")
    s = s.replace(";", "_").replace(":", "_")
    s = s.replace(",", "_")
    s = s.replace("|", "_")
    s = s.replace(" ", "_")
    s = _bad.sub("_", s)
    s = _multi_us.sub("_", s).strip("_")
    if not s:
        return "NA"
    if len(s) > maxlen:
        s = s[:maxlen].rstrip("_")
    return s


def normalize_eggnog_query(q: str) -> str:
    """
    EggNOG query IDs often end with _0, _1 ... when --translate is used.
    Remove ONLY the last _<digits> suffix if present.
    """
    if "_" not in q:
        return q
    left, right = q.rsplit("_", 1)
    if right.isdigit():
        return left
    return q


# ----------------------------
# EggNOG parsing
# ----------------------------
def parse_eggnog_annotations(annotations_path: Path) -> Dict[str, Dict[str, str]]:
    """
    Returns dict keyed by normalized query id, with:
      Description, Preferred_name, PFAMs
    Keeps the first hit per query id.
    """
    header_line = None
    with annotations_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#query\t"):
                header_line = line.strip()[1:]  # remove leading "#"
                break
    if header_line is None:
        raise SystemExit("[ERROR] Cannot find '#query' header line in EggNOG annotations.")

    cols = header_line.split("\t")

    data: Dict[str, Dict[str, str]] = {}
    with annotations_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1:
                continue
            row = {cols[i]: parts[i] if i < len(parts) else "" for i in range(len(cols))}
            qraw = row.get("query", "").strip()
            if not qraw:
                continue
            qid = normalize_eggnog_query(qraw)
            if qid in data:
                continue
            data[qid] = {
                "Description": row.get("Description", "").strip(),
                "Preferred_name": row.get("Preferred_name", "").strip(),
                "PFAMs": row.get("PFAMs", "").strip(),
            }
    return data


def pick_gene_function(meta: Dict[str, str]) -> str:
    pref = (meta.get("Preferred_name") or "").strip()
    desc = (meta.get("Description") or "").strip()
    if pref and pref != "-" and pref.lower() != "nan":
        return pref
    if desc and desc != "-" and desc.lower() != "nan":
        return desc
    return "NA"


def pick_pfam_function(meta: Dict[str, str]) -> str:
    pf = (meta.get("PFAMs") or "").strip()
    if not pf or pf == "-" or pf.lower() == "nan":
        return "NA"
    first = pf.split(",")[0].strip()
    return first if first else "NA"


# ----------------------------
# BLAST utilities
# ----------------------------
def ensure_blast_db(fasta: Path) -> None:
    """
    Make blast db if missing (uses fasta as -out prefix).
    Checks .nin/.nsq/.nhr existence.
    """
    prefix = str(fasta)
    if Path(prefix + ".nin").exists() or Path(prefix + ".nsq").exists() or Path(prefix + ".nhr").exists():
        return
    print(f"[INFO] makeblastdb for: {fasta.name}")
    try:
        subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(fasta)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        raise SystemExit("[ERROR] makeblastdb not found (install NCBI BLAST+).")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[ERROR] makeblastdb failed for {fasta} (code {e.returncode})")


def run_blast_to_file(
    query_fasta: Path,
    db_fasta: Path,
    out_tsv: Path,
    evalue: float,
    max_target_seqs: int,
    threads: int,
) -> None:
    """
    Runs BLAST and writes tabular output to out_tsv.
    """
    outfmt = "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore"
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_fasta),
        "-evalue", str(evalue),
        "-outfmt", outfmt,
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-out", str(out_tsv),
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise SystemExit("[ERROR] blastn not found (install NCBI BLAST+).")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        raise SystemExit(f"[ERROR] blastn failed on DB {db_fasta} (code {e.returncode})")


def read_blast_lines(tsv: Path) -> List[str]:
    if not tsv.exists():
        return []
    with tsv.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.strip() for ln in f if ln.strip()]


def extract_hits_from_genome_index(
    blast_lines: List[str],
    genome_index,
    min_identity: float,
    min_coverage: float,
    min_scaffold_length: int,
) -> Dict[str, List[Tuple[str, int, int, str, str]]]:
    """
    hits_by_query[qid] = list of (sid, start, end, strand, seq)
    Filters:
      pident >= min_identity
      alen >= min_coverage * qlen
      scaffold length >= min_scaffold_length
    """
    hits_by_query = defaultdict(list)

    for line in blast_lines:
        parts = line.split("\t")
        if len(parts) != 10:
            continue

        qid, sid, pident, alen, qlen, qstart, qend, sstart, send, bitscore = parts
        try:
            pident = float(pident)
            alen = int(alen)
            qlen = int(qlen)
            s1, s2 = int(sstart), int(send)
        except ValueError:
            continue

        if pident < min_identity:
            continue
        if alen < int(min_coverage * qlen):
            continue

        if sid not in genome_index:
            continue
        try:
            subj_len = len(genome_index[sid].seq)
        except Exception:
            continue
        if subj_len < min_scaffold_length:
            continue

        strand = "+" if s1 <= s2 else "-"
        start, end = (s1, s2) if s1 <= s2 else (s2, s1)

        seq = genome_index[sid].seq[start - 1:end]
        if strand == "-":
            seq = seq.reverse_complement()

        hits_by_query[qid].append((sid, start, end, strand, str(seq)))

    return hits_by_query


# ----------------------------
# MAIN
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Step 9: Rename filtered HT candidates using EggNOG info, then BLAST vs plant+fungi genomes and extract homologs."
    )

    ap.add_argument("--candidates", required=True, type=Path,
                    help="Filtered candidates FASTA (output of Step 8)")
    ap.add_argument("--annotations", required=True, type=Path,
                    help="EggNOG .emapper.annotations (from Step 7)")

    ap.add_argument("--fungi-dir", required=True, type=Path,
                    help="Directory containing fungi genome FASTAs (*.fasta/*.fa/*.fna)")
    ap.add_argument("--plant-dir", required=True, type=Path,
                    help="Directory containing plant genome FASTAs (*.fasta/*.fa/*.fna)")

    ap.add_argument("--outdir", required=True, type=Path,
                    help="Main output directory (e.g., Result_HT)")
    ap.add_argument("--homologs-dir", default="homologs", type=str,
                    help="Subdirectory name inside --outdir for per-candidate homolog multifastas (default: homologs)")

    # NEW: BLAST results saving
    ap.add_argument("--blast-results-dir", default="blast_results", type=str,
                    help="Subdirectory inside --outdir to store BLAST TSV per genome (default: blast_results)")
    ap.add_argument("--no-reuse-blast", action="store_true",
                    help="Do NOT reuse existing BLAST TSV files; rerun BLAST unless --force-blast is also false (default: reuse)")
    ap.add_argument("--force-blast", action="store_true",
                    help="Force rerun BLAST even if TSV exists (overwrites TSV).")

    # Renaming outputs
    ap.add_argument("--renamed-fasta", default="ht_candidates.renamed.fasta", type=str,
                    help="Renamed candidates FASTA filename (in --outdir)")
    ap.add_argument("--id-map", default="ht_candidates.id_map.tsv", type=str,
                    help="ID mapping TSV filename (in --outdir)")

    # BLAST filtering
    ap.add_argument("--identity", type=float, default=80.0, help="Minimum % identity (default: 80)")
    ap.add_argument("--coverage", type=float, default=0.8, help="Minimum query coverage (0-1) (default: 0.8)")
    ap.add_argument("--evalue", type=float, default=1e-20, help="BLAST e-value (default: 1e-20)")
    ap.add_argument("--max-seqs", type=int, default=10, help="blastn -max_target_seqs per query (default: 10)")
    ap.add_argument("--threads", type=int, default=8, help="blastn threads (default: 8)")
    ap.add_argument("--min-scaffold-length", type=int, default=20000,
                    help="Minimum scaffold/chr length in subject genome to keep hit (default: 20000)")

    args = ap.parse_args()

    # Checks
    for p in [args.candidates, args.annotations]:
        if not p.exists():
            raise SystemExit(f"[ERROR] File not found: {p}")
    for d in [args.fungi_dir, args.plant_dir]:
        if not d.is_dir():
            raise SystemExit(f"[ERROR] Directory not found: {d}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    homologs_dir = args.outdir / args.homologs_dir
    homologs_dir.mkdir(parents=True, exist_ok=True)

    blast_dir = args.outdir / args.blast_results_dir
    blast_dir.mkdir(parents=True, exist_ok=True)

    reuse_blast = not args.no_reuse_blast

    # Parse EggNOG metadata
    eggnog = parse_eggnog_annotations(args.annotations)

    # Rename candidates
    renamed_path = unique_path(args.outdir / args.renamed_fasta)
    idmap_path = unique_path(args.outdir / args.id_map)

    print("[INFO] Renaming candidate IDs...")
    old_to_new: Dict[str, str] = {}
    new_to_meta: Dict[str, Dict[str, str]] = {}

    candidate_records = list(SeqIO.parse(str(args.candidates), "fasta"))
    if not candidate_records:
        raise SystemExit("[ERROR] No records found in candidates FASTA.")

    with idmap_path.open("w", encoding="utf-8") as m:
        m.write("old_id\tnew_id\tgene_function\tpfam_first\tdescription\tpreferred_name\tpfams\n")

        for i, rec in enumerate(candidate_records, start=1):
            old_id = rec.id
            meta = eggnog.get(old_id, {"Description": "NA", "Preferred_name": "NA", "PFAMs": "NA"})

            gene_fun_raw = pick_gene_function(meta)
            pfam_raw = pick_pfam_function(meta)

            gene_fun = sanitize_token(gene_fun_raw, maxlen=35)
            pfam_fun = sanitize_token(pfam_raw, maxlen=25)

            new_id = f"HTcandidate_{i:05d}_{gene_fun}_{pfam_fun}"

            # guarantee uniqueness even if function texts collide
            if new_id in new_to_meta:
                suffix = 2
                base = new_id
                while f"{base}_{suffix}" in new_to_meta:
                    suffix += 1
                new_id = f"{base}_{suffix}"

            old_to_new[old_id] = new_id
            new_to_meta[new_id] = meta

            m.write(
                f"{old_id}\t{new_id}\t{gene_fun_raw}\t{pfam_raw}\t"
                f"{meta.get('Description','')}\t{meta.get('Preferred_name','')}\t{meta.get('PFAMs','')}\n"
            )

            rec.id = new_id
            rec.name = new_id
            rec.description = ""

    with renamed_path.open("w", encoding="utf-8") as out:
        SeqIO.write(candidate_records, out, "fasta")

    print(f"[SUCCESS] Renamed FASTA: {renamed_path}")
    print(f"[SUCCESS] ID map TSV  : {idmap_path}")

    # Prepare genomes list
    def list_fastas(d: Path) -> List[Path]:
        out = []
        for ext in (".fasta", ".fa", ".fna"):
            out.extend(sorted(d.glob(f"*{ext}")))
        return out

    fungi_fastas = list_fastas(args.fungi_dir)
    plant_fastas = list_fastas(args.plant_dir)

    if not fungi_fastas and not plant_fastas:
        raise SystemExit("[ERROR] No genome FASTA found in fungi-dir/plant-dir.")

    genomes: List[Tuple[Path, str]] = [(p, "fungi") for p in fungi_fastas] + [(p, "plant") for p in plant_fastas]
    print(f"[INFO] Genomes found: fungi={len(fungi_fastas)} plant={len(plant_fastas)} total={len(genomes)}")

    # BLAST and extract homologs
    print("[INFO] Running BLAST and extracting homologs...")
    for genome_fa, gtype in genomes:
        gname = genome_fa.name
        print(f"[INFO] -> {gtype}: {gname}")

        ensure_blast_db(genome_fa)

        # BLAST output file for this genome
        # include gtype to avoid name collisions if same filename exists in both dirs
        blast_tsv = blast_dir / f"{gtype}__{genome_fa.stem}.blast.tsv"

        need_run = True
        if blast_tsv.exists():
            if args.force_blast:
                need_run = True
            elif reuse_blast:
                need_run = False
            else:
                need_run = True

        if need_run:
            if args.force_blast and blast_tsv.exists():
                # overwrite on purpose
                blast_tsv.unlink()
            print(f"[INFO]    blastn -> {blast_tsv.name}")
            run_blast_to_file(
                query_fasta=renamed_path,
                db_fasta=genome_fa,
                out_tsv=blast_tsv,
                evalue=args.evalue,
                max_target_seqs=args.max_seqs,
                threads=args.threads,
            )
        else:
            print(f"[INFO]    reuse BLAST TSV: {blast_tsv.name}")

        blast_lines = read_blast_lines(blast_tsv)
        if not blast_lines:
            continue

        # Index genome (lazy access; avoids loading full genome into RAM)
        try:
            genome_index = SeqIO.index(str(genome_fa), "fasta")
        except Exception as e:
            print(f"[WARNING] Could not index {genome_fa}: {e}", file=sys.stderr)
            continue

        hits = extract_hits_from_genome_index(
            blast_lines=blast_lines,
            genome_index=genome_index,
            min_identity=args.identity,
            min_coverage=args.coverage,
            min_scaffold_length=args.min_scaffold_length,
        )

        # write homolog multifastas per candidate (append)
        for qid, regions in hits.items():
            out_fa = homologs_dir / f"{qid}.fasta"
            with out_fa.open("a", encoding="utf-8") as out:
                for sid, start, end, strand, seq in regions:
                    header = f"{qid}__{gtype}__{gname}__{sid}:{start}-{end}({strand})"
                    out.write(f">{header}\n{seq}\n")

        # close index (important on some systems)
        try:
            genome_index.close()
        except Exception:
            pass

    print(f"[DONE] BLAST TSV written to: {blast_dir}")
    print(f"[DONE] Homolog multifastas written to: {homologs_dir}")


if __name__ == "__main__":
    main()

