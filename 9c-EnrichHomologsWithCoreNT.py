#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
From homologs/: build queries per multifasta, BLAST core_nt, write enriched results to homologs_nt/
WITH original homologs sequences (INCLUDING the query sequences) + added core_nt sequences.
"""

import argparse
import sys
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from shutil import which


# -------------------------
# Utils
# -------------------------

def eprint(*args):
    print(*args, file=sys.stderr)


def check_tool(name: str):
    if which(name) is None:
        sys.exit(f"[ERROR] Required tool not found in PATH: {name}")


def run_capture(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


def parse_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


def parse_int(x: str) -> Optional[int]:
    try:
        return int(float(x))
    except Exception:
        return None


def sanitize_species_name(s: str) -> str:
    """Make a safe species label for FASTA headers."""
    if not s:
        return "unknown"
    s = s.strip().split(";")[0].strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s or "unknown"


def load_taxid_set(path: Path) -> Set[str]:
    """Load one taxid per line into a set of strings."""
    tax = set()
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            t = t.split()[0]
            tax.add(t)
    return tax


def normalize_id_unique(base_id: str, existing: Set[str]) -> str:
    if base_id not in existing:
        return base_id
    i = 2
    while True:
        cand = f"{base_id}_v{i}"
        if cand not in existing:
            return cand
        i += 1


def extract_taxid_from_header(rec: SeqRecord) -> str:
    """
    Try to extract taxid from input record header without changing anything.
    Supports patterns: taxid=12345 or taxid_12345
    """
    txt = (rec.id or "") + " " + (rec.description or "")
    m = re.search(r"(?:taxid=|taxid_)(\d+)", txt)
    return m.group(1) if m else ""


# -------------------------
# BLAST parsing / selection
# -------------------------

BLAST_FIELDS = [
    "qseqid", "sseqid",
    "pident", "length", "qlen",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
    "staxids", "sscinames"
]


def run_blast_multiquery(
    query_fasta: Path,
    core_db: str,
    blast_program: str,
    threads: int,
    evalue: float,
    extra_args: List[str],
    out_tsv: Path
) -> None:
    outfmt = "6 " + " ".join(BLAST_FIELDS)
    cmd = [
        blast_program,
        "-query", str(query_fasta),
        "-db", core_db,
        "-outfmt", outfmt,
        "-num_threads", str(threads),
        "-evalue", str(evalue),
        "-max_hsps", "1",
    ] + extra_args + [
        "-out", str(out_tsv)
    ]
    rc, _, err = run_capture(cmd)
    if rc != 0:
        eprint(err)
        sys.exit(f"[ERROR] BLAST failed (exit={rc}).\nCommand:\n  {' '.join(cmd)}")


def load_blast_hits(tsv_path: Path) -> List[Dict[str, str]]:
    hits: List[Dict[str, str]] = []
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return hits
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(BLAST_FIELDS):
                continue
            hits.append({BLAST_FIELDS[i]: parts[i] for i in range(len(BLAST_FIELDS))})
    return hits


def kingdom_from_taxid(taxid: str, plants: Set[str], fungi: Set[str]) -> str:
    if not taxid:
        return "other"
    if taxid in fungi:
        return "fungi"
    if taxid in plants:
        return "plant"
    return "other"


def hit_qcov_pct(h: Dict[str, str]) -> Optional[float]:
    length = parse_float(h.get("length", ""))
    qlen = parse_float(h.get("qlen", ""))
    if length is None or qlen is None or qlen == 0:
        return None
    return (length / qlen) * 100.0


def sort_hits_by_quality(hits: List[Dict[str, str]]) -> List[Dict[str, str]]:
    def keyfun(hh):
        ev = parse_float(hh.get("evalue", "1e9"))
        bs = parse_float(hh.get("bitscore", "0"))
        return (ev if ev is not None else 1e9, -(bs if bs is not None else 0))
    return sorted(hits, key=keyfun)


def filter_and_select_hits_by_query_kind(
    hits: List[Dict[str, str]],
    plants_taxids: Set[str],
    fungi_taxids: Set[str],
    q_to_qkind: Dict[str, str],
    min_pident: float,
    min_qcov: float,
    max_hits_per_species: int,
    max_hits_plant: int,
    max_hits_fungi: int,
    exclude_self: bool = True
) -> Dict[str, List[Dict[str, str]]]:
    """
    For each query:
      - if query kind = plant => keep only PLANT hits (cap max_hits_plant)
      - if query kind = fungi => keep only FUNGI hits (cap max_hits_fungi)
    plus max_hits_per_species within the target kingdom.
    """
    by_q: Dict[str, List[Dict[str, str]]] = {}

    # Filter
    for h in hits:
        q = h.get("qseqid", "")
        if not q or q not in q_to_qkind:
            continue

        pident = parse_float(h.get("pident", ""))
        qcov = hit_qcov_pct(h)
        if pident is None or qcov is None:
            continue
        if pident < min_pident:
            continue
        if qcov < min_qcov:
            continue
        if exclude_self and h.get("sseqid", "") == q:
            continue

        by_q.setdefault(q, []).append(h)

    selected: Dict[str, List[Dict[str, str]]] = {}

    for q, qhits in by_q.items():
        qhits = sort_hits_by_quality(qhits)
        target = q_to_qkind[q]  # "plant" or "fungi"
        cap = max_hits_plant if target == "plant" else max_hits_fungi

        out: List[Dict[str, str]] = []
        per_sp: Dict[str, int] = {}

        for h in qhits:
            tax = (h.get("staxids", "") or "").split(";")[0].strip()
            kind = kingdom_from_taxid(tax, plants_taxids, fungi_taxids)
            if kind != target:
                continue

            sci = sanitize_species_name(h.get("sscinames", ""))
            if sci == "unknown":
                sci = f"taxid_{tax}" if tax else "unknown"

            if len(out) >= cap:
                break
            if max_hits_per_species > 0 and per_sp.get(sci, 0) >= max_hits_per_species:
                continue

            out.append(h)
            per_sp[sci] = per_sp.get(sci, 0) + 1

        selected[q] = out

    return selected


# -------------------------
# Fetch sequences (blastdbcmd)
# -------------------------

def fetch_subseq_blastdbcmd(core_db: str, sseqid: str, sstart: int, send: int) -> Optional[str]:
    start = min(sstart, send)
    end = max(sstart, send)
    strand = "plus" if sstart <= send else "minus"
    cmd = ["blastdbcmd", "-db", core_db, "-entry", sseqid, "-range", f"{start}-{end}", "-strand", strand]
    rc, out, err = run_capture(cmd)
    if rc != 0 or not out.strip():
        eprint(f"[WARN] blastdbcmd failed for {sseqid}:{start}-{end} ({strand})")
        if err.strip():
            eprint(err.strip())
        return None
    return out


def make_core_nt_record(
    fasta_text: str,
    query_file_stem: str,
    hit_row: Dict[str, str],
    plants_taxids: Set[str],
    fungi_taxids: Set[str]
) -> Optional[SeqRecord]:
    """
    KEEP HEADER STYLE (unchanged):
      >{query_file_stem}__{fungi|plant|other}__{Species}__{sseqid}:{start}-{end}(+/-)
    """
    from io import StringIO
    recs = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    if not recs:
        return None
    rec = recs[0]

    sseqid = hit_row.get("sseqid", "NA")
    sstart_s = hit_row.get("sstart", "NA")
    send_s = hit_row.get("send", "NA")
    pident = hit_row.get("pident", "")
    evalue = hit_row.get("evalue", "")

    tax = (hit_row.get("staxids", "") or "").split(";")[0].strip()
    kind = kingdom_from_taxid(tax, plants_taxids, fungi_taxids)

    sp = sanitize_species_name(hit_row.get("sscinames", ""))
    if sp == "unknown" and tax:
        sp = f"taxid_{tax}"

    ss = parse_int(sstart_s)
    se = parse_int(send_s)
    if ss is None or se is None:
        loc = f"{sstart_s}-{send_s}(?)"
    else:
        strand_sym = "-" if ss > se else "+"
        start = min(ss, se)
        end = max(ss, se)
        loc = f"{start}-{end}({strand_sym})"

    fasta_id = f"{query_file_stem}__{kind}__{sp}__{sseqid}:{loc}"
    rec.id = fasta_id
    rec.name = fasta_id
    rec.description = f"core_nt_hit taxid={tax} species={sp} pident={pident} evalue={evalue}"
    return rec


# -------------------------
# Query picking by kingdom (input homologs)
# -------------------------

def kingdom_of_input_record(rec: SeqRecord, plants: Set[str], fungi: Set[str]) -> str:
    """
    Detect plant/fungi in INPUT records without changing headers:
    - Prefer explicit markers in id/description: __plant__ or __fungi__
    - Fallback to taxid in header (taxid= / taxid_) using provided sets
    """
    txt = ((rec.id or "") + " " + (rec.description or "")).lower()
    if "__fungi__" in txt:
        return "fungi"
    if "__plant__" in txt:
        return "plant"

    tax = extract_taxid_from_header(rec)
    if tax:
        return kingdom_from_taxid(tax, plants, fungi)

    return "other"


def species_from_record_header(rec: SeqRecord) -> str:
    """
    Extract species/genome label from the __plant__Species__ or __fungi__Species__
    header pattern used throughout the pipeline. Falls back to the record ID.
    """
    parts = rec.id.split("__")
    for i, part in enumerate(parts):
        if part in ("plant", "fungi") and i + 1 < len(parts):
            return parts[i + 1]
    return rec.id


def pick_representatives_by_kind(
    records: List[SeqRecord],
    kind: str,
    plants: Set[str],
    fungi: Set[str],
    max_reps: int,
) -> List[SeqRecord]:
    """
    Select up to max_reps sequences of the given kingdom, one per species
    (longest sequence wins per species), sorted by length descending.
    Replaces pick_longest_record_by_kind to improve BLAST query diversity.
    """
    by_species: Dict[str, SeqRecord] = {}
    for r in records:
        if kingdom_of_input_record(r, plants, fungi) != kind:
            continue
        sp = species_from_record_header(r)
        if sp not in by_species or len(r.seq) > len(by_species[sp].seq):
            by_species[sp] = r
    reps = sorted(by_species.values(), key=lambda r: len(r.seq), reverse=True)
    return reps[:max_reps]


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="From homologs/: build queries (longest plant + longest fungi per multifasta), BLAST core_nt, "
                    "write enriched results to homologs_nt/ WITH original homologs + added core_nt."
    )

    ap.add_argument("--homologs-dir", required=True, type=Path,
                    help="Directory with existing homolog multifastas (e.g. Result_HT/homologs)")
    ap.add_argument("--pattern", default="*.fasta", help="Pattern (default: *.fasta)")

    ap.add_argument("--core-db", required=True, type=str, help="BLAST DB name/path for core_nt (e.g. core_nt)")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="Output directory (e.g. Result_HT/homologs_nt). Files keep same names as in homologs-dir.")

    # Taxonomy labeling for headers
    ap.add_argument("--plants-taxids", required=True, type=Path,
                    help="Taxid list for plant species (one taxid per line)")
    ap.add_argument("--fungi-taxids", required=True, type=Path,
                    help="Taxid list for fungi species (one taxid per line)")

    # Optional restriction for BLAST
    ap.add_argument("--taxidlist", type=Path, default=None,
                    help="Optional taxidlist passed to BLAST -taxidlist (e.g. plants_fungi_species.taxids)")

    # BLAST parameters
    ap.add_argument("--blast-program", default="blastn", help="BLAST program (default: blastn)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-50)
    ap.add_argument("--min-pident", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)

    ap.add_argument("--max-hits-per-species", type=int, default=3,
                    help="Max hits per species WITHIN the target kingdom (default: 3). Use 0 to disable.")
    ap.add_argument("--max-hits-plant", type=int, default=50,
                    help="Max selected PLANT hits per plant-query (default: 50)")
    ap.add_argument("--max-hits-fungi", type=int, default=50,
                    help="Max selected FUNGI hits per fungi-query (default: 50)")
    ap.add_argument("--max-hits-other", type=int, default=50,
                    help="Kept for compatibility (not used by this query strategy).")

    ap.add_argument("--max-queries-per-kingdom", type=int, default=5,
                    help="Max representative sequences per kingdom used as BLAST queries "
                         "(one per species, longest first; default: 5)")

    ap.add_argument("--exclude-self", action="store_true", default=True,
                    help="Exclude hits where sseqid == qseqid (default: ON)")
    ap.add_argument("--include-self", dest="exclude_self", action="store_false")

    ap.add_argument("--summary", type=Path, default=None,
                    help="Summary TSV (default: <outdir>/core_nt_enrichment_summary.tsv)")
    ap.add_argument("--keep-temp", action="store_true",
                    help="Keep temporary query fasta / blast tsv")

    args = ap.parse_args()

    if not args.homologs_dir.is_dir():
        sys.exit(f"[ERROR] homologs-dir not found: {args.homologs_dir}")

    if not args.plants_taxids.exists():
        sys.exit(f"[ERROR] plants-taxids not found: {args.plants_taxids}")
    if not args.fungi_taxids.exists():
        sys.exit(f"[ERROR] fungi-taxids not found: {args.fungi_taxids}")
    if args.taxidlist is not None and not args.taxidlist.exists():
        sys.exit(f"[ERROR] taxidlist not found: {args.taxidlist}")

    check_tool(args.blast_program)
    check_tool("blastdbcmd")

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary or (args.outdir / "core_nt_enrichment_summary.tsv")

    plants_set = load_taxid_set(args.plants_taxids)
    fungi_set = load_taxid_set(args.fungi_taxids)

    fasta_files = sorted(args.homologs_dir.glob(args.pattern))
    if not fasta_files:
        sys.exit(f"[ERROR] No FASTA files found in {args.homologs_dir} with pattern {args.pattern}")

    # Temp files
    tmp_query_fasta = args.outdir / "_tmp_representative_queries.fasta"
    tmp_blast_tsv = args.outdir / "_tmp_core_nt_blast.tsv"

    # Query mappings: qseqid -> source file / kingdom
    q_to_file: Dict[str, Path] = {}
    q_to_qkind: Dict[str, str] = {}

    # All query IDs grouped by stem, for hit aggregation
    stem_to_qids: Dict[str, List[str]] = {}

    # stem -> source file, for output writing
    stem_to_file: Dict[str, Path] = {}

    queries: List[SeqRecord] = []

    for fp in fasta_files:
        stem = fp.stem
        stem_to_file[stem] = fp
        stem_to_qids[stem] = []

        recs = list(SeqIO.parse(str(fp), "fasta"))
        if not recs:
            eprint(f"[WARN] Empty fasta, skipping: {fp.name}")
            continue

        for kind in ("plant", "fungi"):
            reps = pick_representatives_by_kind(
                recs, kind, plants_set, fungi_set, args.max_queries_per_kingdom
            )
            for i, rep in enumerate(reps):
                qid = f"{stem}__Q{kind}_{i}"
                qrec = rep[:]
                qrec.id = qid
                qrec.name = qid
                qrec.description = ""
                queries.append(qrec)
                q_to_file[qid] = fp
                q_to_qkind[qid] = kind
                stem_to_qids[stem].append(qid)

    if not queries:
        sys.exit("[ERROR] No queries extracted. (No plant/fungi records detected in inputs?)")

    SeqIO.write(queries, str(tmp_query_fasta), "fasta")
    print(f"[INFO] Queries built (up to {args.max_queries_per_kingdom} representatives per kingdom "
          f"per file, one per species): {len(queries)}")
    print(f"[INFO] Temp query FASTA: {tmp_query_fasta}")

    extra_args: List[str] = []
    if args.taxidlist is not None:
        extra_args += ["-taxidlist", str(args.taxidlist)]

    # Run BLAST
    print("[INFO] Running BLAST against core_nt ...")
    run_blast_multiquery(
        query_fasta=tmp_query_fasta,
        core_db=args.core_db,
        blast_program=args.blast_program,
        threads=args.threads,
        evalue=args.evalue,
        extra_args=extra_args,
        out_tsv=tmp_blast_tsv
    )

    hits = load_blast_hits(tmp_blast_tsv)
    print(f"[INFO] Loaded raw hits: {len(hits)}")

    # Select hits by query kind (plant-query => plant hits; fungi-query => fungi hits)
    selected_by_q = filter_and_select_hits_by_query_kind(
        hits=hits,
        plants_taxids=plants_set,
        fungi_taxids=fungi_set,
        q_to_qkind=q_to_qkind,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        max_hits_per_species=args.max_hits_per_species,
        max_hits_plant=args.max_hits_plant,
        max_hits_fungi=args.max_hits_fungi,
        exclude_self=args.exclude_self
    )

    # Write per-candidate outputs (keep same filename as homologs/)
    with summary_path.open("w", encoding="utf-8") as s:
        s.write(
            "candidate_file\tqseqid\tstatus\tn_in_homologs\t"
            "n_selected_total\tn_selected_plant\tn_selected_fungi\tn_selected_other\t"
            "n_added_core_nt\treason\n"
        )

        for stem in sorted(stem_to_file.keys()):
            src_fp = stem_to_file[stem]
            out_fp = args.outdir / src_fp.name  # keep same filename

            recs = list(SeqIO.parse(str(src_fp), "fasta"))
            if not recs:
                s.write(f"{src_fp.name}\t{stem}\tEMPTY\t0\t0\t0\t0\t0\t0\tempty_source_fasta\n")
                continue

            # KEEP ALL original records (do NOT remove the longest query)
            homolog_recs = recs[:]

            # Dedup sets against existing homologs
            existing_ids = {r.id for r in homolog_recs}
            existing_seqs = {str(r.seq) for r in homolog_recs}

            # Collect selected hits from all representative queries for this file
            sel_hits: List[Dict[str, str]] = []
            for qid in stem_to_qids.get(stem, []):
                sel_hits.extend(selected_by_q.get(qid, []))

            # Count selected by kingdom (should be plant+fungi only with this strategy)
            sel_pl = sel_fu = sel_ot = 0
            for h in sel_hits:
                tax = (h.get("staxids", "") or "").split(";")[0].strip()
                k = kingdom_from_taxid(tax, plants_set, fungi_set)
                if k == "plant":
                    sel_pl += 1
                elif k == "fungi":
                    sel_fu += 1
                else:
                    sel_ot += 1

            added: List[SeqRecord] = []
            fetch_fail = 0

            for h in sel_hits:
                sseqid = h.get("sseqid", "")
                ss = parse_int(h.get("sstart", "0"))
                se = parse_int(h.get("send", "0"))
                if not sseqid or ss is None or se is None:
                    fetch_fail += 1
                    continue

                fasta_text = fetch_subseq_blastdbcmd(args.core_db, sseqid, ss, se)
                if fasta_text is None:
                    fetch_fail += 1
                    continue

                # IMPORTANT: keep headers unchanged => use original file stem
                rec = make_core_nt_record(
                    fasta_text=fasta_text,
                    query_file_stem=stem,
                    hit_row=h,
                    plants_taxids=plants_set,
                    fungi_taxids=fungi_set
                )
                if rec is None:
                    fetch_fail += 1
                    continue

                # Avoid duplicates
                if rec.id in existing_ids:
                    continue
                if str(rec.seq) in existing_seqs:
                    continue

                rec.id = normalize_id_unique(rec.id, existing_ids)
                rec.name = rec.id

                existing_ids.add(rec.id)
                existing_seqs.add(str(rec.seq))
                added.append(rec)

            # Write final output: ALL original homologs + added core_nt
            SeqIO.write(homolog_recs + added, str(out_fp), "fasta")

            status = "OK"
            reason = "-"
            if not sel_hits:
                status = "NO_HITS"
                reason = "no_core_nt_hits_after_filtering"
            elif len(added) == 0:
                status = "NO_NEW_SEQS"
                reason = "all_hits_failed_or_duplicates"
            elif fetch_fail > 0:
                status = "OK_WITH_WARNINGS"
                reason = f"fetch_failures={fetch_fail}"

            s.write(
                f"{src_fp.name}\t{stem}\t{status}\t{len(recs)}\t"
                f"{len(sel_hits)}\t{sel_pl}\t{sel_fu}\t{sel_ot}\t"
                f"{len(added)}\t{reason}\n"
            )

    print("[SUCCESS] Done.")
    print(f"[INFO] Output directory: {args.outdir}")
    print(f"[INFO] Summary TSV     : {summary_path}")

    if not args.keep_temp:
        try:
            tmp_query_fasta.unlink(missing_ok=True)
            tmp_blast_tsv.unlink(missing_ok=True)
        except Exception:
            pass


if __name__ == "__main__":
    main()
