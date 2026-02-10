#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
9b-FetchHomologs.py

- Input : dossier homologs/ avec 1 multifasta par candidat
- Query : la séquence la plus longue de chaque multifasta (1 query par fichier)
- BLAST : contre core_nt (ou core_nt_plants_fungi alias)
- Output: dossier outdir/ avec les MÊMES noms de fichiers que homologs/
          et on GARDE toutes les séquences d'origine (y compris la query).
- Enrichissement : on ajoute les séquences core_nt extraites via blastdbcmd,
                   avec headers formatés comme tes homologs.

NOUVEAU:
- Sélection BALANCÉE par espèces:
    --max-plant-species (def 50) / --max-fungi-species (def 50) / --max-other-species (def 50)
  Et par espèce:
    --max-hits-per-species (def 3)

Defaults demandés:
- --min-pident 70
- --evalue 1e-50
- 50 espèces plant / 50 espèces fungi / 50 espèces other
- hits per species = 3
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
    max_target_seqs: int,
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
        "-max_target_seqs", str(max_target_seqs),
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


def filter_and_select_hits_balanced_by_species(
    hits: List[Dict[str, str]],
    plants_taxids: Set[str],
    fungi_taxids: Set[str],
    min_pident: float,
    min_qcov: float,
    max_hits_per_species: int,
    max_plant_species: int,
    max_fungi_species: int,
    max_other_species: int,
    exclude_self: bool = True
) -> Dict[str, List[Dict[str, str]]]:
    """
    Pour chaque query:
      - filtre pident/qcov
      - trie par qualité (evalue asc, bitscore desc)
      - sélectionne au plus:
          * max_plant_species espèces plantes
          * max_fungi_species espèces champignons
          * max_other_species espèces other
        et au plus max_hits_per_species hits par espèce (dans chaque royaume)
    """
    by_q: Dict[str, List[Dict[str, str]]] = {}

    # Filtrage
    for h in hits:
        q = h.get("qseqid", "")
        if not q:
            continue
        if exclude_self and h.get("sseqid", "") == q:
            continue

        pident = parse_float(h.get("pident", ""))
        qcov = hit_qcov_pct(h)
        if pident is None or qcov is None:
            continue
        if pident < min_pident or qcov < min_qcov:
            continue

        by_q.setdefault(q, []).append(h)

    selected: Dict[str, List[Dict[str, str]]] = {}

    for q, qhits in by_q.items():
        qhits = sort_hits_by_quality(qhits)

        # espèces déjà “ouvertes” par royaume
        plant_species: Set[str] = set()
        fungi_species: Set[str] = set()
        other_species: Set[str] = set()

        # compte hits par espèce (par royaume)
        plant_counts: Dict[str, int] = {}
        fungi_counts: Dict[str, int] = {}
        other_counts: Dict[str, int] = {}

        out: List[Dict[str, str]] = []

        for h in qhits:
            tax = (h.get("staxids", "") or "").split(";")[0].strip()
            kind = kingdom_from_taxid(tax, plants_taxids, fungi_taxids)

            sci = sanitize_species_name(h.get("sscinames", ""))
            if sci == "unknown":
                sci = f"taxid_{tax}" if tax else "unknown"

            if kind == "plant":
                if sci not in plant_species and len(plant_species) >= max_plant_species:
                    continue
                if max_hits_per_species > 0 and plant_counts.get(sci, 0) >= max_hits_per_species:
                    continue
                plant_species.add(sci)
                plant_counts[sci] = plant_counts.get(sci, 0) + 1
                out.append(h)

            elif kind == "fungi":
                if sci not in fungi_species and len(fungi_species) >= max_fungi_species:
                    continue
                if max_hits_per_species > 0 and fungi_counts.get(sci, 0) >= max_hits_per_species:
                    continue
                fungi_species.add(sci)
                fungi_counts[sci] = fungi_counts.get(sci, 0) + 1
                out.append(h)

            else:
                if max_other_species <= 0:
                    continue
                if sci not in other_species and len(other_species) >= max_other_species:
                    continue
                if max_hits_per_species > 0 and other_counts.get(sci, 0) >= max_hits_per_species:
                    continue
                other_species.add(sci)
                other_counts[sci] = other_counts.get(sci, 0) + 1
                out.append(h)

            # arrêt rapide si tout est rempli et qu'on ne peut plus ouvrir de nouvelles espèces
            if (len(plant_species) >= max_plant_species and
                len(fungi_species) >= max_fungi_species and
                (len(other_species) >= max_other_species or max_other_species <= 0)):
                # on pourrait encore ajouter des hits pour espèces déjà ouvertes,
                # donc on ne break pas "agressivement" ici.
                pass

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
    Format comme tes homologs:
      >{query_file_stem}__{plant|fungi|other}__{Species}__{sseqid}:{start}-{end}(+/-)
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
# Longest query extraction
# -------------------------

def pick_longest_record(records: List[SeqRecord]) -> SeqRecord:
    if not records:
        raise ValueError("Empty FASTA records")
    best = records[0]
    best_len = len(best.seq)
    for r in records[1:]:
        L = len(r.seq)
        if L > best_len:
            best = r
            best_len = L
    return best


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Enrich homologs with core_nt from longest-per-file queries. Keeps ALL original sequences (including query) in final output."
    )

    ap.add_argument("--homologs-dir", required=True, type=Path,
                    help="Directory with existing homolog multifastas (e.g. Result_HT/homologs)")
    ap.add_argument("--pattern", default="*.fasta", help="Pattern (default: *.fasta)")

    ap.add_argument("--core-db", required=True, type=str,
                    help="BLAST DB name/path (e.g. core_nt_plants_fungi)")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="Output directory (files keep same names as in homologs-dir)")

    ap.add_argument("--plants-taxids", required=True, type=Path,
                    help="Taxid list for plant species (one taxid per line)")
    ap.add_argument("--fungi-taxids", required=True, type=Path,
                    help="Taxid list for fungi species (one taxid per line)")

    # BLAST parameters (defaults demandés)
    ap.add_argument("--blast-program", default="blastn")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-50)
    ap.add_argument("--max-target-seqs", type=int, default=5000,
                    help="blastn -max_target_seqs (default: 5000). Needs to be high enough to capture 50 species/kingdom.")

    # filters
    ap.add_argument("--min-pident", type=float, default=70.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)

    # selection by species (defaults demandés)
    ap.add_argument("--max-hits-per-species", type=int, default=3,
                    help="Max hits per species within each kingdom (default: 3). Use 0 to disable.")
    ap.add_argument("--max-plant-species", type=int, default=50)
    ap.add_argument("--max-fungi-species", type=int, default=50)
    ap.add_argument("--max-other-species", type=int, default=50)

    ap.add_argument("--exclude-self", action="store_true", default=True)
    ap.add_argument("--include-self", dest="exclude_self", action="store_false")

    ap.add_argument("--summary", type=Path, default=None,
                    help="Summary TSV (default: <outdir>/core_nt_from_longest_summary.tsv)")
    ap.add_argument("--keep-temp", action="store_true",
                    help="Keep temporary query fasta / blast tsv")

    args = ap.parse_args()

    if not args.homologs_dir.is_dir():
        sys.exit(f"[ERROR] homologs-dir not found: {args.homologs_dir}")
    if not args.plants_taxids.exists():
        sys.exit(f"[ERROR] plants-taxids not found: {args.plants_taxids}")
    if not args.fungi_taxids.exists():
        sys.exit(f"[ERROR] fungi-taxids not found: {args.fungi_taxids}")

    check_tool(args.blast_program)
    check_tool("blastdbcmd")

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary or (args.outdir / "core_nt_from_longest_summary.tsv")

    plants_set = load_taxid_set(args.plants_taxids)
    fungi_set = load_taxid_set(args.fungi_taxids)

    fasta_files = sorted(args.homologs_dir.glob(args.pattern))
    if not fasta_files:
        sys.exit(f"[ERROR] No FASTA files found in {args.homologs_dir} with pattern {args.pattern}")

    tmp_query_fasta = args.outdir / "_tmp_longest_queries.fasta"
    tmp_blast_tsv = args.outdir / "_tmp_core_nt_blast.tsv"

    q_to_file: Dict[str, Path] = {}
    queries: List[SeqRecord] = []

    for fp in fasta_files:
        stem = fp.stem
        recs = list(SeqIO.parse(str(fp), "fasta"))
        if not recs:
            eprint(f"[WARN] Empty fasta, skipping: {fp.name}")
            continue
        longest = pick_longest_record(recs)

        qrec = longest[:]  # copy
        qrec.id = stem
        qrec.name = stem
        qrec.description = ""
        queries.append(qrec)
        q_to_file[stem] = fp

    if not queries:
        sys.exit("[ERROR] No queries extracted (all files empty?).")

    SeqIO.write(queries, str(tmp_query_fasta), "fasta")
    print(f"[INFO] Queries built from longest sequences: {len(queries)}")
    print(f"[INFO] Temp query FASTA: {tmp_query_fasta}")

    # run BLAST (one single run)
    print("[INFO] Running BLAST against core DB ...")
    run_blast_multiquery(
        query_fasta=tmp_query_fasta,
        core_db=args.core_db,
        blast_program=args.blast_program,
        threads=args.threads,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        extra_args=[],
        out_tsv=tmp_blast_tsv
    )

    hits = load_blast_hits(tmp_blast_tsv)
    print(f"[INFO] Loaded raw hits: {len(hits)}")

    selected = filter_and_select_hits_balanced_by_species(
        hits=hits,
        plants_taxids=plants_set,
        fungi_taxids=fungi_set,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        max_hits_per_species=args.max_hits_per_species,
        max_plant_species=args.max_plant_species,
        max_fungi_species=args.max_fungi_species,
        max_other_species=args.max_other_species,
        exclude_self=args.exclude_self
    )

    with summary_path.open("w", encoding="utf-8") as s:
        s.write(
            "candidate_file\tqseqid\tstatus\tn_in_homologs\t"
            "n_selected_total\tn_added_core_nt\t"
            "plant_species\tfungi_species\tother_species\t"
            "reason\n"
        )

        for qseqid in sorted(q_to_file.keys()):
            src_fp = q_to_file[qseqid]
            out_fp = args.outdir / src_fp.name

            base_recs = list(SeqIO.parse(str(src_fp), "fasta"))
            if not base_recs:
                s.write(f"{src_fp.name}\t{qseqid}\tEMPTY\t0\t0\t0\t0\t0\t0\tempty_source_fasta\n")
                continue

            # ✅ On garde toutes les séquences d'origine (y compris la query)
            existing_ids = {r.id for r in base_recs}
            existing_seqs = {str(r.seq) for r in base_recs}

            sel_hits = selected.get(qseqid, [])

            # compter espèces par royaume dans la sélection
            pl_sp = set()
            fu_sp = set()
            ot_sp = set()
            for h in sel_hits:
                tax = (h.get("staxids", "") or "").split(";")[0].strip()
                k = kingdom_from_taxid(tax, plants_set, fungi_set)
                sp = sanitize_species_name(h.get("sscinames", ""))
                if sp == "unknown":
                    sp = f"taxid_{tax}" if tax else "unknown"
                if k == "plant":
                    pl_sp.add(sp)
                elif k == "fungi":
                    fu_sp.add(sp)
                else:
                    ot_sp.add(sp)

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

                rec = make_core_nt_record(
                    fasta_text=fasta_text,
                    query_file_stem=qseqid,
                    hit_row=h,
                    plants_taxids=plants_set,
                    fungi_taxids=fungi_set
                )
                if rec is None:
                    fetch_fail += 1
                    continue

                if rec.id in existing_ids:
                    continue
                if str(rec.seq) in existing_seqs:
                    continue

                rec.id = normalize_id_unique(rec.id, existing_ids)
                rec.name = rec.id
                existing_ids.add(rec.id)
                existing_seqs.add(str(rec.seq))
                added.append(rec)

            SeqIO.write(base_recs + added, str(out_fp), "fasta")

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
                f"{src_fp.name}\t{qseqid}\t{status}\t{len(base_recs)}\t"
                f"{len(sel_hits)}\t{len(added)}\t"
                f"{len(pl_sp)}\t{len(fu_sp)}\t{len(ot_sp)}\t"
                f"{reason}\n"
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
