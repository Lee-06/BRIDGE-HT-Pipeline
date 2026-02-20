#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

TARGET_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
DELIMS = {":", ",", ")", ";"}

FALLBACK_NAME_CLASSES = {
    "synonym",
    "equivalent name",
    "includes",
    "authority",
    "genbank synonym",
    "genbank common name",
    "common name",
}

GENUS_PLACEHOLDERS = {"sp", "sp.", "spp", "spp."}
GROUPS = {"fungi", "plant", "animal", "bacteria", "archaea", "protist", "virus"}


# -------------------------
# Newick leaf extraction
# -------------------------
def extract_leaf_labels(newick: str) -> List[str]:
    labels: List[str] = []
    i, n = 0, len(newick)
    while i < n:
        ch = newick[i]
        if ch in "(,":
            i += 1
            start = i
            while i < n and newick[i] not in DELIMS:
                i += 1
            lab = newick[start:i].strip()
            if lab:
                labels.append(lab)
            continue
        i += 1
    return labels


def normalize_spaces(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip())


# -------------------------
# Extract taxon from your IDs
# -------------------------
def taxon_from_leaf_id(label: str) -> str:
    """
    From:
      ...__fungi__Achaetomium_strumarium__scaffold...
      ...__plant__Chlorella_variabilis__gi|...|...
      ...__fungi__Xylariaceae_sp__scaffold...
      ...__plant__Magnoliophyta_environmental_sample__gi|...|...
    return a query:
      - "Genus species" or "Genus" (if sp/spp)
      - or "" if cannot parse
    """
    parts = label.split("__")
    if len(parts) < 3:
        return ""

    idx = None
    for i, p in enumerate(parts):
        if p in GROUPS:
            idx = i
            break
    if idx is None or idx + 1 >= len(parts):
        return ""

    token = parts[idx + 1].strip()
    if not token:
        return ""

    token_sp = normalize_spaces(token.replace("_", " "))
    fields = token_sp.split(" ")
    if not fields:
        return ""

    # handle "Magnoliophyta environmental sample"
    if len(fields) >= 2 and fields[1].lower() in {"environmental", "environmental_sample", "sample"}:
        return fields[0]

    genus = fields[0]
    if len(fields) == 1:
        return genus

    epithet = fields[1]
    el = epithet.lower()

    # sp / spp cases: "Ophiostoma_sp", "Xylariaceae_sp", "Diaporthaceae_sp"
    if el in GENUS_PLACEHOLDERS or el.startswith("sp-") or el.startswith("spp-"):
        return genus

    # strain-like suffix inside epithet: subellipsoidea-C-169, etc.
    if "-" in epithet and (epithet.count("-") >= 2 or re.search(r"-\d", epithet)):
        base = epithet.split("-")[0]
        return f"{genus} {base}"

    return f"{genus} {epithet}"


def candidate_queries(raw: str) -> List[str]:
    """
    Lookup fallbacks in names.dmp:
      - exact
      - if looks strain-like: base epithet
      - genus-only
    """
    s = normalize_spaces(raw)
    out = [s]
    toks = s.split(" ")
    if len(toks) >= 2:
        genus, second = toks[0], toks[1]
        if second.lower().startswith("sp"):
            out.append(genus)
            return dedup(out)
        if "-" in second and (second.count("-") >= 2 or re.search(r"-\d", second)):
            out.append(f"{genus} {second.split('-')[0]}")
        out.append(genus)
        return dedup(out)
    return dedup(out)


def dedup(items: List[str]) -> List[str]:
    seen = set()
    out = []
    for x in items:
        k = x.lower()
        if x and k not in seen:
            seen.add(k)
            out.append(x)
    return out


# -------------------------
# Taxdump parsing + lookup
# -------------------------
def parse_names_dmp_multi(names_path: Path) -> Dict[str, List[Tuple[int, str]]]:
    idx: Dict[str, List[Tuple[int, str]]] = {}
    with names_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            taxid_s, name_txt, _, name_class = parts[0], parts[1], parts[2], parts[3]
            try:
                taxid = int(taxid_s)
            except ValueError:
                continue
            idx.setdefault(name_txt.lower(), []).append((taxid, name_class))
    return idx


def parse_nodes_dmp(nodes_path: Path) -> Tuple[Dict[int, int], Dict[int, str]]:
    parent: Dict[int, int] = {}
    rank: Dict[int, str] = {}
    with nodes_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[0]); par = int(parts[1])
            except ValueError:
                continue
            parent[taxid] = par
            rank[taxid] = parts[2]
    return parent, rank


def parse_merged_dmp(merged_path: Path) -> Dict[int, int]:
    mp: Dict[int, int] = {}
    if not merged_path.exists():
        return mp
    with merged_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 2:
                continue
            try:
                old = int(parts[0]); new = int(parts[1])
            except ValueError:
                continue
            mp[old] = new
    return mp


def resolve_merged(taxid: int, merged: Dict[int, int]) -> int:
    seen = set()
    cur = taxid
    while cur in merged and cur not in seen:
        seen.add(cur)
        cur = merged[cur]
    return cur


def lineage_table(taxid: int, parent: Dict[int, int], rank: Dict[int, str]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    cur = taxid
    seen = set()
    while True:
        if cur in seen:
            break
        seen.add(cur)
        r = rank.get(cur, "")
        if r in TARGET_RANKS and r not in out:
            out[r] = cur
        par = parent.get(cur, cur)
        if par == cur:
            break
        cur = par
    return out


def closest_ancestor_with_rank(taxid: int, target_rank: str, parent: Dict[int, int], rank: Dict[int, str]) -> Optional[int]:
    cur = taxid
    seen = set()
    while True:
        if cur in seen:
            return None
        seen.add(cur)
        if rank.get(cur, "") == target_rank:
            return cur
        par = parent.get(cur, cur)
        if par == cur:
            return None
        cur = par


def score_taxid(taxid: int, parent: Dict[int, int], rank: Dict[int, str]) -> int:
    r = rank.get(taxid, "")
    sc = 0
    if r == "species":
        sc += 100
    elif closest_ancestor_with_rank(taxid, "species", parent, rank) is not None:
        sc += 70
    lin = lineage_table(taxid, parent, rank)
    if "genus" in lin:
        sc += 20
    if "family" in lin:
        sc += 10
    return sc


def choose_best(candidates: List[int], parent: Dict[int, int], rank: Dict[int, str]) -> int:
    best = candidates[0]
    best_score = -1
    for t in candidates:
        sc = score_taxid(t, parent, rank)
        if sc > best_score:
            best = t
            best_score = sc
    return best


def find_taxid_for_query(
    raw_query: str,
    name_index: Dict[str, List[Tuple[int, str]]],
    parent: Dict[int, int],
    rank: Dict[int, str],
    merged: Dict[int, int],
    allow_fallback_classes: bool = True,
) -> Tuple[str, str]:
    for q in candidate_queries(raw_query):
        hits = name_index.get(q.lower(), [])
        if not hits:
            continue

        sci = [resolve_merged(tid, merged) for (tid, cls) in hits if cls == "scientific name"]
        if sci:
            return (q, str(choose_best(sci, parent, rank)))

        if allow_fallback_classes:
            fb = [resolve_merged(tid, merged) for (tid, cls) in hits if cls in FALLBACK_NAME_CLASSES]
            if fb:
                return (q, str(choose_best(fb, parent, rank)))

    return ("", "")


def build_taxid_to_sciname(names_path: Path, needed: Set[int]) -> Dict[int, str]:
    out: Dict[int, str] = {}
    if not needed:
        return out
    with names_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            if parts[3] != "scientific name":
                continue
            try:
                tid = int(parts[0])
            except ValueError:
                continue
            if tid in needed and tid not in out:
                out[tid] = parts[1]
    return out


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="Extract unique taxa from tree leaves and build NCBI taxonomy TSV (local taxdump).")
    ap.add_argument("-i", "--trees-dir", required=True)
    ap.add_argument("-d", "--taxdump-dir", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--pattern", default="*.treefile")
    ap.add_argument("--no-fallback-names", action="store_true")
    ap.add_argument("--no-lift-to-species", action="store_true")
    ap.add_argument("--debug", action="store_true", help="Print a few parsed leaves/taxa for sanity-check")
    args = ap.parse_args()

    trees_dir = Path(args.trees_dir).expanduser().resolve()
    taxdir = Path(args.taxdump_dir).expanduser().resolve()
    out_path = Path(args.output).expanduser().resolve()

    names_path = taxdir / "names.dmp"
    nodes_path = taxdir / "nodes.dmp"
    merged_path = taxdir / "merged.dmp"
    if not trees_dir.is_dir():
        raise SystemExit(f"Trees dir not found: {trees_dir}")
    if not names_path.exists() or not nodes_path.exists():
        raise SystemExit("taxdump dir must contain names.dmp and nodes.dmp")

    files = sorted(trees_dir.glob(args.pattern))
    if not files:
        raise SystemExit(f"No files matched pattern '{args.pattern}' in {trees_dir}")

    print("Indexing names.dmp (all name classes; can take a bit)...")
    name_index = parse_names_dmp_multi(names_path)
    print("Reading nodes.dmp...")
    parent, rank = parse_nodes_dmp(nodes_path)
    merged = parse_merged_dmp(merged_path)

    allow_fallback = not args.no_fallback_names
    lift_to_species = not args.no_lift_to_species

    # 1) collect unique taxa (species or genus for sp/spp)
    taxa_set: Set[str] = set()
    debug_examples = 0

    for fp in files:
        txt = fp.read_text(encoding="utf-8", errors="replace")
        newick = "".join(line.strip() for line in txt.splitlines() if line.strip())
        labels = extract_leaf_labels(newick)

        if args.debug and labels and debug_examples < 5:
            print("EXAMPLE LEAF:", labels[0])
            debug_examples += 1

        for leaf in labels:
            taxon = taxon_from_leaf_id(leaf)
            if taxon:
                taxa_set.add(taxon)

    taxa = sorted(taxa_set, key=lambda x: x.lower())
    print(f"Found {len(taxa)} unique taxa (species or genus).")

    # 2) resolve taxonomy
    needed_taxids: Set[int] = set()

    # output columns: taxon (from leaves) + ncbi species + ranks (names) ; NA if missing
    rows: List[Dict[str, str]] = []

    for taxon in taxa:
        matched_q, matched_taxid_s = find_taxid_for_query(taxon, name_index, parent, rank, merged, allow_fallback_classes=allow_fallback)

        # Defaults = NA
        out = {
            "leaf_taxon": taxon,
            "ncbi_species": "NA",
            "taxid": "NA",
            "matched_query": matched_q if matched_q else "NA",
        }
        for r in TARGET_RANKS:
            out[r] = "NA"

        if matched_taxid_s.isdigit():
            matched_taxid = int(matched_taxid_s)
            base = matched_taxid

            if lift_to_species and rank.get(matched_taxid, "") != "species":
                cs = closest_ancestor_with_rank(matched_taxid, "species", parent, rank)
                if cs:
                    base = cs

            lin = lineage_table(base, parent, rank)

            # Fill rank names where possible
            for r in TARGET_RANKS:
                tid = lin.get(r)
                if tid:
                    out[r] = str(tid)  # temporarily store taxid; convert to name later
                    needed_taxids.add(tid)

            # species output
            sp_tid = lin.get("species")
            if sp_tid:
                out["taxid"] = str(sp_tid)
                needed_taxids.add(sp_tid)

        rows.append(out)

    # 3) taxid -> scientific name
    taxid_to_name = build_taxid_to_sciname(names_path, needed_taxids)

    # 4) write TSV
    header = ["leaf_taxon", "ncbi_species", "taxid"] + TARGET_RANKS

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            # Convert stored rank taxids to names; keep NA otherwise
            line = [row["leaf_taxon"]]

            # ncbi_species = name for species taxid
            if row["taxid"].isdigit():
                row["ncbi_species"] = taxid_to_name.get(int(row["taxid"]), "NA") or "NA"
            else:
                row["ncbi_species"] = "NA"

            line += [row["ncbi_species"], row["taxid"]]

            for r in TARGET_RANKS:
                v = row[r]
                if v.isdigit():
                    line.append(taxid_to_name.get(int(v), "NA") or "NA")
                else:
                    line.append("NA")

            f.write("\t".join(line) + "\n")

    print(f"Done: {out_path}")


if __name__ == "__main__":
    main()
