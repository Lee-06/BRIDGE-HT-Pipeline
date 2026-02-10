#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from Bio import Phylo, SeqIO


# ---------------------------
# Helpers
# ---------------------------

def get_candidate_id_from_treefile(path: Path) -> str:
    name = path.name
    for suf in [
        ".trimmed.aln.contree", ".trimmed.aln.treefile",
        ".aln.contree", ".aln.treefile",
        ".contree", ".treefile",
    ]:
        if name.endswith(suf):
            return name[: -len(suf)]
    return path.stem


def get_candidate_id_from_homolog_fasta(path: Path) -> str:
    # In Step 10, cid = mf.stem, so we mirror that here.
    return path.stem


def confidence_of_clade(clade) -> Optional[float]:
    c = getattr(clade, "confidence", None)
    if c is None:
        return None
    try:
        return float(c)
    except Exception:
        return None


def try_midpoint_root(tree) -> bool:
    try:
        tree.root_at_midpoint()
        return True
    except Exception:
        return False


def sanitize_species(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\.(fasta|fa|fna|fas)$", "", s, flags=re.IGNORECASE)
    return s


def parse_species_from_leaf(name: str) -> str:
    """
    Extract species label from leaf header.
    Uses the token after "__plant__" or "__fungi__" up to next "__"
    and strips fasta extensions.
    """
    if not name:
        return "unknown"
    n = str(name)

    m = re.search(r"__(plant|fungi)__([^_].*?)__", n)
    if m:
        return sanitize_species(m.group(2))

    m2 = re.search(r"__(plant|fungi)__([^_].+)$", n)
    if m2:
        chunk = m2.group(2).split("__")[0]
        return sanitize_species(chunk)

    return "unknown"


def classify_leaf(name: str, candidate_kingdom: str) -> str:
    if not name:
        return "other"
    n = str(name)
    if "__plant__" in n:
        return "plant"
    if "__fungi__" in n:
        return "fungi"
    if n.startswith("HTcandidate_"):
        return candidate_kingdom
    return "other"


def count_and_species_in_clade(clade, leaf_kind: Dict[str, str]) -> Tuple[int, int, int]:
    plants = fungi = other = 0
    for t in clade.get_terminals():
        k = leaf_kind.get(t.name, "other")
        if k == "plant":
            plants += 1
        elif k == "fungi":
            fungi += 1
        else:
            other += 1
    return plants, fungi, other


def species_sets_from_terminals(terminals, leaf_kind: Dict[str, str]) -> Tuple[Set[str], Set[str], Set[str]]:
    plant_sp, fungi_sp, other_sp = set(), set(), set()
    for t in terminals:
        k = leaf_kind.get(t.name, "other")
        sp = parse_species_from_leaf(t.name)
        if k == "plant":
            plant_sp.add(sp)
        elif k == "fungi":
            fungi_sp.add(sp)
        else:
            other_sp.add(sp)
    plant_sp.discard("unknown")
    fungi_sp.discard("unknown")
    other_sp.discard("unknown")
    return plant_sp, fungi_sp, other_sp


def find_mrca(tree, leaf_names: List[str]):
    if not leaf_names:
        return None
    name_to_term = {t.name: t for t in tree.get_terminals()}
    terms = [name_to_term[n] for n in leaf_names if n in name_to_term]
    if not terms:
        return None
    if len(terms) == 1:
        return terms[0]
    return tree.common_ancestor(terms)


def clade_contains_kind(clade, leaf_kind: Dict[str, str], kind: str) -> bool:
    for t in clade.get_terminals():
        if leaf_kind.get(t.name) == kind:
            return True
    return False


def supported_mixed_clade_stats(tree, leaf_kind: Dict[str, str], thr: float) -> Tuple[int, Optional[float]]:
    """
    returns (count_supported_mixed_nodes, max_bootstrap)
    A mixed node has >=1 plant and >=1 fungi.
    """
    count = 0
    max_bs = None
    for clade in tree.get_nonterminals():
        bs = confidence_of_clade(clade)
        if bs is None or bs < thr:
            continue
        p, f, _ = count_and_species_in_clade(clade, leaf_kind)
        if p > 0 and f > 0:
            count += 1
            if max_bs is None or bs > max_bs:
                max_bs = bs
    return count, max_bs


def closest_plant_fungi_pair(tree, leaf_kind: Dict[str, str]) -> Tuple[str, str, str, str, str, str]:
    """
    Find closest plant leaf to fungi leaf by patristic distance.
    Returns:
      (plant_species, fungi_species, plant_leaf, fungi_leaf, min_distance_str, mrca_bootstrap_str)
    If not possible => "NA"
    """
    terminals = tree.get_terminals()
    plant_terms = [t for t in terminals if leaf_kind.get(t.name) == "plant"]
    fungi_terms = [t for t in terminals if leaf_kind.get(t.name) == "fungi"]

    if not plant_terms or not fungi_terms:
        return ("NA", "NA", "NA", "NA", "NA", "NA")

    best = None  # (dist, plant_term, fungi_term, mrca_bs)
    for p in plant_terms:
        for f in fungi_terms:
            try:
                d = tree.distance(p, f)
            except Exception:
                continue
            if d is None:
                continue
            mrca = tree.common_ancestor([p, f])
            bs = confidence_of_clade(mrca)
            if best is None or d < best[0]:
                best = (d, p, f, bs)

    if best is None:
        return ("NA", "NA", "NA", "NA", "NA", "NA")

    d, p, f, bs = best
    psp = parse_species_from_leaf(p.name)
    fsp = parse_species_from_leaf(f.name)
    bs_str = "NA" if bs is None else f"{bs:.2f}"
    return (psp, fsp, p.name, f.name, f"{d:.6f}", bs_str)


def fmt_species_ids(sps: Set[str]) -> str:
    if not sps:
        return ""
    return ",".join(sorted(sps))


# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Step 11: Summarize trees; optionally also report candidates with no phylogeny (<min-seqs) using homolog FASTAs."
    )
    ap.add_argument("--phylo-dir", required=True, type=Path, help="Directory with IQ-TREE outputs")
    ap.add_argument("--pattern", default="*.contree", help="Tree file pattern (default: *.contree)")
    ap.add_argument("--bootstrap-threshold", type=float, default=70.0, help="Bootstrap threshold (default: 70)")
    ap.add_argument("--candidate-kingdom", choices=["plant", "fungi", "other"], default="plant",
                    help="How to treat 'HTcandidate_*' leaf if present (default: plant)")
    ap.add_argument("--out", default="step11_phylo_summary.tsv", help="Output TSV")
    ap.add_argument("--keep-list", default="HT_keep_ids.txt", help="IDs to keep as HT candidates")
    ap.add_argument("--exclude-list", default="HT_exclude_ids.txt", help="IDs to exclude (conserved false positives)")
    ap.add_argument("--write-lists", action="store_true", default=False,
                    help="Write keep/exclude ID list files (disabled by default).")
    ap.add_argument("--midpoint-root", action="store_true", default=True, help="Apply midpoint rooting (default: ON)")
    ap.add_argument("--no-midpoint-root", dest="midpoint_root", action="store_false", help="Disable midpoint rooting")

    ap.add_argument("--min-fungi-species", type=int, default=2, help="Min fungi species for HT_SPECIFIC_PLANT (default: 2)")
    ap.add_argument("--min-plant-species", type=int, default=2, help="Min plant species for HT_SPECIFIC_FUNGI (default: 2)")

    # NEW: include candidates without trees (too few seqs)
    ap.add_argument("--homologs-dir", type=Path, default=None,
                    help="Directory with homolog multi-FASTAs (Step 9). If provided, candidates without trees are added to TSV.")
    ap.add_argument("--homologs-pattern", default="*.fasta",
                    help="Pattern for homolog multifastas in --homologs-dir (default: *.fasta)")
    ap.add_argument("--min-seqs", type=int, default=4,
                    help="Minimum sequences required to have a phylogeny (default: 4, same as Step 10).")

    args = ap.parse_args()

    if not args.phylo_dir.is_dir():
        sys.exit(f"[ERROR] phylo-dir not found: {args.phylo_dir}")

    thr = float(args.bootstrap_threshold)

    # Load tree files
    tree_files = sorted(args.phylo_dir.glob(args.pattern))
    if not tree_files:
        sys.exit(f"[ERROR] No tree files matching {args.pattern} in {args.phylo_dir}")

    # Map: candidate_id -> treefile path
    cid_to_tree: Dict[str, Path] = {}
    for fp in tree_files:
        cid = get_candidate_id_from_treefile(fp)
        cid_to_tree[cid] = fp

    rows = []
    keep, exclude = [], []

    # ---------- Process treefiles ----------
    for cid, fp in sorted(cid_to_tree.items()):
        try:
            tree = Phylo.read(str(fp), "newick")
        except Exception as e:
            rows.append({
                "candidate_id": cid,
                "status": "ERROR",
                "category": "ERROR",
                "reason": f"cannot_read_tree: {e}",
                "n_plants": 0, "n_fungi": 0, "n_other": 0,
                "n_plants_species": 0, "n_fungi_species": 0, "n_other_species": 0,
                "plant_species_ids": "",
                "fungi_species_ids": "",
                "other_species_ids": "",
                "plant_monophyly": "na",
                "fungi_monophyly": "na",
                "max_supported_mixed_bootstrap": "",
                "midpoint_rooted": "no",
                "closest_plant_species": "NA",
                "closest_fungi_species": "NA",
                "closest_plant_leaf": "NA",
                "closest_fungi_leaf": "NA",
                "closest_pair_distance": "NA",
                "closest_pair_mrca_bootstrap": "NA",
                "file": fp.name
            })
            continue

        midpoint_done = "no"
        if args.midpoint_root:
            if try_midpoint_root(tree):
                midpoint_done = "yes"

        terminals = tree.get_terminals()
        leaf_kind = {t.name: classify_leaf(t.name, args.candidate_kingdom) for t in terminals}

        plant_leaves = [n for n, k in leaf_kind.items() if k == "plant"]
        fungi_leaves = [n for n, k in leaf_kind.items() if k == "fungi"]
        other_leaves = [n for n, k in leaf_kind.items() if k == "other"]

        n_plants, n_fungi, n_other = len(plant_leaves), len(fungi_leaves), len(other_leaves)

        plant_sp, fungi_sp, other_sp = species_sets_from_terminals(terminals, leaf_kind)
        n_pl_sp, n_fu_sp, n_ot_sp = len(plant_sp), len(fungi_sp), len(other_sp)

        mixed_count, max_mixed_bs = supported_mixed_clade_stats(tree, leaf_kind, thr)

        plant_mrca = find_mrca(tree, plant_leaves)
        fungi_mrca = find_mrca(tree, fungi_leaves)

        plant_mono = (n_plants >= 2) and (plant_mrca is not None) and (not clade_contains_kind(plant_mrca, leaf_kind, "fungi"))
        fungi_mono = (n_fungi >= 2) and (fungi_mrca is not None) and (not clade_contains_kind(fungi_mrca, leaf_kind, "plant"))

        (cp_sp, cf_sp, cp_leaf, cf_leaf, cp_dist, cp_bs) = closest_plant_fungi_pair(tree, leaf_kind)

        # -------------------------
        # Decision logic (as-is)
        # -------------------------
        status = "AMBIGUOUS"
        category = "AMBIGUOUS_LOW_SUPPORT"
        reason = "default"

        if n_plants == 0 or n_fungi == 0:
            status = "INSUFFICIENT_TAXA"
            category = "INSUFFICIENT_TAXA"
            reason = "missing_plants_or_fungi"
        else:
            if n_pl_sp == 1 and n_fu_sp >= args.min_fungi_species:
                embedded = (plant_mrca is not None) and clade_contains_kind(plant_mrca, leaf_kind, "fungi")
                if embedded:
                    status = "KEEP"
                    category = "HT_SPECIFIC_PLANT"
                    reason = "single_plant_species_embedded_in_fungi_clade"
                else:
                    status = "AMBIGUOUS"
                    category = "AMBIGUOUS_LOW_SUPPORT"
                    reason = "single_plant_species_but_not_clearly_embedded"

            elif n_fu_sp == 1 and n_pl_sp >= args.min_plant_species:
                embedded = (fungi_mrca is not None) and clade_contains_kind(fungi_mrca, leaf_kind, "plant")
                if embedded:
                    status = "KEEP"
                    category = "HT_SPECIFIC_FUNGI"
                    reason = "single_fungi_species_embedded_in_plant_clade"
                else:
                    status = "AMBIGUOUS"
                    category = "AMBIGUOUS_LOW_SUPPORT"
                    reason = "single_fungi_species_but_not_clearly_embedded"

            elif mixed_count > 0:
                status = "KEEP"
                category = "HT_SUPPORTED_CLASSICAL"
                reason = f"mixed_plant_fungi_clade_supported_bootstrap>={thr}"

            elif plant_mono and fungi_mono:
                status = "CONGRUENT"
                category = "CONGRUENT_SEPARATED"
                reason = "plants_and_fungi_monophyletic_no_supported_mixing"

            elif (n_pl_sp >= 3 and n_fu_sp >= 3 and mixed_count >= 2):
                status = "REJECT"
                category = "CONSERVE_FALSE_POSITIVE"
                reason = f"multiple_supported_mixed_nodes_bootstrap>={thr}"

            else:
                status = "AMBIGUOUS"
                category = "AMBIGUOUS_LOW_SUPPORT"
                reason = f"no_supported_mixing_bootstrap>={thr}"

        if status == "KEEP":
            keep.append(cid)
        if status == "REJECT":
            exclude.append(cid)

        rows.append({
            "candidate_id": cid,
            "status": status,
            "category": category,
            "reason": reason,
            "n_plants": n_plants,
            "n_fungi": n_fungi,
            "n_other": n_other,
            "n_plants_species": n_pl_sp,
            "n_fungi_species": n_fu_sp,
            "n_other_species": n_ot_sp,
            "plant_species_ids": fmt_species_ids(plant_sp),
            "fungi_species_ids": fmt_species_ids(fungi_sp),
            "other_species_ids": fmt_species_ids(other_sp),
            "plant_monophyly": "yes" if plant_mono else "no" if n_plants >= 2 else "na",
            "fungi_monophyly": "yes" if fungi_mono else "no" if n_fungi >= 2 else "na",
            "max_supported_mixed_bootstrap": "" if max_mixed_bs is None else f"{max_mixed_bs:.2f}",
            "midpoint_rooted": midpoint_done,
            "closest_plant_species": cp_sp,
            "closest_fungi_species": cf_sp,
            "closest_plant_leaf": cp_leaf,
            "closest_fungi_leaf": cf_leaf,
            "closest_pair_distance": cp_dist,
            "closest_pair_mrca_bootstrap": cp_bs,
            "file": fp.name
        })

    # ---------- Add "no phylogeny" candidates from homolog FASTAs ----------
    if args.homologs_dir is not None:
        if not args.homologs_dir.is_dir():
            sys.exit(f"[ERROR] homologs-dir not found: {args.homologs_dir}")

        homolog_fastas = sorted(args.homologs_dir.glob(args.homologs_pattern))
        if not homolog_fastas:
            sys.exit(f"[ERROR] No homolog FASTAs matching {args.homologs_pattern} in {args.homologs_dir}")

        already = {r["candidate_id"] for r in rows}

        for mf in homolog_fastas:
            cid = get_candidate_id_from_homolog_fasta(mf)
            if cid in already:
                continue  # has tree summary already

            records = list(SeqIO.parse(str(mf), "fasta"))
            n_records = len(records)

            # Only add if it's below threshold OR if you want to add ALL missing-tree cases
            if n_records >= args.min_seqs:
                # This means: no treefile found even though enough seqs -> could be failure.
                status = "NO_PHYLOGENY"
                category = "NO_PHYLOGENY"
                reason = f"no_treefile_found_but_nseqs={n_records}"
            else:
                status = "NO_PHYLOGENY"
                category = "NO_PHYLOGENY"
                reason = f"no_phylogeny_less_than_{args.min_seqs}_sequences"

            # Build fake "terminals" from FASTA headers
            # We reuse the same parsing rules as trees (based on header strings).
            leaf_names = [r.id for r in records]
            leaf_kind = {n: classify_leaf(n, args.candidate_kingdom) for n in leaf_names}

            plant_leaves = [n for n, k in leaf_kind.items() if k == "plant"]
            fungi_leaves = [n for n, k in leaf_kind.items() if k == "fungi"]
            other_leaves = [n for n, k in leaf_kind.items() if k == "other"]

            n_plants, n_fungi, n_other = len(plant_leaves), len(fungi_leaves), len(other_leaves)

            # Species sets from FASTA "terminals"
            plant_sp, fungi_sp, other_sp = set(), set(), set()
            for n in leaf_names:
                k = leaf_kind.get(n, "other")
                sp = parse_species_from_leaf(n)
                if k == "plant":
                    plant_sp.add(sp)
                elif k == "fungi":
                    fungi_sp.add(sp)
                else:
                    other_sp.add(sp)
            plant_sp.discard("unknown")
            fungi_sp.discard("unknown")
            other_sp.discard("unknown")

            rows.append({
                "candidate_id": cid,
                "status": status,
                "category": category,
                "reason": reason,
                "n_plants": n_plants,
                "n_fungi": n_fungi,
                "n_other": n_other,
                "n_plants_species": len(plant_sp),
                "n_fungi_species": len(fungi_sp),
                "n_other_species": len(other_sp),
                "plant_species_ids": fmt_species_ids(plant_sp),
                "fungi_species_ids": fmt_species_ids(fungi_sp),
                "other_species_ids": fmt_species_ids(other_sp),
                "plant_monophyly": "na",
                "fungi_monophyly": "na",
                "max_supported_mixed_bootstrap": "",
                "midpoint_rooted": "no",
                "closest_plant_species": "NA",
                "closest_fungi_species": "NA",
                "closest_plant_leaf": "NA",
                "closest_fungi_leaf": "NA",
                "closest_pair_distance": "NA",
                "closest_pair_mrca_bootstrap": "NA",
                "file": mf.name
            })

    out = Path(args.out)

    headers = [
        "candidate_id", "status", "category", "reason",
        "n_plants", "n_fungi", "n_other",
        "n_plants_species", "n_fungi_species", "n_other_species",
        "plant_species_ids", "fungi_species_ids", "other_species_ids",
        "plant_monophyly", "fungi_monophyly",
        "max_supported_mixed_bootstrap",
        "midpoint_rooted",
        "closest_plant_species", "closest_fungi_species",
        "closest_plant_leaf", "closest_fungi_leaf",
        "closest_pair_distance", "closest_pair_mrca_bootstrap",
        "file",
    ]

    # Sort rows by candidate_id for readability
    rows_sorted = sorted(rows, key=lambda r: r.get("candidate_id", ""))

    with out.open("w", encoding="utf-8") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows_sorted:
            f.write("\t".join(str(r.get(h, "")) for h in headers) + "\n")

    # Write keep/exclude lists only if explicitly requested
    if args.write_lists:
        Path(args.keep_list).write_text("\n".join(sorted(set(keep))) + "\n", encoding="utf-8")
        Path(args.exclude_list).write_text("\n".join(sorted(set(exclude))) + "\n", encoding="utf-8")
    

    counts = {}
    for r in rows_sorted:
        counts[r["category"]] = counts.get(r["category"], 0) + 1

    print(f"[SUCCESS] Summary TSV: {out}")
    if args.write_lists:
        print(f"[SUCCESS] KEEP IDs   : {args.keep_list}")
        print(f"[SUCCESS] REJECT IDs : {args.exclude_list}")
    else:
        print("[INFO] keep/reject list files were not written (use --write-lists to enable).")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
