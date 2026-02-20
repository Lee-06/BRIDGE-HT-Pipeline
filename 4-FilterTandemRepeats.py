#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import re
from pathlib import Path
from typing import Dict, List, Tuple


# ---------------- FASTA IO ----------------

def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_id = None
    buf: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(buf)
    return seqs


def write_fasta(records: List[Tuple[str, str]], out: Path, wrap: int = 60) -> None:
    with out.open("w", encoding="utf-8") as f:
        for sid, seq in records:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")


def concat_selected_fastas(selected_dir: Path, out_fasta: Path) -> int:
    fasta_exts = {".fasta", ".fa", ".fna"}
    used = 0
    with out_fasta.open("w", encoding="utf-8") as out:
        for f in sorted(selected_dir.iterdir()):
            if not f.is_file():
                continue
            if f.suffix.lower() not in fasta_exts:
                continue
            try:
                text = f.read_text(encoding="utf-8", errors="replace").strip()
            except Exception:
                continue
            if not text:
                continue
            out.write(text + "\n")
            used += 1
    return used


# ---------------- TRF UTILS ----------------

def run_trf_in_dir(fasta_in_outdir_name: str, trf_path: str, params: List[str], cwd: Path) -> None:
    """
    Run TRF with cwd forcing outputs (.mask/.dat) into cwd.
    IMPORTANT: fasta arg must be relative or just the filename inside cwd.
    """
    cmd = [trf_path, fasta_in_outdir_name] + params
    print("[INFO] Running TRF in:", str(cwd))
    print("       " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd))


def normalize_mask(seq: str) -> str:
    seq = re.sub(r"[acgt]", "N", seq)
    seq = seq.replace("n", "N")
    return seq


def split_on_N(seq: str, min_run: int) -> List[str]:
    splitter = re.compile(rf"N{{{min_run},}}", re.IGNORECASE)
    return [p for p in splitter.split(seq) if len(p) > 0 and p.upper() != "N"]


# ---------------- MAIN ----------------

def main():
    ap = argparse.ArgumentParser(
        description="Script 4: Clean candidate sequences (output of Script 3) using TRF masking."
    )

    ap.add_argument("--selected-dir", required=True,
                    help="Directory produced by Script 3 (selected_sequences/) containing selected_*.fasta")
    ap.add_argument("--outdir", default="trf_clean",
                    help="Output directory (default: trf_clean)")

    ap.add_argument("--mode", choices=["hardmask", "remove", "split_longest"], default="split_longest")
    ap.add_argument("--min-len", type=int, default=200, help="Minimum length kept (default: 200)")
    ap.add_argument("--min-N-run", type=int, default=10, help="N-run size to split on (default: 10)")
    ap.add_argument("--max-mask-pct", type=float, default=90.0, help="Drop seqs masked >= this % (default: 90.0)")

    ap.add_argument("--trf-path", default="trf", help="Path to TRF binary (default: trf)")
    ap.add_argument("--trf-params", default="2,7,7,80,10,50,500,-f,-d,-m,-h",
                    help="TRF parameters as comma-separated list")
    ap.add_argument("--keep-trf-files", action="store_true",
                    help="Do not delete TRF intermediate files (.dat/.mask)")

    args = ap.parse_args()

    selected_dir = Path(args.selected_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not selected_dir.is_dir():
        raise SystemExit(f"[ERROR] selected-dir not found: {selected_dir}")

    # Master FASTA INSIDE outdir
    input_fasta = outdir / "ht_candidates.fasta"
    used_files = concat_selected_fastas(selected_dir, input_fasta)
    if used_files == 0:
        raise SystemExit(f"[ERROR] No FASTA files found (or all empty) in {selected_dir}")

    print(f"[INFO] Consolidated from {used_files} FASTA files into: {input_fasta}")

    output_fasta = outdir / "ht_candidates.cleaned.fasta"
    report_tsv = outdir / "ht_candidates.trf_report.tsv"

    trf_params = [x.strip() for x in args.trf_params.split(",") if x.strip()]

    # Infer signature from first 7 numeric params (matches TRF naming)
    numeric = []
    for x in trf_params:
        if x.startswith("-"):
            continue
        if re.fullmatch(r"\d+", x):
            numeric.append(x)
        if len(numeric) >= 7:
            break
    if len(numeric) < 7:
        raise SystemExit("[ERROR] Need first 7 numeric TRF params to infer output signature.")

    signature = "." + ".".join(numeric[:7])

    # Run TRF with cwd=outdir so outputs land IN outdir
    try:
        run_trf_in_dir(input_fasta.name, args.trf_path, trf_params, cwd=outdir)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] TRF binary not found: {args.trf_path}")
    except subprocess.CalledProcessError:
        raise SystemExit("[ERROR] TRF failed (check TRF installation/params).")

    mask_fasta = outdir / f"{input_fasta.name}{signature}.mask"
    dat_file = outdir / f"{input_fasta.name}{signature}.dat"

    if not mask_fasta.exists():
        raise SystemExit(f"[ERROR] TRF mask file not found in outdir: {mask_fasta}\n"
                         f"Hint: your TRF outputs are currently in the directory where you launched the script. "
                         f"This fixed version forces TRF outputs into outdir.")

    # Process masked FASTA
    seqs = read_fasta(mask_fasta)

    out_records: List[Tuple[str, str]] = []
    report_lines = ["id\tlen_in\tmasked_bp\tmasked_pct\tmode\tlen_out\tdecision"]

    for sid, raw in seqs.items():
        s = normalize_mask(raw)
        L = len(s)
        masked = s.count("N")
        pct = (masked / L * 100) if L else 0.0

        if pct >= args.max_mask_pct:
            report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\t{args.mode}\t0\tDROP_masked")
            continue

        if args.mode == "hardmask":
            if L >= args.min_len:
                out_records.append((sid, s))
                report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\thardmask\t{L}\tKEEP")
            else:
                report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\thardmask\t0\tDROP_len")

        elif args.mode == "remove":
            cleaned = s.replace("N", "")
            if len(cleaned) >= args.min_len:
                out_records.append((sid, cleaned))
                report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\tremove\t{len(cleaned)}\tKEEP")
            else:
                report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\tremove\t0\tDROP_len")

        else:  # split_longest
            parts = split_on_N(s, args.min_N_run)
            parts = [p for p in parts if len(p) >= args.min_len]
            if not parts:
                report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\tsplit_longest\t0\tDROP_no_fragment")
                continue
            best = max(parts, key=len)
            out_records.append((sid, best))
            report_lines.append(f"{sid}\t{L}\t{masked}\t{pct:.2f}\tsplit_longest\t{len(best)}\tKEEP")

    write_fasta(out_records, output_fasta)

    with report_tsv.open("w", encoding="utf-8") as f:
        f.write("\n".join(report_lines) + "\n")

    # Cleanup
    if not args.keep_trf_files:
        for p in [dat_file, mask_fasta]:
            try:
                if p.exists():
                    p.unlink()
            except Exception:
                pass

    print(f"[SUCCESS] Cleaned FASTA : {output_fasta}")
    print(f"[SUCCESS] Report       : {report_tsv}")
    print(f"[SUCCESS] Sequences kept: {len(out_records)}")


if __name__ == "__main__":
    main()
