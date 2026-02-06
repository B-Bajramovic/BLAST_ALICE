#!/usr/bin/env python3
"""
make_blastp_sbatch.py

Create one Slurm sbatch file per FASTA file in an input directory, for BLASTP vs nr,
optionally restricted to a taxid.

Run:
  python make_blastp_sbatch.py -h
to see usage examples.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shlex
import subprocess
import sys


DEFAULT_DB = "/zfsstore/databases/NCBI/nr/nr"
DEFAULT_OUTFMT = "6 qseqid sacc pident length qcovs evalue bitscore stitle"

FASTA_SUFFIXES = {
    ".fa", ".faa", ".fasta", ".fsa", ".fas", ".pep"
}


def is_fasta_file(p: Path) -> bool:
    if not p.is_file():
        return False
    return p.suffix.lower() in FASTA_SUFFIXES


def sanitize_job_name(name: str) -> str:
    # Slurm job names are happier with simple characters
    keep = []
    for ch in name:
        if ch.isalnum() or ch in ("_", ".", "+", "="):
            keep.append(ch)
        else:
            keep.append("_")
    s = "".join(keep).strip("_")
    return s[:200] if s else "blastp_job"


def format_mem_gb(threads: int, mem_per_cpu_gb: int) -> int:
    return int(threads) * int(mem_per_cpu_gb)


def parse_extra_blast_args(raw: list[str]) -> list[str]:
    """
    Accept either repeated flags like:
      --extra -qcov_hsp_perc 30 --extra -word_size 2
    or a single quoted chunk:
      --extra "-qcov_hsp_perc 30 -word_size 2"
    """
    out: list[str] = []
    for item in raw:
        out.extend(shlex.split(item))
    return out


def write_sbatch(
    sbatch_path: Path,
    *,
    job_name: str,
    partition: str,
    time_limit: str,
    threads: int,
    mem_gb: int,
    mail_user: str | None,
    mail_type: str | None,
    query_fasta: Path,
    out_tsv: Path,
    db: str,
    taxid: int | None,
    evalue: str,
    max_target_seqs: int,
    outfmt: str,
    extra_blast_args: list[str],
    logs_dir: Path,
) -> None:
    blast_cmd = [
        "blastp",
        "-query", str(query_fasta),
        "-db", db,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-out", str(out_tsv),
        "-outfmt", outfmt,
    ]
    if taxid is not None:
        blast_cmd.extend(["-taxids", str(taxid)])
    blast_cmd.extend(extra_blast_args)

    blast_cmd_str = " ".join(shlex.quote(x) for x in blast_cmd)

    lines: list[str] = []
    lines.append("#!/bin/bash")
    lines.append(f"#SBATCH --job-name={job_name}")
    lines.append(f"#SBATCH --partition={partition}")
    lines.append(f"#SBATCH --time={time_limit}")
    lines.append(f"#SBATCH --cpus-per-task={threads}")
    lines.append(f"#SBATCH --mem={mem_gb}G")
    lines.append(f"#SBATCH --output={logs_dir.as_posix()}/%x_%j.out")
    lines.append(f"#SBATCH --error={logs_dir.as_posix()}/%x_%j.err")

    if mail_user:
        lines.append(f"#SBATCH --mail-user={mail_user}")
        if mail_type:
            lines.append(f"#SBATCH --mail-type={mail_type}")

    lines.append("")
    lines.append("set -euo pipefail")
    lines.append("")
    lines.append(f"mkdir -p {shlex.quote(logs_dir.as_posix())}")
    lines.append(f"mkdir -p {shlex.quote(out_tsv.parent.as_posix())}")
    lines.append("")
    lines.append('echo "Host: $(hostname)"')
    lines.append('echo "Start: $(date)"')
    lines.append('echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"')
    lines.append("")
    lines.append("# If your cluster uses environment modules, uncomment and adjust:")
    lines.append("# module load BLAST+/2.14.1")
    lines.append("")
    lines.append("# NCBI BLAST may require BLASTDB or taxdb availability depending on setup.")
    lines.append("# If needed, uncomment and adjust one of these patterns:")
    lines.append("# export BLASTDB=/zfsstore/databases/NCBI/nr")
    lines.append("# export BLASTDB=/zfsstore/databases/NCBI")
    lines.append("")
    lines.append(blast_cmd_str)
    lines.append("")
    lines.append('echo "End: $(date)"')
    lines.append("")

    sbatch_path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Generate sbatch files for BLASTP (one per FASTA file) in a directory.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Usage examples:

# basic run (defaults: taxid 28211, threads 32, mem per cpu 2GB, cpu-short, 04:00:00)
python make_blastp_sbatch.py /path/to/input_fastas

# choose taxid
python make_blastp_sbatch.py /path/to/input_fastas --taxid 9606

# disable taxid filtering
python make_blastp_sbatch.py /path/to/input_fastas --taxid 0

# change resources
python make_blastp_sbatch.py /path/to/input_fastas --threads 16 --mem-per-cpu-gb 2 --time 02:00:00

# add extra blast flags (repeatable)
python make_blastp_sbatch.py /path/to/input_fastas --taxid 28211 \\
  --extra "-qcov_hsp_perc 30" \\
  --extra "-word_size 2 -matrix BLOSUM45"

# generate and submit immediately
python make_blastp_sbatch.py /path/to/input_fastas --taxid 28211 --submit
"""
    )

    ap.add_argument(
        "inputdir",
        type=Path,
        help="Directory containing FASTA files (*.fa, *.faa, *.fasta, etc.)",
    )
    ap.add_argument(
        "--taxid",
        type=int,
        default=28211,
        help="NCBI taxid to restrict results (default: 28211). Use 0 to disable filtering.",
    )
    ap.add_argument(
        "--db",
        type=str,
        default=DEFAULT_DB,
        help=f"BLAST database prefix (default: {DEFAULT_DB})",
    )
    ap.add_argument(
        "--partition",
        type=str,
        default="cpu-short",
        help="Slurm partition (default: cpu-short)",
    )
    ap.add_argument(
        "--time",
        dest="time_limit",
        type=str,
        default="04:00:00",
        help="Walltime limit (default: 04:00:00)",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=32,
        help="Threads for BLASTP and Slurm cpus per task (default: 32)",
    )
    ap.add_argument(
        "--mem-per-cpu-gb",
        type=int,
        default=2,
        help="Memory per cpu in GB; total mem is threads * this (default: 2)",
    )
    ap.add_argument(
        "--evalue",
        type=str,
        default="1e-5",
        help="BLASTP evalue (default: 1e-5)",
    )
    ap.add_argument(
        "--max-target-seqs",
        type=int,
        default=20,
        help="BLASTP max_target_seqs (default: 20)",
    )
    ap.add_argument(
        "--outfmt",
        type=str,
        default=DEFAULT_OUTFMT,
        help=f'BLAST outfmt string (default: "{DEFAULT_OUTFMT}")',
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory for TSV and sbatch files (default: inputdir/blastp_out)",
    )
    ap.add_argument(
        "--logs-dirname",
        type=str,
        default="logs",
        help="Logs folder name inside outdir (default: logs)",
    )
    ap.add_argument(
        "--sbatch-dirname",
        type=str,
        default="sbatch",
        help="Sbatch folder name inside outdir (default: sbatch)",
    )
    ap.add_argument(
        "--mail-user",
        type=str,
        default=None,
        help="Email address for Slurm notifications (default: none)",
    )
    ap.add_argument(
        "--mail-type",
        type=str,
        default="BEGIN,END,FAIL",
        help="Slurm mail types if mail user is set (default: BEGIN,END,FAIL)",
    )
    ap.add_argument(
        "--extra",
        action="append",
        default=[],
        help='Extra BLAST arguments. Repeatable. Example: --extra "-qcov_hsp_perc 30" --extra "-word_size 2"',
    )
    ap.add_argument(
        "--submit",
        action="store_true",
        help="If set, submit each generated sbatch file via sbatch",
    )

    return ap


def main() -> int:
    ap = build_parser()
    args = ap.parse_args()

    inputdir: Path = args.inputdir
    if not inputdir.exists() or not inputdir.is_dir():
        print(f"ERROR: inputdir not found or not a directory: {inputdir}", file=sys.stderr)
        return 2

    threads = int(args.threads)
    if threads < 1:
        print("ERROR: --threads must be >= 1", file=sys.stderr)
        return 2

    mem_per_cpu_gb = int(args.mem_per_cpu_gb)
    if mem_per_cpu_gb < 1:
        print("ERROR: --mem-per-cpu-gb must be >= 1", file=sys.stderr)
        return 2

    taxid_val = int(args.taxid)
    if taxid_val == 0:
        taxid: int | None = None
    elif taxid_val < 0:
        print("ERROR: --taxid must be positive, or 0 to disable filtering", file=sys.stderr)
        return 2
    else:
        taxid = taxid_val

    outdir: Path = args.outdir if args.outdir is not None else (inputdir / "blastp_out")
    logs_dir = outdir / args.logs_dirname
    sbatch_dir = outdir / args.sbatch_dirname
    tsv_dir = outdir / "tsv"

    sbatch_dir.mkdir(parents=True, exist_ok=True)
    tsv_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    extra_blast_args = parse_extra_blast_args(args.extra)

    fasta_files = sorted([p for p in inputdir.iterdir() if is_fasta_file(p)])
    if not fasta_files:
        print(f"ERROR: no FASTA files found in {inputdir}", file=sys.stderr)
        print(f"Recognized suffixes: {', '.join(sorted(FASTA_SUFFIXES))}", file=sys.stderr)
        return 2

    mem_gb = format_mem_gb(threads, mem_per_cpu_gb)

    created: list[Path] = []
    for q in fasta_files:
        base = q.stem
        job_name = sanitize_job_name(base)

        if taxid is None:
            out_tsv = tsv_dir / f"{base}_vs_nr.tsv"
        else:
            out_tsv = tsv_dir / f"{base}_vs_taxid{taxid}.tsv"

        sbatch_path = sbatch_dir / f"{base}.sbatch"

        write_sbatch(
            sbatch_path,
            job_name=job_name,
            partition=args.partition,
            time_limit=args.time_limit,
            threads=threads,
            mem_gb=mem_gb,
            mail_user=args.mail_user,
            mail_type=args.mail_type,
            query_fasta=q,
            out_tsv=out_tsv,
            db=args.db,
            taxid=taxid,
            evalue=args.evalue,
            max_target_seqs=args.max_target_seqs,
            outfmt=args.outfmt,
            extra_blast_args=extra_blast_args,
            logs_dir=logs_dir,
        )
        created.append(sbatch_path)

    print(f"Created {len(created)} sbatch files in: {sbatch_dir}")
    if args.submit:
        for sb in created:
            print(f"Submitting: {sb.name}")
            subprocess.run(["sbatch", str(sb)], check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
