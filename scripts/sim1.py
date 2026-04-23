from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

# Allow running as: python3 scripts/sim1.py ...
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.blast_runner import make_blast_db, run_blastn
from scripts.parse_blast import parse_blast_output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run self-alignment of a nucleotide FASTA sequence with local BLAST+."
    )
    parser.add_argument("--input", required=True, help="Path to input FASTA file.")
    parser.add_argument("--outdir", default="results/sim1", help="Base output directory.")
    parser.add_argument("--blast-task", default="blastn", help="BLAST task/program setting.")
    parser.add_argument("--word-size", type=int, default=None, help="Optional BLAST word size.")
    parser.add_argument("--reward", type=int, default=None, help="Optional BLAST reward.")
    parser.add_argument("--penalty", type=int, default=None, help="Optional BLAST penalty.")
    parser.add_argument("--gapopen", type=int, default=None, help="Optional gap open penalty.")
    parser.add_argument("--gapextend", type=int, default=None, help="Optional gap extend penalty.")
    parser.add_argument("--keep-db", action="store_true", help="Keep BLAST database files.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    seq_id = input_path.stem
    run_dir = Path(args.outdir) / seq_id
    db_dir = run_dir / "db"
    aln_file = run_dir / "aln.txt"
    metrics_file = run_dir / "metrics.tsv"

    if run_dir.exists() and any(run_dir.iterdir()) and not args.force:
        raise FileExistsError(
            f"Output directory already exists and is not empty: {run_dir}\n"
            "Use --force to overwrite."
        )

    run_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] sim1 started for: {seq_id}")
    print(f"[INFO] Input FASTA: {input_path}")
    print(f"[INFO] Output dir : {run_dir}")

    db_prefix = make_blast_db(input_path, db_dir)

    run_blastn(
        query_fasta=input_path,
        db_prefix=db_prefix,
        out_file=aln_file,
        task=args.blast_task,
        word_size=args.word_size,
        reward=args.reward,
        penalty=args.penalty,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
    )

    metrics = parse_blast_output(aln_file)

    with open(metrics_file, "w", encoding="utf-8") as handle:
        header = ["seq_id"] + list(metrics.keys())
        values = [seq_id] + [metrics[key] for key in metrics.keys()]
        handle.write("\t".join(map(str, header)) + "\n")
        handle.write("\t".join(map(str, values)) + "\n")

    if not args.keep_db and db_dir.exists():
        shutil.rmtree(db_dir)

    print(f"[DONE] Alignment report : {aln_file}")
    print(f"[DONE] Metrics table    : {metrics_file}")


if __name__ == "__main__":
    main()