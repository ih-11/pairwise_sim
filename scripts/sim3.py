from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

# allow running as script
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.blast_runner import make_blast_db, run_blastn
from scripts.parse_blast import parse_blast_output
from scripts.mutate_sequence import (
    read_single_fasta,
    mutate_sequence,
    write_fasta,
    write_mutation_log,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run mutation + alignment across a gradual error range, then aggregate metrics."
    )
    parser.add_argument("--input", required=True, help="Path to input FASTA file.")
    parser.add_argument("--min-error", required=True, type=float, help="Minimum error rate in percent.")
    parser.add_argument("--max-error", required=True, type=float, help="Maximum error rate in percent.")
    parser.add_argument("--step", required=True, type=float, help="Error-rate increment.")

    parser.add_argument("--sub", required=True, type=float, help="Substitution fraction.")
    parser.add_argument("--ins", required=True, type=float, help="Insertion fraction.")
    parser.add_argument("--del_", dest="del_frac", required=True, type=float, help="Deletion fraction.")

    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--outdir", default="results/sim3", help="Base output directory.")

    parser.add_argument("--blast-task", default="blastn", help="BLAST task/program setting.")
    parser.add_argument("--word-size", type=int, default=None, help="Optional BLAST word size.")
    parser.add_argument("--reward", type=int, default=None, help="Optional BLAST reward.")
    parser.add_argument("--penalty", type=int, default=None, help="Optional BLAST penalty.")
    parser.add_argument("--gapopen", type=int, default=None, help="Optional gap open penalty.")
    parser.add_argument("--gapextend", type=int, default=None, help="Optional gap extend penalty.")

    parser.add_argument("--keep-db", action="store_true", help="Keep BLAST database files.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser


def frange(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("step must be > 0")
    values = []
    x = start
    while x <= stop + 1e-9:
        values.append(round(x, 10))
        x += step
    return values


def error_to_tag(error_rate: float) -> str:
    return f"e{int(round(error_rate)):03d}"


def main() -> None:
    args = build_parser().parse_args()

    if args.min_error > args.max_error:
        raise ValueError("--min-error must be <= --max-error")

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    seq_id = input_path.stem
    run_root = Path(args.outdir) / seq_id / f"seed{args.seed}"
    combined_metrics_file = run_root / "metrics_all.tsv"

    if run_root.exists() and any(run_root.iterdir()) and not args.force:
        raise FileExistsError(f"{run_root} exists. Use --force to overwrite.")

    run_root.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] sim3 started: {seq_id}")
    print(f"[INFO] error range: {args.min_error} to {args.max_error} (step {args.step})")
    print(f"[INFO] seed: {args.seed}")

    header, seq = read_single_fasta(input_path)

    db_dir = run_root / "db"
    db_prefix = make_blast_db(input_path, db_dir)

    rows: list[dict[str, object]] = []

    for error_rate in frange(args.min_error, args.max_error, args.step):
        tag = error_to_tag(error_rate)
        run_dir = run_root / tag
        run_dir.mkdir(parents=True, exist_ok=True)

        print(f"[INFO] running {tag} ({error_rate}%)")

        mut_fasta = run_dir / "query_mut.fa"
        mut_log = run_dir / "mut_log.tsv"
        aln_file = run_dir / "aln.txt"

        mut_seq, mut_log_rows, counts = mutate_sequence(
            sequence=seq,
            error_rate=error_rate,
            sub_frac=args.sub,
            ins_frac=args.ins,
            del_frac=args.del_frac,
            seed=args.seed,
        )

        write_fasta(f"{seq_id}_mut_{tag}", mut_seq, mut_fasta)
        write_mutation_log(mut_log_rows, mut_log)

        run_blastn(
            query_fasta=mut_fasta,
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

        row = {
            "seq_id": seq_id,
            "error_rate": error_rate,
            "seed": args.seed,
            "len_ref": len(seq),
            "len_query": len(mut_seq),
            "n_sub": counts["sub"],
            "n_ins": counts["ins"],
            "n_del": counts["del"],
            "bitscore": metrics["bitscore"],
            "evalue": metrics["evalue"],
            "align_len": metrics["align_len"],
            "pct_ident": metrics["pct_ident"],
            "gaps": metrics["gaps"],
            "run_dir": str(run_dir),
        }
        rows.append(row)

    fieldnames = [
        "seq_id",
        "error_rate",
        "seed",
        "len_ref",
        "len_query",
        "n_sub",
        "n_ins",
        "n_del",
        "bitscore",
        "evalue",
        "align_len",
        "pct_ident",
        "gaps",
        "run_dir",
    ]

    with open(combined_metrics_file, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    if not args.keep_db and db_dir.exists():
        shutil.rmtree(db_dir)

    print(f"[DONE] combined metrics: {combined_metrics_file}")
    print(f"[DONE] run root: {run_root}")


if __name__ == "__main__":
    main()