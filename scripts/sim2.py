from __future__ import annotations

import argparse
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
        description="Mutate a nucleotide FASTA sequence, then align it to the original with BLAST+."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--error", required=True, type=float)

    parser.add_argument("--sub", required=True, type=float)
    parser.add_argument("--ins", required=True, type=float)
    parser.add_argument("--del_", dest="del_frac", required=True, type=float)

    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default="results/sim2")

    parser.add_argument("--blast-task", default="blastn")
    parser.add_argument("--word-size", type=int, default=None)
    parser.add_argument("--reward", type=int, default=None)
    parser.add_argument("--penalty", type=int, default=None)
    parser.add_argument("--gapopen", type=int, default=None)
    parser.add_argument("--gapextend", type=int, default=None)

    parser.add_argument("--keep-db", action="store_true")
    parser.add_argument("--force", action="store_true")

    return parser


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    seq_id = input_path.stem
    run_name = f"e{int(args.error):03d}_s{args.seed}"
    run_dir = Path(args.outdir) / seq_id / run_name

    db_dir = run_dir / "db"
    aln_file = run_dir / "aln.txt"
    metrics_file = run_dir / "metrics.tsv"
    mut_fasta = run_dir / "query_mut.fa"
    mut_log = run_dir / "mut_log.tsv"

    if run_dir.exists() and any(run_dir.iterdir()) and not args.force:
        raise FileExistsError(f"{run_dir} exists. Use --force to overwrite.")

    run_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] sim2 started: {seq_id}")
    print(f"[INFO] error={args.error}, seed={args.seed}")

    # ---- read original ----
    header, seq = read_single_fasta(input_path)

    # ---- mutate ----
    mut_seq, mut_log_rows, counts = mutate_sequence(
        sequence=seq,
        error_rate=args.error,
        sub_frac=args.sub,
        ins_frac=args.ins,
        del_frac=args.del_frac,
        seed=args.seed,
    )

    print(f"[INFO] edit counts: {counts}")

    write_fasta(f"{seq_id}_mut", mut_seq, mut_fasta)
    write_mutation_log(mut_log_rows, mut_log)

    # ---- BLAST ----
    db_prefix = make_blast_db(input_path, db_dir)

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

    # ---- save metrics ----
    with open(metrics_file, "w") as f:
        header = [
            "seq_id",
            "error_rate",
            "seed",
            "len_ref",
            "len_query",
            "n_sub",
            "n_ins",
            "n_del",
        ] + list(metrics.keys())

        values = [
            seq_id,
            args.error,
            args.seed,
            len(seq),
            len(mut_seq),
            counts["sub"],
            counts["ins"],
            counts["del"],
        ] + list(metrics.values())

        f.write("\t".join(map(str, header)) + "\n")
        f.write("\t".join(map(str, values)) + "\n")

    if not args.keep_db:
        shutil.rmtree(db_dir)

    print(f"[DONE] results in {run_dir}")


if __name__ == "__main__":
    main()