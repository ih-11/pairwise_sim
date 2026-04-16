"""Phase-1 CLI placeholder for one-error-level mutation + alignment workflow."""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Mutate a nucleotide FASTA sequence, then align it to the original with local BLAST+."
    )
    parser.add_argument("--input", required=True, help="Path to input FASTA file.")
    parser.add_argument("--error", required=True, type=float, help="Total error rate in percent.")
    parser.add_argument("--sub", required=True, type=float, help="Substitution fraction.")
    parser.add_argument("--ins", required=True, type=float, help="Insertion fraction.")
    parser.add_argument("--del", dest="delete", required=True, type=float, help="Deletion fraction.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--outdir", default="results/sim2", help="Base output directory.")
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
    print("[pairwise_sim] sim2 placeholder")
    print(f"input      : {Path(args.input)}")
    print(f"error      : {args.error}%")
    print(f"fractions  : sub={args.sub}, ins={args.ins}, del={args.delete}")
    print(f"seed       : {args.seed}")
    print(f"outdir     : {Path(args.outdir)}")
    print("Status     : repository skeleton ready; implementation starts in phase 3.")


if __name__ == "__main__":
    main()
