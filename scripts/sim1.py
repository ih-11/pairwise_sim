"""Phase-1 CLI placeholder for self-alignment workflow."""

from __future__ import annotations

import argparse
from pathlib import Path


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
    print("[pairwise_sim] sim1 placeholder")
    print(f"input      : {Path(args.input)}")
    print(f"outdir     : {Path(args.outdir)}")
    print(f"blast_task : {args.blast_task}")
    print("Status     : repository skeleton ready; implementation starts in phase 2.")


if __name__ == "__main__":
    main()
