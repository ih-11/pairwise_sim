"""Phase-1 CLI placeholder for gradual-error simulation workflow."""

from __future__ import annotations

import argparse
from pathlib import Path


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
    parser.add_argument("--del", dest="delete", required=True, type=float, help="Deletion fraction.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--outdir", default="results/sim3", help="Base output directory.")
    parser.add_argument("--blast-task", default="blastn", help="BLAST task/program setting.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    print("[pairwise_sim] sim3 placeholder")
    print(f"input      : {Path(args.input)}")
    print(f"error      : {args.min_error}% to {args.max_error}% step {args.step}")
    print(f"fractions  : sub={args.sub}, ins={args.ins}, del={args.delete}")
    print(f"seed       : {args.seed}")
    print(f"outdir     : {Path(args.outdir)}")
    print("Status     : repository skeleton ready; implementation starts in phase 4.")


if __name__ == "__main__":
    main()
