from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable


class BlastError(RuntimeError):
    """Raised when a BLAST-related command fails."""


def _ensure_executable(name: str) -> None:
    """Check that an external executable is available on PATH."""
    if shutil.which(name) is None:
        raise BlastError(
            f"Required executable '{name}' was not found on PATH. "
            "Make sure the conda environment is activated and BLAST+ is installed."
        )


def _run_command(cmd: Iterable[str], cwd: Path | None = None) -> None:
    """Run a shell command and raise a readable error if it fails."""
    try:
        subprocess.run(
            list(cmd),
            check=True,
            cwd=str(cwd) if cwd else None,
            text=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else "No stderr output."
        stdout = exc.stdout.strip() if exc.stdout else "No stdout output."
        raise BlastError(
            "Command failed.\n"
            f"CMD: {' '.join(exc.cmd)}\n"
            f"STDOUT:\n{stdout}\n\n"
            f"STDERR:\n{stderr}"
        ) from exc


def make_blast_db(fasta_path: str | Path, db_dir: str | Path, db_name: str = "db") -> Path:
    """
    Create a nucleotide BLAST database from a FASTA file.

    Parameters
    ----------
    fasta_path : str or Path
        Path to the input FASTA file.
    db_dir : str or Path
        Directory where BLAST database files will be written.
    db_name : str
        Base name for the BLAST database.

    Returns
    -------
    Path
        Path prefix of the created BLAST database.
        Example: results/sim1/test/db/db
    """
    _ensure_executable("makeblastdb")

    fasta_path = Path(fasta_path)
    db_dir = Path(db_dir)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / db_name

    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ]
    _run_command(cmd)

    return db_prefix


def run_blastn(
    query_fasta: str | Path,
    db_prefix: str | Path,
    out_file: str | Path,
    *,
    outfmt: int = 0,
    task: str = "blastn",
    word_size: int | None = None,
    reward: int | None = None,
    penalty: int | None = None,
    gapopen: int | None = None,
    gapextend: int | None = None,
) -> Path:
    """
    Run blastn against a nucleotide BLAST database.
    """
    _ensure_executable("blastn")

    query_fasta = Path(query_fasta)
    db_prefix = Path(db_prefix)
    out_file = Path(out_file)

    if not query_fasta.exists():
        raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

    out_file.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "blastn",
        "-task", task,
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-out", str(out_file),
        "-outfmt", str(outfmt),
    ]

    if word_size is not None:
        cmd.extend(["-word_size", str(word_size)])
    if reward is not None:
        cmd.extend(["-reward", str(reward)])
    if penalty is not None:
        cmd.extend(["-penalty", str(penalty)])
    if gapopen is not None:
        cmd.extend(["-gapopen", str(gapopen)])
    if gapextend is not None:
        cmd.extend(["-gapextend", str(gapextend)])

    _run_command(cmd)

    return out_file