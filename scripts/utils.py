"""Utility helpers for pairwise_sim.

These functions are intentionally minimal in phase 1.
They will be expanded in later phases.
"""

from __future__ import annotations

from pathlib import Path
import re


def sanitize_seqid(seq_id: str) -> str:
    """Convert a sequence ID into a filesystem-safe name."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", seq_id.strip())
    return cleaned or "sequence"


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if it does not already exist."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
