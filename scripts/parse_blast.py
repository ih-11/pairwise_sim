import re
from pathlib import Path


def parse_blast_output(blast_file):
    """
    Parse key metrics from BLAST pairwise output (outfmt 0).

    Parameters
    ----------
    blast_file : str or Path
        Path to BLAST output file.

    Returns
    -------
    dict
        Dictionary with parsed metrics.
    """

    blast_file = Path(blast_file)

    if not blast_file.exists():
        raise FileNotFoundError(f"BLAST output not found: {blast_file}")

    text = blast_file.read_text()

    metrics = {
        "bitscore": None,
        "evalue": None,
        "align_len": None,
        "pct_ident": None,
        "gaps": None,
    }

    # ---- Score and evalue ----
    m = re.search(r"Score =\s+([\d\.]+) bits.*Expect = ([\deE\.\-]+)", text)
    if m:
        metrics["bitscore"] = float(m.group(1))
        metrics["evalue"] = m.group(2)

    # ---- Identities (alignment length + % identity) ----
    m = re.search(r"Identities = (\d+)/(\d+) \((\d+)%\)", text)
    if m:
        metrics["align_len"] = int(m.group(2))
        metrics["pct_ident"] = float(m.group(3))

    # ---- Gaps ----
    m = re.search(r"Gaps = (\d+)/(\d+)", text)
    if m:
        metrics["gaps"] = int(m.group(1))

    return metrics