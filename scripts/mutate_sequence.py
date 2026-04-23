from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List, Tuple

DNA_ALPHABET = ("A", "C", "G", "T")


def read_single_fasta(fasta_path: str | Path) -> Tuple[str, str]:
    """
    Read a single-sequence FASTA file.

    Returns
    -------
    tuple[str, str]
        (header_without_>, sequence)
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    header = None
    seq_parts: List[str] = []

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError("FASTA contains more than one sequence; expected exactly one.")
                header = line[1:]
            else:
                seq_parts.append(line.upper())

    if header is None:
        raise ValueError("Invalid FASTA: missing header line starting with '>'.")

    sequence = "".join(seq_parts)
    if not sequence:
        raise ValueError("FASTA sequence is empty.")

    invalid = set(sequence) - set(DNA_ALPHABET)
    if invalid:
        raise ValueError(
            f"Sequence contains unsupported characters: {sorted(invalid)}. "
            "Version 1 supports only A/C/G/T."
        )

    return header, sequence


def allocate_edit_counts(
    seq_len: int,
    error_rate: float,
    sub_frac: float,
    ins_frac: float,
    del_frac: float,
) -> Dict[str, int]:
    """
    Convert error rate and edit fractions into integer counts.

    error_rate is percent of the original sequence length.
    """
    if seq_len <= 0:
        raise ValueError("Sequence length must be > 0.")
    if error_rate < 0:
        raise ValueError("error_rate must be >= 0.")

    total_frac = sub_frac + ins_frac + del_frac
    if abs(total_frac - 1.0) > 1e-9:
        raise ValueError("sub_frac + ins_frac + del_frac must sum to 1.0.")

    total_edits = round(seq_len * (error_rate / 100.0))

    raw_sub = total_edits * sub_frac
    raw_ins = total_edits * ins_frac
    raw_del = total_edits * del_frac

    sub_n = int(round(raw_sub))
    ins_n = int(round(raw_ins))
    del_n = int(round(raw_del))

    diff = total_edits - (sub_n + ins_n + del_n)
    if diff != 0:
        buckets = [
            ("sub", raw_sub - int(raw_sub)),
            ("ins", raw_ins - int(raw_ins)),
            ("del", raw_del - int(raw_del)),
        ]
        if diff > 0:
            buckets.sort(key=lambda x: x[1], reverse=True)
            for i in range(diff):
                if buckets[i % 3][0] == "sub":
                    sub_n += 1
                elif buckets[i % 3][0] == "ins":
                    ins_n += 1
                else:
                    del_n += 1
        else:
            buckets.sort(key=lambda x: x[1])
            for i in range(abs(diff)):
                if buckets[i % 3][0] == "sub" and sub_n > 0:
                    sub_n -= 1
                elif buckets[i % 3][0] == "ins" and ins_n > 0:
                    ins_n -= 1
                elif buckets[i % 3][0] == "del" and del_n > 0:
                    del_n -= 1

    if del_n > seq_len:
        raise ValueError(
            f"Deletion count ({del_n}) exceeds sequence length ({seq_len}). "
            "Reduce error rate or deletion fraction."
        )

    return {
        "total": total_edits,
        "sub": sub_n,
        "ins": ins_n,
        "del": del_n,
    }


def _random_alt_base(ref_base: str, rng: random.Random) -> str:
    choices = [b for b in DNA_ALPHABET if b != ref_base]
    return rng.choice(choices)


def mutate_sequence(
    sequence: str,
    error_rate: float,
    sub_frac: float,
    ins_frac: float,
    del_frac: float,
    seed: int,
) -> Tuple[str, List[Dict[str, object]], Dict[str, int]]:
    """
    Apply single-nucleotide substitutions, insertions, and deletions.

    Returns
    -------
    tuple
        mutated_sequence, mutation_log, edit_counts
    """
    rng = random.Random(seed)
    seq_len = len(sequence)
    counts = allocate_edit_counts(seq_len, error_rate, sub_frac, ins_frac, del_frac)

    seq_list = list(sequence)

    deletion_positions = set(rng.sample(range(seq_len), counts["del"])) if counts["del"] > 0 else set()

    remaining_positions = [i for i in range(seq_len) if i not in deletion_positions]
    if counts["sub"] > len(remaining_positions):
        raise ValueError(
            "Substitution count exceeds number of non-deleted positions. "
            "Reduce error rate or change proportions."
        )

    substitution_positions = (
        set(rng.sample(remaining_positions, counts["sub"])) if counts["sub"] > 0 else set()
    )

    insertion_slots = [rng.randrange(seq_len + 1) for _ in range(counts["ins"])]
    insertion_events_by_slot: Dict[int, List[str]] = {}
    for slot in insertion_slots:
        insertion_events_by_slot.setdefault(slot, []).append(rng.choice(DNA_ALPHABET))

    mutation_log: List[Dict[str, object]] = []
    mutated_chars: List[str] = []
    event_id = 1

    for slot in range(seq_len + 1):
        if slot in insertion_events_by_slot:
            for ins_base in insertion_events_by_slot[slot]:
                mutation_log.append(
                    {
                        "event_id": event_id,
                        "event_type": "ins",
                        "position_ref": slot,
                        "ref_base": "-",
                        "alt_base": ins_base,
                    }
                )
                mutated_chars.append(ins_base)
                event_id += 1

        if slot == seq_len:
            continue

        ref_base = seq_list[slot]

        if slot in deletion_positions:
            mutation_log.append(
                {
                    "event_id": event_id,
                    "event_type": "del",
                    "position_ref": slot + 1,
                    "ref_base": ref_base,
                    "alt_base": "-",
                }
            )
            event_id += 1
            continue

        if slot in substitution_positions:
            alt_base = _random_alt_base(ref_base, rng)
            mutation_log.append(
                {
                    "event_id": event_id,
                    "event_type": "sub",
                    "position_ref": slot + 1,
                    "ref_base": ref_base,
                    "alt_base": alt_base,
                }
            )
            mutated_chars.append(alt_base)
            event_id += 1
        else:
            mutated_chars.append(ref_base)

    mutated_sequence = "".join(mutated_chars)
    return mutated_sequence, mutation_log, counts


def write_fasta(header: str, sequence: str, out_path: str | Path, width: int = 80) -> Path:
    """Write a FASTA file with fixed line width."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i:i + width] + "\n")

    return out_path


def write_mutation_log(log_rows: List[Dict[str, object]], out_path: str | Path) -> Path:
    """Write mutation log as TSV."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["event_id", "event_type", "position_ref", "ref_base", "alt_base"]

    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in log_rows:
            handle.write("\t".join(str(row[col]) for col in columns) + "\n")

    return out_path