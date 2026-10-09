#!/usr/bin/env python3
"""Extract a 1-based inclusive interval from a FASTA record."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    # Stream one named FASTA record, validate 1-based inclusive coordinates and
    # write only the requested subsequence.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("record")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    if args.start < 1 or args.end < args.start:
        raise ValueError("Coordinates must satisfy 1 <= start <= end")

    # Avoid loading unrelated records from a potentially large assembly FASTA.
    collecting = False
    parts: list[str] = []
    with args.fasta.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                current = raw_line[1:].split(maxsplit=1)[0]
                if collecting:
                    break
                collecting = current == args.record
                continue
            if collecting:
                parts.append(raw_line.strip())

    if not parts:
        raise ValueError(f"FASTA record not found: {args.record}")
    sequence = "".join(parts)
    if args.end > len(sequence):
        raise ValueError(
            f"Requested end {args.end} exceeds record length {len(sequence)}"
        )
    subsequence = sequence[args.start - 1 : args.end]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(f">{args.record}:{args.start}-{args.end}\n")
        for offset in range(0, len(subsequence), 80):
            handle.write(subsequence[offset : offset + 80] + "\n")

    print(f"Wrote {len(subsequence)} bp to {args.output}")


if __name__ == "__main__":
    main()
