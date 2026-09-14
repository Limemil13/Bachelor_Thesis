#!/usr/bin/env python3
"""Remove the duplicated lamprey LUM-like sequence from the FMOD panel.

XP_075930353.1 was historically selected for FMOD and was independently
selected for LUM.  It is identical in both panels, reverse-BLASTs to human
LUM, and is annotated lumican-like. The original FMOD FASTA remains intact;
the curated output includes a deduplicated FMOD panel and an audit record.
"""

from __future__ import annotations

import csv
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "candidates" / "FMOD_domain_input.faa"
OUTPUT = BASE / "candidates" / "FMOD_domain_input_lamprey_deduplicated.faa"
REPORT = BASE / "candidates" / "class_ii_lamprey_duplicate_curation.tsv"
ACCESSION = "XP_075930353.1"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def main() -> None:
    records = read_fasta(SOURCE)
    removed = [(header, seq) for header, seq in records if ACCESSION in header]
    kept = [(header, seq) for header, seq in records if ACCESSION not in header]
    if len(records) != 15 or len(removed) != 1 or len(kept) != 14:
        raise ValueError(
            "Expected to remove one XP_075930353.1 record from 15 FMOD records"
        )

    with OUTPUT.open("w", encoding="utf-8") as handle:
        for header, sequence in kept:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    fields = ["accession", "historical_panel", "final_panel", "decision", "evidence"]
    row = {
        "accession": ACCESSION,
        "historical_panel": "FMOD and LUM (duplicated)",
        "final_panel": "LUM",
        "decision": "exclude from canonical FMOD; retain original source as provenance",
        "evidence": "identical sequence; reverse hit NP_002336.1; NCBI lumican-like annotation",
    }
    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    print(f"Wrote {OUTPUT} ({len(kept)} FMOD proteins)")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
