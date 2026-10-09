#!/usr/bin/env python3
"""Create the canonical SLRP tree input without two partial BGN records."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXCLUDED_ACCESSIONS = {
    "XP_038638277.1": "catshark partial BGN model; 68% query coverage",
    "XP_048475905.1": "whale-shark partial BGN model; 39% query coverage",
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            elif header is None:
                raise ValueError(f"Sequence found before a FASTA header in {path}")
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def parse_header(header: str) -> tuple[str, str, str]:
    token = header.split()[0]
    fields = token.split("|")
    if len(fields) != 3:
        raise ValueError(f"Expected gene|species|accession header, found: {header}")
    return fields[0], fields[1], fields[2]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--excluded-table", required=True, type=Path)
    parser.add_argument("--count-table", required=True, type=Path)
    args = parser.parse_args()

    records = read_fasta(args.input)
    if len(records) != 103:
        raise ValueError(f"Expected 103 canonical records, found {len(records)}")

    kept: list[tuple[str, str]] = []
    excluded_rows: list[dict[str, object]] = []
    before: Counter[str] = Counter()
    after: Counter[str] = Counter()
    for header, sequence in records:
        gene, species, accession = parse_header(header)
        before[gene] += 1
        if accession in EXCLUDED_ACCESSIONS:
            if gene != "BGN":
                raise ValueError(f"Excluded accession is not labelled BGN: {header}")
            excluded_rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "accession": accession,
                    "protein_length_aa": len(sequence),
                    "sensitivity_exclusion_reason": EXCLUDED_ACCESSIONS[accession],
                }
            )
            continue
        kept.append((header, sequence))
        after[gene] += 1

    found = {str(row["accession"]) for row in excluded_rows}
    if found != set(EXCLUDED_ACCESSIONS):
        missing = sorted(set(EXCLUDED_ACCESSIONS) - found)
        raise ValueError(f"Expected partial BGN records were not found: {missing}")
    if len(kept) != 101:
        raise ValueError(f"Expected 101 retained records, found {len(kept)}")

    write_fasta(args.output, kept)
    args.excluded_table.parent.mkdir(parents=True, exist_ok=True)
    with args.excluded_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(excluded_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(excluded_rows)

    genes = sorted(before)
    with args.count_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("gene", "canonical_count", "sensitivity_count", "removed"),
            delimiter="\t",
        )
        writer.writeheader()
        for gene in genes:
            writer.writerow(
                {
                    "gene": gene,
                    "canonical_count": before[gene],
                    "sensitivity_count": after[gene],
                    "removed": before[gene] - after[gene],
                }
            )

    print(f"Retained {len(kept)}/{len(records)} proteins; excluded 2 partial BGN records")


if __name__ == "__main__":
    main()
