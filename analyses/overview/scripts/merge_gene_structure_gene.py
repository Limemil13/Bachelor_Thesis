#!/usr/bin/env python3
"""Replace one gene's rows in a combined gene-structure extraction table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing", required=True, type=Path)
    parser.add_argument("--gene-table", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    existing = read_tsv(args.existing)
    added = read_tsv(args.gene_table)
    if not added or {row["query_gene"] for row in added} != {args.gene}:
        raise ValueError(f"Replacement table must contain only {args.gene}")
    expected_species = {"human", "mouse", "cow", "chicken", "zebrafish"}
    if {row["species_or_file"] for row in added} != expected_species:
        raise ValueError(f"{args.gene} table does not cover all five reference species")

    retained = [row for row in existing if row["query_gene"] != args.gene]
    rows = retained + added
    fields = list(existing[0])
    if set(fields) != set(added[0]):
        raise ValueError("Combined and replacement gene-structure columns differ")
    rows.sort(
        key=lambda row: (
            row["species_or_file"],
            row["query_gene"],
            row["transcript_id"],
        )
    )
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Wrote {args.output} ({len(rows)} transcript rows; {len(added)} {args.gene})"
    )


if __name__ == "__main__":
    main()
