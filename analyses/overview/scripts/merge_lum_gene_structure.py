#!/usr/bin/env python3
"""Merge a five-species LUM extraction into the established gene-structure table."""

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
    parser.add_argument("--lum", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    existing = [row for row in read_tsv(args.existing) if row["query_gene"] != "LUM"]
    lum = read_tsv(args.lum)
    context_genes = {"BGN", "EPYC", "FMOD", "OGN", "OMD", "PRELP"}
    context_species = {"human", "mouse", "cow", "chicken", "zebrafish"}
    context_pairs = {(row["query_gene"], row["species_or_file"]) for row in existing}
    if {row["query_gene"] for row in existing} != context_genes or context_pairs != {
        (gene, species) for gene in context_genes for species in context_species
    }:
        raise ValueError("Expected complete six-gene x five-species context coverage")
    if not lum or {row["query_gene"] for row in lum} != {"LUM"}:
        raise ValueError("LUM extraction contains non-LUM or no rows")
    if {row["species_or_file"] for row in lum} != {
        "human",
        "mouse",
        "cow",
        "chicken",
        "zebrafish",
    }:
        raise ValueError("LUM extraction does not cover all five reference species")

    rows = existing + lum
    fields = list(rows[0])
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Wrote {args.output} ({len(rows)} transcript rows before representative selection)"
    )


if __name__ == "__main__":
    main()
