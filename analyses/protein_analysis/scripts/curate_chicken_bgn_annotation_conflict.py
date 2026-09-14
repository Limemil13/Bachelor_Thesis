#!/usr/bin/env python3
"""Exclude the ASPN-like chicken BGN-locus protein from the BGN protein panel."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CANDIDATES = ROOT / "analyses" / "protein_analysis" / "candidates"
SOURCE = CANDIDATES / "BGN_domain_input.faa"
OUTPUT = CANDIDATES / "BGN_domain_input_chicken_asporin_excluded.faa"
REPORT = CANDIDATES / "BGN_chicken_annotation_conflict_curation.tsv"
TARGET = "XP_414298.2"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    chunks = []
    for raw in path.read_text(encoding="utf-8").splitlines():
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
    excluded = [record for record in records if TARGET in record[0]]
    retained = [record for record in records if TARGET not in record[0]]
    if len(records) != 16 or len(excluded) != 1 or len(retained) != 15:
        raise SystemExit(
            f"Expected 16 source / 1 excluded / 15 retained; got "
            f"{len(records)} / {len(excluded)} / {len(retained)}"
        )

    with OUTPUT.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in retained:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    row = {
        "gene": "BGN",
        "species": "chicken",
        "accession": TARGET,
        "action": "excluded_from_canonical_BGN_protein_panel",
        "ncbi_gene_locus": "BGN (GeneID 415954; MODEL)",
        "ncbi_product": "asporin",
        "reciprocal_human_hit": "ASPN NP_060150.4",
        "reference_tree_result": "sister_to_human_ASPN",
        "reference_tree_support": "99.3_SH-aLRT/100_UFBoot",
        "synvoy_result": "no_retained_HIGH_or_MEDIUM_chicken_BGN_candidate",
        "interpretation": (
            "Retain as annotation-conflict provenance; do not claim it as a "
            "chicken BGN protein ortholog without independent locus/synteny resolution."
        ),
    }
    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(f"Wrote {OUTPUT} ({len(retained)} retained BGN proteins)")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
