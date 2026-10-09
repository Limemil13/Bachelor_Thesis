#!/usr/bin/env python3
"""Prepare full and SignalP-trimmed canonical protein sets for motif discovery."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
EXPECTED_CANONICAL_COUNT = 103


def read_fasta(path: Path) -> dict[str, str]:
    # Read canonical proteins keyed by their compact project identifiers.
    records: dict[str, list[str]] = {}
    header: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            header = line[1:].split()[0]
            records[header] = []
        elif header is not None:
            records[header].append(line.strip())
    return {key: "".join(value) for key, value in records.items()}


def write_fasta(path: Path, records: dict[str, str]) -> None:
    # Write deterministic FASTA files in dictionary insertion order.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records.items():
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def main() -> None:
    # Validate all 103 canonical proteins, then write full-length and
    # SignalP-cleavage-trimmed sets for global and per-gene motif searches.
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--protein-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    with args.protein_table.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_accession = {row["Protein accession"]: row for row in rows}

    mature: dict[str, str] = {}
    metadata: list[dict[str, str | int]] = []
    for header, sequence in records.items():
        gene, species, accession = header.split("|", 2)
        row = by_accession[accession]
        removed = 0
        cleavage = row["SignalP cleavage site"]
        if row["SignalP prediction"] == "SP" and cleavage:
            match = re.match(r"(\d+)-(\d+)", cleavage)
            if match:
                removed = int(match.group(1))
        mature[header] = sequence[removed:]
        metadata.append(
            {
                "gene": gene,
                "species": species,
                "accession": accession,
                "full_length_aa": len(sequence),
                "signal_peptide_removed_aa": removed,
                "mature_length_aa": len(sequence) - removed,
                "signalp_prediction": row["SignalP prediction"],
                "signalp_quality": row["SignalP quality"],
            }
        )

    if len(mature) != EXPECTED_CANONICAL_COUNT:
        raise SystemExit(
            f"Expected {EXPECTED_CANONICAL_COUNT} canonical proteins, found {len(mature)}"
        )
    write_fasta(args.output_dir / "slrp_canonical_103_mature.faa", mature)
    write_fasta(args.output_dir / "slrp_canonical_103_full.faa", records)

    for gene in GENES:
        positives = {h: s for h, s in mature.items() if h.startswith(f"{gene}|")}
        controls = {h: s for h, s in mature.items() if not h.startswith(f"{gene}|")}
        write_fasta(args.output_dir / "per_gene" / f"{gene}_mature.faa", positives)
        write_fasta(
            args.output_dir / "per_gene" / f"{gene}_other_genes_control.faa", controls
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(metadata[0])
    with (args.output_dir / "protein_motif_input_metadata.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata)
    print(
        "Prepared 103 full and mature canonical proteins for seven-gene MEME analysis"
    )


if __name__ == "__main__":
    main()
