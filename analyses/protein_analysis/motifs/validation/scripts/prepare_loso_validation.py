#!/usr/bin/env python3
"""Prepare leave-one-species-out inputs for protein motif validation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
EXCLUDED = {
    "XP_075929720.1": (
        "SLRP3-like",
        "broader reference tree reassigned the former lamprey BGN candidate",
    ),
    "XP_075930353.1": (
        "SLRP6-like",
        "broader reference tree reassigned the former lamprey LUM candidate",
    ),
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = line[1:].split()[0]
            sequence = []
        elif header is not None:
            sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def split_header(header: str) -> tuple[str, str, str]:
    parts = header.split("|", 2)
    if len(parts) != 3 or parts[0] not in GENES:
        raise ValueError(f"Unexpected FASTA identifier: {header}")
    return parts[0], parts[1], parts[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mature-fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    records = read_fasta(args.mature_fasta)
    if len(records) != 103:
        raise SystemExit(f"Expected 103 mature-protein records, found {len(records)}")

    parsed = []
    for header, sequence in records:
        gene, species, accession = split_header(header)
        parsed.append(
            {
                "header": header,
                "sequence": sequence,
                "gene": gene,
                "species": species,
                "accession": accession,
            }
        )

    excluded_rows: list[dict[str, object]] = []
    core = []
    challenges = []
    for row in parsed:
        accession = str(row["accession"])
        if accession in EXCLUDED:
            revised_label, reason = EXCLUDED[accession]
            excluded_rows.append(
                {
                    "original_gene": row["gene"],
                    "species": row["species"],
                    "accession": accession,
                    "revised_label": revised_label,
                    "reason": reason,
                }
            )
            challenges.append((str(row["header"]), str(row["sequence"])))
        else:
            core.append(row)

    if len(core) != 101 or len(excluded_rows) != 2:
        raise SystemExit(
            f"Expected 101 validation records and two excluded lamprey controls; "
            f"found {len(core)} and {len(excluded_rows)}"
        )

    species = sorted({str(row["species"]) for row in core})
    manifest_rows: list[dict[str, object]] = []
    for held_out in species:
        fold_dir = args.output_dir / "folds" / held_out
        test = [row for row in core if row["species"] == held_out]
        train = [row for row in core if row["species"] != held_out]
        write_fasta(
            fold_dir / "test.faa",
            [(str(row["header"]), str(row["sequence"])) for row in test],
        )
        for gene in GENES:
            gene_train = [row for row in train if row["gene"] == gene]
            if len(gene_train) < 2:
                raise SystemExit(f"Too few {gene} training records in fold {held_out}")
            write_fasta(
                fold_dir / "train" / f"{gene}.faa",
                [(str(row["header"]), str(row["sequence"])) for row in gene_train],
            )

        training_sequences = {str(row["sequence"]) for row in train}
        for row in test:
            manifest_rows.append(
                {
                    "fold": held_out,
                    "gene": row["gene"],
                    "species": row["species"],
                    "accession": row["accession"],
                    "sequence_length_aa": len(str(row["sequence"])),
                    "exact_sequence_present_in_training": (
                        "yes" if str(row["sequence"]) in training_sequences else "no"
                    ),
                    "training_record_count": len(train),
                    "test_record_count": len(test),
                }
            )

    for gene in GENES:
        gene_records = [row for row in core if row["gene"] == gene]
        write_fasta(
            args.output_dir / "final_models" / "train" / f"{gene}.faa",
            [(str(row["header"]), str(row["sequence"])) for row in gene_records],
        )
    write_fasta(args.output_dir / "challenge_controls.faa", challenges)

    write_tsv(
        args.output_dir / "validation_manifest.tsv",
        manifest_rows,
        list(manifest_rows[0]),
    )
    write_tsv(
        args.output_dir / "excluded_from_target_models.tsv",
        excluded_rows,
        list(excluded_rows[0]),
    )
    design_rows = [
        {"item": "input_records", "value": len(records)},
        {"item": "validated_target_records", "value": len(core)},
        {"item": "held_out_species_folds", "value": len(species)},
        {"item": "target_genes", "value": len(GENES)},
        {"item": "reclassified_challenge_controls", "value": len(challenges)},
    ]
    write_tsv(args.output_dir / "validation_design.tsv", design_rows, ["item", "value"])
    print(
        f"Prepared {len(species)} leave-one-species-out folds for {len(core)} proteins; "
        "two reclassified lamprey proteins were reserved as challenge controls"
    )


if __name__ == "__main__":
    main()
