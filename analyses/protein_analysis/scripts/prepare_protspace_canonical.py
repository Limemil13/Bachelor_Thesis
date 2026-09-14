"""Prepare the canonical SLRP + four-control ProtSpace input.

The previous ProtSpace input supplies only the four non-SLRP control
sequences and the stable ProtSpace identifiers.  Every SLRP sequence is
replaced from the canonical manifest/FASTA, including corrected dog OGN.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from canonical_dataset import COMBINED_FASTA, MANIFEST

BASE = Path(__file__).resolve().parents[1]
OUT_DIR = BASE / "protspace" / "input"
OUT_FASTA = OUT_DIR / "slrp_embedding_canonical_with_controls.faa"
OUT_ANNOTATIONS = OUT_DIR / "slrp_annotations_canonical_with_controls.csv"


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if header in records:
                        raise ValueError(
                            f"Duplicate FASTA identifier in {path}: {header}"
                        )
                    records[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        if header in records:
            raise ValueError(f"Duplicate FASTA identifier in {path}: {header}")
        records[header] = "".join(chunks)
    return records


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--previous-fasta", required=True, type=Path)
    parser.add_argument("--previous-annotations", required=True, type=Path)
    args = parser.parse_args()

    canonical_sequences = read_fasta(COMBINED_FASTA)
    previous_sequences = read_fasta(args.previous_fasta)
    manifest_rows = read_rows(MANIFEST, "\t")
    annotation_rows = read_rows(args.previous_annotations, ",")

    manifest = {row["canonical_id"]: row for row in manifest_rows}
    expected_slrp_count = len(manifest_rows)
    if (
        len(canonical_sequences) != expected_slrp_count
        or len(manifest) != expected_slrp_count
    ):
        raise ValueError("Canonical FASTA and manifest counts differ")

    previous_slrp = {
        row["canonical_id"]: row
        for row in annotation_rows
        if row["analysis_group"] == "SLRP_candidate" and row.get("canonical_id")
    }
    control_rows = [
        row for row in annotation_rows if row["analysis_group"] != "SLRP_candidate"
    ]
    if len(control_rows) != 4:
        raise ValueError(f"Expected four retained controls, found {len(control_rows)}")

    output_records: list[tuple[str, str]] = []
    output_annotations: list[dict[str, str]] = []
    seen_identifiers: set[str] = set()
    reserved_identifiers = {
        row["identifier"]
        for canonical_id, row in previous_slrp.items()
        if canonical_id in manifest
    } | {row["identifier"] for row in control_rows}
    slrp_count = 0
    control_count = 0

    class_by_gene = {
        "BGN": "SLRP_class_I",
        "DCN": "SLRP_class_I",
        "EPYC": "SLRP_class_III",
        "FMOD": "SLRP_class_II",
        "OGN": "SLRP_class_III",
        "PRELP": "SLRP_class_II",
        "LUM": "SLRP_class_II",
    }
    gene_indices: dict[str, int] = {}

    for record in manifest_rows:
        canonical_id = record["canonical_id"]
        gene = record["gene"]
        gene_indices[gene] = gene_indices.get(gene, 0) + 1
        prior = previous_slrp.get(canonical_id)
        if prior:
            identifier = prior["identifier"]
        else:
            candidate_index = 1
            identifier = f"{gene}_{candidate_index:03d}"
            while identifier in reserved_identifiers or identifier in seen_identifiers:
                candidate_index += 1
                identifier = f"{gene}_{candidate_index:03d}"
        if identifier in seen_identifiers:
            raise ValueError(f"Duplicate ProtSpace identifier: {identifier}")
        seen_identifiers.add(identifier)
        sequence = canonical_sequences[canonical_id]
        output_records.append((identifier, sequence))
        output_annotations.append(
            {
                "identifier": identifier,
                "gene": gene,
                "SLRP_class": class_by_gene[gene],
                "analysis_group": "SLRP_candidate",
                "canonical_id": canonical_id,
                "canonical_source": record["source_file"],
                "sequence_correction": record["sequence_correction"],
                "old_header": record["original_header"],
            }
        )
        slrp_count += 1

    for annotation in control_rows:
        identifier = annotation["identifier"]
        if identifier in seen_identifiers:
            raise ValueError(f"Duplicate ProtSpace identifier: {identifier}")
        if identifier not in previous_sequences:
            raise KeyError(f"Missing control sequence in previous FASTA: {identifier}")
        seen_identifiers.add(identifier)
        output_records.append((identifier, previous_sequences[identifier]))
        output_annotations.append(
            {
                "identifier": identifier,
                "gene": annotation["gene"],
                "SLRP_class": annotation["SLRP_class"],
                "analysis_group": annotation["analysis_group"],
                "canonical_id": "",
                "canonical_source": "retained_from_previous_control_set",
                "sequence_correction": "not_applicable",
                "old_header": annotation["old_header"],
            }
        )
        control_count += 1

    expected_total = expected_slrp_count + 4
    if (slrp_count, control_count, len(output_records)) != (
        expected_slrp_count,
        4,
        expected_total,
    ):
        raise ValueError(
            f"Unexpected ProtSpace composition: {slrp_count} SLRPs, "
            f"{control_count} controls, {len(output_records)} total"
        )

    dog_index = next(
        i
        for i, row in enumerate(output_annotations)
        if row["canonical_id"] == "OGN|canis_lupus_familiaris|XP_038383338.1"
    )
    if len(output_records[dog_index][1]) != 297:
        raise ValueError("Canonical ProtSpace dog OGN must be 297 aa")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_FASTA.open("w", encoding="utf-8") as handle:
        for identifier, sequence in output_records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    fields = list(output_annotations[0])
    with OUT_ANNOTATIONS.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_annotations)

    print(f"Wrote {OUT_FASTA}")
    print(f"Wrote {OUT_ANNOTATIONS}")
    print(
        f"Sequences: {len(output_records)} ({slrp_count} SLRPs + {control_count} controls)"
    )
    print(f"Dog OGN length: {len(output_records[dog_index][1])} aa")


if __name__ == "__main__":
    main()
