#!/usr/bin/env python3
"""Build the two focused SLRP phylogeny datasets.

The class-I extension reuses the finished thesis BGN/DCN/ASPN tree dataset and
adds the shark DCN2 and cyclostome SLRP1--4 proteins published by Gil et al.
The class-II extension uses the published FMOD/LUM/LumL and cyclostome
SLRP6--8 sequences, plus the current thesis lamprey LUM and zebrafish FMOD
records that need to be checked.
"""

from __future__ import annotations

import argparse
import csv
import re
import textwrap
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


def read_fasta(path: Path, *, remove_gaps: bool = False) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks)
                    if remove_gaps:
                        sequence = sequence.replace("-", "").replace(".", "")
                    records.append(FastaRecord(header, sequence))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        sequence = "".join(chunks)
        if remove_gaps:
            sequence = sequence.replace("-", "").replace(".", "")
        records.append(FastaRecord(header, sequence))
    return records


def write_fasta(path: Path, records: list[FastaRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            handle.write("\n".join(textwrap.wrap(record.sequence, 80)) + "\n")


def normalize_accession(value: str) -> str:
    token = value.strip().split()[0]
    token = token.split(".")[0]
    return re.sub(r"[^A-Za-z0-9]", "", token).upper()


def safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")


def normalized_group(value: str) -> str:
    value = value.strip()
    slrp = re.search(r"SLRP(\d+)", value, flags=re.IGNORECASE)
    if slrp:
        return f"SLRP{slrp.group(1)}"
    return value.upper()


def published_alignment_index(path: Path) -> dict[str, FastaRecord]:
    index: dict[str, FastaRecord] = {}
    for record in read_fasta(path, remove_gaps=True):
        accession_token = record.header.split("-", 1)[0]
        key = normalize_accession(accession_token)
        if key:
            index[key] = record
    return index


def read_target_rows(path: Path, analysis: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row for row in rows if row["analysis"] == analysis]


def accession_from_local_header(header: str) -> str:
    fields = header.split()[0].split("|")
    return fields[2] if len(fields) >= 3 else fields[-1]


def validate_record(record: FastaRecord) -> None:
    if len(record.sequence) < 100:
        raise ValueError(f"Sequence is unexpectedly short: {record.header} ({len(record.sequence)} aa)")
    invalid = sorted(set(record.sequence.upper()) - set("ABCDEFGHIKLMNPQRSTVWXYZUO*"))
    if invalid:
        raise ValueError(f"Unexpected sequence characters in {record.header}: {invalid}")


def add_record(
    output: list[FastaRecord],
    manifest: list[dict[str, str]],
    seen_headers: set[str],
    record: FastaRecord,
    row: dict[str, str],
) -> None:
    validate_record(record)
    if record.header in seen_headers:
        return
    seen_headers.add(record.header)
    output.append(record)
    manifest.append(row)


def build_class1(root: Path, published: dict[str, FastaRecord], targets: Path) -> tuple[list[FastaRecord], list[dict[str, str]]]:
    output: list[FastaRecord] = []
    manifest: list[dict[str, str]] = []
    seen_headers: set[str] = set()
    local_path = root / "analyses/phylogenetics/combined_trees/rooted_class_i_tree/results/rooted_class_i_dataset.faa"

    for local in read_fasta(local_path):
        gene = local.header.split("|", 1)[0]
        if gene not in {"BGN", "DCN", "ASPN", "FMOD", "PRELP"}:
            continue
        accession = accession_from_local_header(local.header)
        add_record(
            output,
            manifest,
            seen_headers,
            local,
            {
                "analysis": "class1",
                "role": "outgroup" if gene in {"FMOD", "PRELP"} else "thesis_ingroup",
                "group": gene,
                "species": local.header.split("|")[1],
                "accession": accession,
                "length_aa": str(len(local.sequence)),
                "source": str(local_path),
                "note": "existing curated thesis sequence",
            },
        )

    for row in read_target_rows(targets, "class1"):
        key = normalize_accession(row["accession"])
        source = published.get(key)
        if source is None:
            raise KeyError(f"Published alignment does not contain {row['accession']}")
        group = normalized_group(row["group"])
        header = f"{group}|{safe_label(row['species']).lower()}|{row['accession']}|Gil2025"
        record = FastaRecord(header, source.sequence)
        add_record(
            output,
            manifest,
            seen_headers,
            record,
            {
                "analysis": "class1",
                "role": "published_ingroup_extension",
                "group": group,
                "species": row["species"],
                "accession": row["accession"],
                "length_aa": str(len(record.sequence)),
                "source": "Gil et al. 2025 Supplementary File 1/Table 1",
                "note": row["database_annotation"],
            },
        )
    return output, manifest


def build_class2(
    root: Path, published: dict[str, FastaRecord], targets: Path
) -> tuple[list[FastaRecord], list[dict[str, str]], list[dict[str, str]]]:
    output: list[FastaRecord] = []
    manifest: list[dict[str, str]] = []
    exact_matches: list[dict[str, str]] = []
    seen_headers: set[str] = set()

    for row in read_target_rows(targets, "class2"):
        key = normalize_accession(row["accession"])
        source = published.get(key)
        if source is None:
            raise KeyError(f"Published alignment does not contain {row['accession']}")
        group = normalized_group(row["group"])
        header = f"{group}|{safe_label(row['species']).lower()}|{row['accession']}|Gil2025"
        record = FastaRecord(header, source.sequence)
        add_record(
            output,
            manifest,
            seen_headers,
            record,
            {
                "analysis": "class2",
                "role": row["role"],
                "group": group,
                "species": row["species"],
                "accession": row["accession"],
                "length_aa": str(len(record.sequence)),
                "source": "Gil et al. 2025 Supplementary File 1/Table 1",
                "note": row["database_annotation"],
            },
        )

    focus = (
        ("FMOD", "zebrafish", "FMODB_CHECK"),
        ("LUM", "lamprey", "LUM_CHECK"),
    )
    for gene, species, label in focus:
        local_path = root / f"analyses/protein_analysis/candidates/canonical_by_gene/{gene}_canonical.faa"
        candidates = [
            record
            for record in read_fasta(local_path)
            if record.header.startswith(f"{species}|{gene}|")
        ]
        if len(candidates) != 1:
            raise ValueError(f"Expected one {species} {gene} record, found {len(candidates)}")
        local = candidates[0]
        accession = accession_from_local_header(local.header)
        published_matches: list[tuple[FastaRecord, dict[str, str], int]] = []
        for record, row in zip(output, manifest):
            matcher = SequenceMatcher(None, record.sequence, local.sequence, autojunk=False)
            shared = sum(block.size for block in matcher.get_matching_blocks())
            if shared == min(len(record.sequence), len(local.sequence)):
                published_matches.append((record, row, shared))
        if published_matches:
            published_record, published_row, shared = published_matches[0]
            exact_matches.append(
                {
                    "thesis_gene": gene,
                    "thesis_species": species,
                    "thesis_accession": accession,
                    "published_group": published_row["group"],
                    "published_species": published_row["species"],
                    "published_accession": published_row["accession"],
                    "thesis_length_aa": str(len(local.sequence)),
                    "published_length_aa": str(len(published_record.sequence)),
                    "shared_residues": str(shared),
                    "overlap_identity_percent": "100.00",
                    "length_difference_aa": str(abs(len(local.sequence) - len(published_record.sequence))),
                    "tree_record": published_record.header,
                    "interpretation": "same protein core with only inserted residues; one copy retained in tree",
                }
            )
            continue
        record = FastaRecord(f"{label}|{species}|{accession}|thesis", local.sequence)
        add_record(
            output,
            manifest,
            seen_headers,
            record,
            {
                "analysis": "class2",
                "role": "thesis_focus_record",
                "group": label,
                "species": species,
                "accession": accession,
                "length_aa": str(len(record.sequence)),
                "source": str(local_path),
                "note": "current thesis candidate added for direct placement check",
            },
        )
    return output, manifest, exact_matches


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["analysis", "role", "group", "species", "accession", "length_aa", "source", "note"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--analysis", choices=("class1", "class2", "all"), default="all")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    extension_root = root / "analyses/phylogenetics/targeted_extensions"
    reference = extension_root / "reference"
    alignment_path = reference / "Supplemental_File_S1_G3-2024-405575.unk"
    targets_path = reference / "gil2025_target_sequences.tsv"
    published = published_alignment_index(alignment_path)

    analyses = ("class1", "class2") if args.analysis == "all" else (args.analysis,)
    for analysis in analyses:
        if analysis == "class1":
            records, manifest = build_class1(root, published, targets_path)
            exact_matches: list[dict[str, str]] = []
        else:
            records, manifest, exact_matches = build_class2(root, published, targets_path)
        write_fasta(args.output_dir / f"{analysis}_targeted_dataset.faa", records)
        write_manifest(args.output_dir / f"{analysis}_targeted_manifest.tsv", manifest)
        if exact_matches:
            write_tsv_path = args.output_dir / f"{analysis}_exact_sequence_matches.tsv"
            write_tsv_path.parent.mkdir(parents=True, exist_ok=True)
            with write_tsv_path.open("w", encoding="utf-8", newline="") as handle:
                fields = list(exact_matches[0])
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(exact_matches)
        print(f"{analysis}: {len(records)} sequences")


if __name__ == "__main__":
    main()
