#!/usr/bin/env python3
"""Audit exact sequence-ID compatibility between paired FASTA and GFF3 files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_ids(path: Path) -> set[str]:
    # Collect exact sequence identifiers from FASTA headers.
    ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split(maxsplit=1)[0])
    return ids


def gff_ids(path: Path) -> set[str]:
    # Collect sequence IDs used by non-comment GFF feature rows.
    ids: set[str] = set()
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            if len(fields) >= 3 and fields[2] in {"gene", "mRNA", "CDS"}:
                ids.add(fields[0])
    return ids


def unversioned(accession: str) -> str:
    # Remove a final numeric accession version for secondary diagnostics only.
    head, separator, tail = accession.rpartition(".")
    return head if separator and tail.isdigit() else accession


def main() -> None:
    # Compare every FASTA/GFF pair using exact and unversioned identifiers and
    # report whether SynVoy can safely combine the files.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fna-dir", required=True, type=Path)
    parser.add_argument("--gff-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    summaries: list[dict[str, object]] = []
    mismatches: list[dict[str, str]] = []
    for fasta in sorted(args.fna_dir.glob("*.fna")):
        gff = args.gff_dir / f"{fasta.stem}.gff"
        if not gff.is_file():
            summaries.append(
                {
                    "species": fasta.stem,
                    "fasta_seqids": 0,
                    "gff_annotated_seqids": 0,
                    "exact_matches": 0,
                    "version_only_mismatches": 0,
                    "missing_from_fasta": 0,
                    "status": "missing GFF",
                }
            )
            continue

        fasta_set = fasta_ids(fasta)
        gff_set = gff_ids(gff)
        exact = gff_set & fasta_set
        fasta_by_base = {unversioned(item): item for item in fasta_set}
        version_only = 0
        absent = 0
        for seqid in sorted(gff_set - fasta_set):
            base = unversioned(seqid)
            if base in fasta_by_base:
                mismatch_type = "version_only"
                fasta_seqid = fasta_by_base[base]
                version_only += 1
            else:
                mismatch_type = "absent"
                fasta_seqid = ""
                absent += 1
            mismatches.append(
                {
                    "species": fasta.stem,
                    "gff_seqid": seqid,
                    "fasta_seqid": fasta_seqid,
                    "mismatch_type": mismatch_type,
                }
            )
        status = "PASS" if version_only == 0 and absent == 0 else "FAIL"
        summaries.append(
            {
                "species": fasta.stem,
                "fasta_seqids": len(fasta_set),
                "gff_annotated_seqids": len(gff_set),
                "exact_matches": len(exact),
                "version_only_mismatches": version_only,
                "missing_from_fasta": absent,
                "status": status,
            }
        )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = [
        (
            args.output_prefix.with_name(args.output_prefix.name + "_summary.tsv"),
            summaries,
        ),
        (
            args.output_prefix.with_name(args.output_prefix.name + "_mismatches.tsv"),
            mismatches,
        ),
    ]
    for path, rows in outputs:
        if not rows:
            path.write_text("", encoding="utf-8")
            continue
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
