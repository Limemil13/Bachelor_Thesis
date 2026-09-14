#!/usr/bin/env python3
"""Merge retained SignalP jobs and keep only current canonical proteins.

Missing canonical rows are intentional: newly promoted or replaced proteins
remain explicitly pending until a new SignalP run is imported.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing", required=True, type=Path)
    parser.add_argument("--lum", required=True, type=Path)
    parser.add_argument(
        "--supplemental",
        action="append",
        default=[],
        type=Path,
        help="Additional parsed SignalP table; may be supplied more than once",
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    existing = read_tsv(args.existing)
    lum = read_tsv(args.lum)
    if len(existing) != 77 or {row["gene"] for row in existing} != {
        "BGN",
        "EPYC",
        "FMOD",
        "OGN",
        "PRELP",
    }:
        raise ValueError(
            "Existing SignalP summary is not the validated five-gene 77-row table"
        )
    if len(lum) != 16 or {row["gene"] for row in lum} != {"LUM"}:
        raise ValueError("LUM SignalP summary must contain exactly 16 LUM rows")

    canonical = read_tsv(args.manifest)
    canonical_keys = {
        (row["gene"], row["species"], row["accession"]) for row in canonical
    }
    supplemental = [
        row for supplemental_path in args.supplemental for row in read_tsv(supplemental_path)
    ]
    historical = existing + lum + supplemental
    rows = [
        row
        for row in historical
        if (row["gene"], row["species"], row["accession"]) in canonical_keys
    ]
    keys = [(row["gene"], row["species"], row["accession"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Merged SignalP summary contains duplicate canonical proteins")
    if not set(keys).issubset(canonical_keys):
        raise ValueError("Merged SignalP summary contains noncanonical proteins")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    write_tsv(args.out_dir / "signalp_summary_clean.tsv", rows, fields)
    write_tsv(
        args.out_dir / "signalp_candidates_to_check.tsv",
        [row for row in rows if row["strict_check"] == "yes"],
        fields,
    )
    write_tsv(
        args.out_dir / "signalp_candidates_watch.tsv",
        [row for row in rows if row["watch"] == "yes"],
        fields,
    )
    print(f"Merged SignalP rows: {len(rows)}")
    print(f"Signal peptide yes: {sum(row['signal_peptide'] == 'yes' for row in rows)}")
    print(f"Signal peptide no: {sum(row['signal_peptide'] == 'no' for row in rows)}")
    print(f"Canonical proteins pending SignalP: {len(canonical_keys) - len(rows)}")


if __name__ == "__main__":
    main()
