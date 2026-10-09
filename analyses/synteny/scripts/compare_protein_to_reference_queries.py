#!/usr/bin/env python3
"""Run one candidate protein against named reference FASTAs with BLASTP."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

OUTFMT = (
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
)


def parse_subject(value: str) -> tuple[str, Path]:
    # Accept readable NAME=FASTA arguments for multiple reference panels.
    if "=" not in value:
        raise argparse.ArgumentTypeError("Use NAME=PATH for --subject")
    name, path = value.split("=", 1)
    return name, Path(path)


def main() -> None:
    # Run the candidate against every named reference FASTA and keep all tabular
    # BLASTP fields needed to compare identity, coverage and score.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--subject", required=True, action="append", type=parse_subject)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    for label, subject in args.subject:
        command = [
            "blastp",
            "-query",
            str(args.query),
            "-subject",
            str(subject),
            "-outfmt",
            OUTFMT,
            "-max_hsps",
            "1",
            "-max_target_seqs",
            "1",
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        line = result.stdout.strip().splitlines()
        if not line:
            rows.append({"reference_gene": label, "status": "no hit"})
            continue
        fields = line[0].split("\t")
        (
            query_id,
            subject_id,
            identity,
            alignment_length,
            query_length,
            subject_length,
            query_start,
            query_end,
            subject_start,
            subject_end,
            evalue,
            bit_score,
        ) = fields
        rows.append(
            {
                "reference_gene": label,
                "status": "hit",
                "candidate_id": query_id,
                "reference_id": subject_id,
                "identity_percent": identity,
                "alignment_length": alignment_length,
                "candidate_length": query_length,
                "reference_length": subject_length,
                "candidate_coverage_percent": round(
                    100 * int(alignment_length) / int(query_length), 2
                ),
                "reference_coverage_percent": round(
                    100 * int(alignment_length) / int(subject_length), 2
                ),
                "candidate_start": query_start,
                "candidate_end": query_end,
                "reference_start": subject_start,
                "reference_end": subject_end,
                "evalue": evalue,
                "bit_score": bit_score,
            }
        )

    rows.sort(key=lambda row: float(row.get("bit_score", 0)), reverse=True)
    fields = [
        "reference_gene",
        "status",
        "candidate_id",
        "reference_id",
        "identity_percent",
        "alignment_length",
        "candidate_length",
        "reference_length",
        "candidate_coverage_percent",
        "reference_coverage_percent",
        "candidate_start",
        "candidate_end",
        "reference_start",
        "reference_end",
        "evalue",
        "bit_score",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Compared {args.query} with {len(rows)} reference queries")


if __name__ == "__main__":
    main()
