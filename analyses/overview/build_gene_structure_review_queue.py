#!/usr/bin/env python3
"""Prioritize manual transcript checks for the five-species gene-structure panel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "analyses/overview/gene_structure_manual_review_queue.tsv",
    )
    parser.add_argument(
        "--ncbi-review",
        type=Path,
        default=ROOT
        / "analyses/overview/gene_structure_predicted_transcript_ncbi_review.tsv",
    )
    args = parser.parse_args()

    with args.representatives.open(encoding="utf-8", newline="") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    if len(source) != 35:
        raise ValueError("Expected 35 representative gene/species rows")

    ncbi_review: dict[str, dict[str, str]] = {}
    if args.ncbi_review.exists():
        with args.ncbi_review.open(encoding="utf-8", newline="") as handle:
            ncbi_review = {
                row["transcript_accession"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }

    rows: list[dict[str, str]] = []
    for row in source:
        transcript = row["transcript_id"]
        if row["structural_flag_clean"]:
            priority, model = "P1", "structural exception"
        elif "NM_" in transcript:
            priority, model = "P3", "curated RefSeq transcript"
        elif "XM_" in transcript:
            priority, model = "P2", "predicted RefSeq transcript"
        else:
            priority, model = "P1", "non-NM/XM annotation model"
        accession = transcript.removeprefix("rna-")
        audit = ncbi_review.get(accession)
        if audit and audit["review_status"] == "unresolved annotation conflict":
            priority, model = "P1", "predicted RefSeq with gene/product conflict"
        rows.append(
            {
                "priority": priority,
                "gene": row["query_gene"],
                "panel_role": "context" if row["query_gene"] == "OMD" else "main",
                "species": row["species_or_file"],
                "transcript_id": transcript,
                "model_type": model,
                "cds_exon_count": row["cds_exon_count"],
                "protein_length_est_aa": row["protein_length_est_aa"],
                "gene_span_bp": row["gene_span_bp"],
                "structural_flag": row["structural_flag_clean"],
                "manual_action": (
                    "Resolve the chicken BGN-locus/asporin-product conflict with neighborhood/synteny evidence; "
                    "do not assert XP_414298.2 as a BGN protein ortholog"
                    if audit
                    and audit["review_status"] == "unresolved annotation conflict"
                    else (
                        "Database annotation review complete; optional NCBI Genome Data Viewer spot-check"
                        if audit
                        else "Confirm gene symbol, transcript accession/status, and CDS exon count in NCBI Gene/RefSeq; "
                        "record whether this is the preferred/canonical model"
                    )
                ),
                "review_status": audit["review_status"] if audit else "",
                "reviewer_notes": audit["review_conclusion"] if audit else "",
            }
        )
    rank = {"P1": 1, "P2": 2, "P3": 3}
    rows.sort(key=lambda row: (rank[row["priority"]], row["gene"], row["species"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts = {
        priority: sum(row["priority"] == priority for row in rows)
        for priority in ("P1", "P2", "P3")
    }
    print(f"Wrote {args.output} ({len(rows)} rows; {counts})")


if __name__ == "__main__":
    main()
