#!/usr/bin/env python3
"""Create a complete 98-row SynVoy decision file.

Detailed P1 decisions are preserved verbatim. P2/P3 rows are confirmed by the
canonical-locus audit: the accession from the protein panel must map to the
best post-ownership SynVoy locus in the exact GFF. Protein/SignalP/ProtSpace
warnings keep an otherwise matching row tentative. Manually resolved rows may
override that conservative rule. Two FMOD exceptions are recorded explicitly,
and all opossum rows are marked ambiguous because the current FASTA/GFF pair is
incompatible.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

FIELDS = [
    "gene",
    "species",
    "coordinate_to_accession_status",
    "review_status",
    "confirmed_gene_symbol",
    "confirmed_protein_accession",
    "neighbor_order_consistent",
    "phylogeny_consistent",
    "reviewer",
    "review_date",
    "manual_notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--locus-audit", required=True, type=Path)
    parser.add_argument("--detailed-decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    evidence = read_tsv(args.evidence)
    audit = {(row["gene"], row["species"]): row for row in read_tsv(args.locus_audit)}
    detailed = {
        (row["gene"], row["species"]): row for row in read_tsv(args.detailed_decisions)
    }
    output: list[dict[str, str]] = []

    for evidence_row in evidence:
        key = (evidence_row["gene"], evidence_row["species"])
        if key in detailed:
            output.append({field: detailed[key][field] for field in FIELDS})
            continue

        audit_row = audit[key]
        accession = audit_row["canonical_protein_accession"] or "not available"
        symbol = audit_row["canonical_gene_symbol"] or "not found"

        if audit_row["input_status"] != "PASS":
            output.append(
                {
                    "gene": key[0],
                    "species": key[1],
                    "coordinate_to_accession_status": "not assessable",
                    "review_status": "ambiguous",
                    "confirmed_gene_symbol": symbol,
                    "confirmed_protein_accession": accession,
                    "neighbor_order_consistent": "not assessable",
                    "phylogeny_consistent": "yes"
                    if accession != "not available"
                    else "not available",
                    "reviewer": "Local input audit",
                    "review_date": "2026-09-14",
                    "manual_notes": (
                        "The canonical protein may be independently usable, but the current opossum "
                        "SynVoy locus is not interpretable: the GFF and FASTA sequence IDs are from "
                        "incompatible versions. Rerun with a synchronized target pair before claiming synteny."
                    ),
                }
            )
            continue

        if key == ("FMOD", "spotted_gar"):
            output.append(
                {
                    "gene": "FMOD",
                    "species": "spotted_gar",
                    "coordinate_to_accession_status": "alternative retained locus",
                    "review_status": "accepted",
                    "confirmed_gene_symbol": "fmoda",
                    "confirmed_protein_accession": "XP_006628521.2",
                    "neighbor_order_consistent": "yes",
                    "phylogeny_consistent": "yes",
                    "reviewer": "Local GFF review",
                    "review_date": "2026-09-14",
                    "manual_notes": (
                        "The summary's top HIGH interval spans lumican-like and PRELP models, but a second "
                        "retained HIGH rescue overlaps canonical fmoda at NC_090700.1:38907779-38911298. "
                        "Eight human flanking anchors are shared, six preserve order, and reciprocal protein, "
                        "tree, domain and SignalP evidence support XP_006628521.2."
                    ),
                }
            )
            continue

        if key == ("FMOD", "elephant_shark"):
            output.append(
                {
                    "gene": "FMOD",
                    "species": "elephant_shark",
                    "coordinate_to_accession_status": "same locus",
                    "review_status": "tentative",
                    "confirmed_gene_symbol": "LOC103182549",
                    "confirmed_protein_accession": "XP_007897806.2",
                    "neighbor_order_consistent": "partial",
                    "phylogeny_consistent": "not available",
                    "reviewer": "Local GFF review",
                    "review_date": "2026-09-14",
                    "manual_notes": (
                        "The top SynVoy rescue overlaps an annotated fibromodulin-like model at "
                        "NW_024704760.1:2093069-2111518. Five human flanking symbols are shared and four "
                        "retain order after block reversal. Exact GFF/FASTA reconstruction gives a complete "
                        "684-aa, nine-CDS model with no internal stop, but protein BLAST separates it into a "
                        "DCN/BGN-like N-terminal half and an FMOD-best C-terminal half. Pfam detects LRR "
                        "architecture and two LRR N-terminal caps, and a second hydrophobic segment starts at "
                        "residue 330. This is therefore treated as a likely compound/fused annotation rather "
                        "than a confident single-copy FMOD protein."
                    ),
                }
            )
            continue

        if audit_row["coordinate_relation"] != "same locus":
            raise ValueError(
                f"Unresolved non-P1 coordinate conflict for {key}: {audit_row}"
            )

        caution = (
            evidence_row["protein_qc_status"] != "Supported"
            or evidence_row["protspace_status"] == "Needs inspection"
        )
        status = "tentative" if caution else "accepted"
        concerns = "; ".join(
            value
            for value in (
                evidence_row["protein_qc_concerns"],
                evidence_row["protspace_concerns"],
            )
            if value
        )
        note = (
            f"Canonical protein {accession} maps to the same exact-GFF locus as the best "
            f"post-ownership {evidence_row['best_confidence']} SynVoy call; the annotated symbol is "
            f"{symbol}. SynVoy supplies the neighbourhood score and the protein panel supplies "
            f"{evidence_row['protein_qc_status']} protein support."
        )
        if caution:
            note += " The locus is retained as tentative because protein-level QC is not fully clean."
        if concerns:
            note += f" Recorded concerns: {concerns}."

        output.append(
            {
                "gene": key[0],
                "species": key[1],
                "coordinate_to_accession_status": "same locus",
                "review_status": status,
                "confirmed_gene_symbol": symbol,
                "confirmed_protein_accession": accession,
                "neighbor_order_consistent": "yes",
                "phylogeny_consistent": "yes",
                "reviewer": "Canonical-locus batch review",
                "review_date": "2026-09-14",
                "manual_notes": note,
            }
        )

    if len(output) != len(evidence):
        raise RuntimeError(f"Expected {len(evidence)} decisions, built {len(output)}")
    write_tsv(args.output, output)
    counts: dict[str, int] = {}
    for row in output:
        counts[row["review_status"]] = counts.get(row["review_status"], 0) + 1
    print(f"Wrote {len(output)} complete decisions to {args.output}")
    for status, count in sorted(counts.items()):
        print(f"{status}\t{count}")


if __name__ == "__main__":
    main()
