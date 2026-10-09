#!/usr/bin/env python3
"""Cross-reference gene-structure representatives with the canonical panel."""

from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[3]
STRUCTURE_FASTA = (
    ROOT
    / "analyses/gene_structure/extended/sequences/representative_proteins_5species.faa"
)
CANONICAL_DIR = ROOT / "analyses/protein_analysis/candidates/canonical_by_gene"
OUTPUT = (
    ROOT / "analyses/gene_structure/extended/tables/structure_protein_crosswalk.tsv"
)
STRUCTURE_TABLE = (
    ROOT / "analyses/gene_structure/tables/gene_structure_representative_5species.tsv"
)

SPECIES_TO_CANONICAL = {
    "human": "human",
    "mouse": "mouse",
    "cow": "bos_taurus",
    "chicken": "chicken",
    "zebrafish": "zebrafish",
}


def canonical_records() -> dict[tuple[str, str], object]:
    # Collect curated protein records under a common (gene, species) key.
    records: dict[tuple[str, str], object] = {}
    for path in sorted(CANONICAL_DIR.glob("*_canonical.faa")):
        for record in SeqIO.parse(path, "fasta"):
            species, gene, accession = record.id.split("|")[:3]
            record.annotations["accession"] = accession
            records[(gene, species)] = record
    return records


def structure_symbols() -> dict[tuple[str, str], str]:
    # Preserve the gene symbol attached to each selected structure transcript.
    with STRUCTURE_TABLE.open(encoding="utf-8", newline="") as handle:
        return {
            (row["query_gene"], row["species_or_file"]): row["matched_gene_name"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def alignment_metrics(sequence_a: str, sequence_b: str) -> tuple[float, float, float]:
    # Global alignment tests full-length agreement rather than a short shared domain.
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(sequence_a, sequence_b)[0]

    matches = 0
    aligned_residues = 0
    for (a_start, a_end), (b_start, b_end) in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        block_a = sequence_a[a_start:a_end]
        block_b = sequence_b[b_start:b_end]
        aligned_residues += len(block_a)
        matches += sum(a == b for a, b in zip(block_a, block_b, strict=True))

    identity = 100.0 * matches / aligned_residues if aligned_residues else 0.0
    coverage_a = 100.0 * aligned_residues / len(sequence_a)
    coverage_b = 100.0 * aligned_residues / len(sequence_b)
    return identity, coverage_a, coverage_b


def main() -> None:
    # Match every structure-derived translation to the curated protein for the
    # same gene and species, then report accession and sequence agreement.
    canonical = canonical_records()
    symbols = structure_symbols()
    # One row is written for every representative protein in the structure panel.
    rows: list[dict[str, object]] = []

    for record in SeqIO.parse(STRUCTURE_FASTA, "fasta"):
        gene, species, structure_accession = record.id.split("|")[:3]
        structure_symbol = symbols[(gene, species)]
        canonical_species = SPECIES_TO_CANONICAL[species]
        comparison = canonical.get((gene, canonical_species))

        if gene == "OMD":
            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "structure_gene_symbol": structure_symbol,
                    "structure_accession": structure_accession,
                    "canonical_accession": "",
                    "accession_relation": "context_only_not_in_canonical_panel",
                    "aligned_identity_percent": "",
                    "structure_sequence_coverage_percent": "",
                    "canonical_sequence_coverage_percent": "",
                    "orthology_use": "structure context only",
                }
            )
            continue

        if comparison is None:
            note = "not orthology support"
            relation = "not_in_canonical_panel"
            if gene == "BGN" and species == "chicken":
                relation = "excluded_ASPN_like_product"
            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "structure_gene_symbol": structure_symbol,
                    "structure_accession": structure_accession,
                    "canonical_accession": "",
                    "accession_relation": relation,
                    "aligned_identity_percent": "",
                    "structure_sequence_coverage_percent": "",
                    "canonical_sequence_coverage_percent": "",
                    "orthology_use": note,
                }
            )
            continue

        canonical_accession = comparison.annotations["accession"]
        identity, structure_coverage, canonical_coverage = alignment_metrics(
            str(record.seq), str(comparison.seq)
        )
        relation = (
            "exact_accession"
            if structure_accession == canonical_accession
            else "different_accession_same_annotated_gene"
        )
        orthology_use = "compatible with canonical protein panel"
        if gene == "BGN" and species == "zebrafish":
            relation = "different_teleost_coortholog"
            orthology_use = "BGNa structural representative; canonical panel uses BGNb"
        rows.append(
            {
                "gene": gene,
                "species": species,
                "structure_gene_symbol": structure_symbol,
                "structure_accession": structure_accession,
                "canonical_accession": canonical_accession,
                "accession_relation": relation,
                "aligned_identity_percent": f"{identity:.2f}",
                "structure_sequence_coverage_percent": f"{structure_coverage:.2f}",
                "canonical_sequence_coverage_percent": f"{canonical_coverage:.2f}",
                "orthology_use": orthology_use,
            }
        )

    # The crosswalk makes missing or discordant pairs visible for manual review.
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    exact = sum(row["accession_relation"] == "exact_accession" for row in rows)
    alternative = sum(
        row["accession_relation"] == "different_accession_same_annotated_gene"
        for row in rows
    )
    print(
        f"Wrote {len(rows)} rows to {OUTPUT}; "
        f"exact accessions={exact}, alternative accessions={alternative}."
    )


if __name__ == "__main__":
    main()
