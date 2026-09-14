#!/usr/bin/env python3
"""Build the compact five-protein manual-review packet.

The packet contains the unaligned target proteins, their integrated QC rows,
and small alignment subsets with appropriate within-gene comparators. It does
not make or overwrite any biological curation decision.
"""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PHYLO = ROOT / "analyses" / "protein_analysis"
TREE = ROOT / "analyses" / "phylogenetics" / "combined_trees" / "six_gene_tree"
OUT = PHYLO / "manual_review"

TARGETS = {
    "XP_048463897.1",
    "XP_056660002.1",
    "NP_001013588.1",
    "XP_006631165.2",
    "XP_066265713.1",
}

ALIGNMENT_COMPARATORS = {
    "EPYC": {
        "XP_048463897.1",  # target: whale shark
        "XP_038636412.1",  # nearest compact-tree neighbor: catshark
        "XP_007893502.1",  # second neighbor: elephant shark
        "NP_004941.2",  # human reference
    },
    "OGN": {
        "XP_056660002.1",  # target: opossum
        "NP_001013588.1",  # target: zebrafish
        "XP_006631165.2",  # target: spotted gar
        "XP_066265713.1",  # target: amphioxus
        "NP_148935.1",  # human reference
        "NP_032786.1",  # mouse reference
        "XP_020377273.1",  # closest OGN tip to amphioxus in compact tree
    },
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def accession(header: str) -> str:
    parts = header.split("|")
    return parts[-1].split()[0]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    canonical = read_fasta(PHYLO / "candidates" / "slrp_candidates_canonical.faa")
    target_records = [record for record in canonical if accession(record[0]) in TARGETS]
    found = {accession(header) for header, _ in target_records}
    if found != TARGETS:
        raise SystemExit(f"Target mismatch: missing={sorted(TARGETS - found)}")
    write_fasta(OUT / "manual_sequence_review_5.faa", target_records)

    summary_path = (
        PHYLO / "tables" / "protein_conservation_domain_msa_signalp_summary.tsv"
    )
    with summary_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["Protein accession"] in TARGETS]
    if len(selected) != 5:
        raise SystemExit(f"Expected 5 QC rows, found {len(selected)}")
    with (OUT / "manual_sequence_review_5.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=selected[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    for gene, wanted in ALIGNMENT_COMPARATORS.items():
        path = PHYLO / "alignments" / "canonical" / f"{gene}_canonical_aligned.faa"
        records = [
            record for record in read_fasta(path) if accession(record[0]) in wanted
        ]
        present = {accession(header) for header, _ in records}
        if present != wanted:
            raise SystemExit(
                f"{gene} comparator mismatch: missing={sorted(wanted - present)}"
            )
        write_fasta(OUT / f"{gene}_manual_review_alignment_subset.faa", records)

    neighbor_path = TREE / "compact_tree_flagged_tip_neighbors.tsv"
    with neighbor_path.open(encoding="utf-8", newline="") as handle:
        neighbor_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected_neighbors = [
        row for row in neighbor_rows if accession(row["query_tip"]) in TARGETS
    ]
    with (OUT / "manual_sequence_review_5_tree_neighbors.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=neighbor_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(selected_neighbors)

    synvoy_path = (
        ROOT / "analyses" / "synteny" / "tables" / "synvoy_gene_species_evidence.tsv"
    )
    with synvoy_path.open(encoding="utf-8", newline="") as handle:
        synvoy_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected_synteny = [
        row for row in synvoy_rows if row["canonical_protein_accession"] in TARGETS
    ]
    if len(selected_synteny) != 5:
        raise SystemExit(
            f"Expected 5 SynVoy context rows, found {len(selected_synteny)}"
        )
    selected_synteny.sort(key=lambda row: (row["gene"], row["species"]))
    synteny_fields = [
        "gene",
        "species",
        "canonical_protein_accession",
        "genome_qc",
        "retained_high_goi_candidates",
        "retained_medium_goi_candidates",
        "best_confidence",
        "best_model_status",
        "best_owning_gene",
        "self_consistency_flags",
        "evidence_grade",
        "review_priority",
        "review_reasons",
    ]
    with (OUT / "manual_sequence_review_5_synteny_context.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=synteny_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {field: row[field] for field in synteny_fields} for row in selected_synteny
        )

    print(f"Wrote manual-review packet for {len(target_records)} proteins to {OUT}")


if __name__ == "__main__":
    main()
