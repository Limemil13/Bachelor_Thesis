#!/usr/bin/env python3
"""Select one biologically plausible transcript per gene/species GFF result.

The selector is the reproducible, parameterized replacement for the historical
six-gene script in the SynVoy checkout.  It keeps the same output columns and
quality flags while including LUM explicitly.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED_PROTEIN_LENGTH = {
    "BGN": 370,
    "DCN": 359,
    "EPYC": 323,
    "FMOD": 377,
    "LUM": 338,
    "OGN": 300,
    "OMD": 423,
    "PRELP": 382,
}


def as_int(value: str | None) -> int | None:
    try:
        return int(float(value)) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def transcript_rank(transcript_id: str) -> int:
    if "NM_" in transcript_id:
        return 2
    if "XM_" in transcript_id:
        return 1
    return 0


def candidate_rank(
    row: dict[str, str], expected_length: int
) -> tuple[int, int, int, int]:
    protein_length = as_int(row.get("protein_length_est_aa"))
    distance = (
        abs(protein_length - expected_length) if protein_length is not None else 9999
    )
    return (
        distance,
        -transcript_rank(row.get("transcript_id", "")),
        -(as_int(row.get("cds_exon_count")) or 0),
        -(as_int(row.get("cds_length_bp")) or 0),
    )


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No gene-structure rows in {args.input}")

    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["query_gene"], row["species_or_file"])].append(row)

    representatives: list[dict[str, str]] = []
    for (gene, _species), candidates in groups.items():
        expected = EXPECTED_PROTEIN_LENGTH[gene]
        representatives.append(
            min(
                candidates,
                key=lambda row, expected_length=expected: candidate_rank(
                    row, expected_length
                ),
            )
        )

    exon_counts: dict[str, list[int]] = defaultdict(list)
    for row in representatives:
        value = as_int(row.get("cds_exon_count"))
        if value is not None:
            exon_counts[row["query_gene"]].append(value)
    modal_exons = {
        gene: Counter(values).most_common(1)[0][0]
        for gene, values in exon_counts.items()
    }

    for row in representatives:
        gene = row["query_gene"]
        cds_exons = as_int(row.get("cds_exon_count"))
        protein_length = as_int(row.get("protein_length_est_aa"))
        gene_span = as_int(row.get("gene_span_bp"))
        mode = modal_exons.get(gene)
        flags: list[str] = []
        if cds_exons is not None and mode is not None and cds_exons != mode:
            flags.append("cds_exon_count_differs_from_gene_mode")
        if protein_length is not None and not 250 <= protein_length <= 500:
            flags.append("unusual_slrp_protein_length")
        if gene_span is not None and gene_span > 100_000:
            flags.append("very_large_gene_span")
        row["gene_modal_cds_exon_count"] = str(mode) if mode is not None else ""
        row["cds_exon_count_deviation"] = (
            str(cds_exons - mode) if cds_exons is not None and mode is not None else ""
        )
        row["structural_flag_clean"] = ";".join(flags)

    representatives.sort(key=lambda row: (row["query_gene"], row["species_or_file"]))
    fields = [
        "query_gene",
        "species_or_file",
        "matched_gene_name",
        "transcript_id",
        "seqid",
        "strand",
        "gene_start",
        "gene_end",
        "gene_span_bp",
        "cds_start",
        "cds_end",
        "cds_span_bp",
        "exon_count",
        "cds_exon_count",
        "cds_length_bp",
        "protein_length_est_aa",
        "intron_count_between_cds",
        "median_cds_intron_length_bp",
        "gene_modal_cds_exon_count",
        "cds_exon_count_deviation",
        "structural_flag_clean",
        "note",
        "gff_file",
    ]
    write_tsv(
        args.out_dir / "main_slrp_gene_structure_5species_representative_clean.tsv",
        representatives,
        fields,
    )

    summary: list[dict[str, str | int]] = []
    for gene in sorted({row["query_gene"] for row in representatives}):
        gene_rows = [row for row in representatives if row["query_gene"] == gene]
        proteins = [
            value
            for row in gene_rows
            if (value := as_int(row.get("protein_length_est_aa"))) is not None
        ]
        spans = [
            value
            for row in gene_rows
            if (value := as_int(row.get("gene_span_bp"))) is not None
        ]
        counts = sorted(
            {
                value
                for row in gene_rows
                if (value := as_int(row.get("cds_exon_count"))) is not None
            }
        )
        summary.append(
            {
                "query_gene": gene,
                "n_species": len({row["species_or_file"] for row in gene_rows}),
                "species": ",".join(
                    sorted({row["species_or_file"] for row in gene_rows})
                ),
                "cds_exon_counts": ",".join(map(str, counts)),
                "modal_cds_exon_count": modal_exons.get(gene, ""),
                "protein_length_min": min(proteins) if proteins else "",
                "protein_length_max": max(proteins) if proteins else "",
                "gene_span_min_bp": min(spans) if spans else "",
                "gene_span_max_bp": max(spans) if spans else "",
                "structural_outlier_species": ",".join(
                    row["species_or_file"]
                    for row in gene_rows
                    if row["structural_flag_clean"]
                ),
            }
        )
    summary_fields = [
        "query_gene",
        "n_species",
        "species",
        "cds_exon_counts",
        "modal_cds_exon_count",
        "protein_length_min",
        "protein_length_max",
        "gene_span_min_bp",
        "gene_span_max_bp",
        "structural_outlier_species",
    ]
    write_tsv(
        args.out_dir / "main_slrp_gene_structure_5species_summary_clean.tsv",
        summary,
        summary_fields,
    )
    write_tsv(
        args.out_dir
        / "main_slrp_gene_structure_5species_structural_outliers_clean.tsv",
        [row for row in representatives if row["structural_flag_clean"]],
        fields,
    )

    expected_pairs = {
        (gene, species)
        for gene in EXPECTED_PROTEIN_LENGTH
        for species in {row["species_or_file"] for row in rows}
    }
    observed_pairs = {
        (row["query_gene"], row["species_or_file"]) for row in representatives
    }
    missing = [
        {"query_gene": gene, "species_or_file": species}
        for gene, species in sorted(expected_pairs - observed_pairs)
    ]
    write_tsv(
        args.out_dir / "main_slrp_gene_structure_5species_missing.tsv",
        missing,
        ["query_gene", "species_or_file"],
    )

    print(f"Representatives: {len(representatives)}")
    print(f"Genes: {len(summary)}")
    print(
        f"Structural outliers: {sum(bool(row['structural_flag_clean']) for row in representatives)}"
    )
    print(f"Missing gene/species pairs: {len(missing)}")


if __name__ == "__main__":
    main()
