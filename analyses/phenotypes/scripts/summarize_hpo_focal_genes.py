#!/usr/bin/env python3
"""Subset the pinned HPO human gene–phenotype and gene–disease files."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-phenotype", required=True, type=Path)
    parser.add_argument("--gene-disease", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    phenotypes = [
        row for row in read_tsv(args.gene_phenotype) if row["gene_symbol"] in GENES
    ]
    diseases = [
        row for row in read_tsv(args.gene_disease) if row["gene_symbol"] in GENES
    ]
    write_tsv(
        args.output_dir / "hpo_human_focal_gene_phenotypes.tsv",
        phenotypes,
        list(phenotypes[0])
        if phenotypes
        else [
            "ncbi_gene_id",
            "gene_symbol",
            "hpo_id",
            "hpo_name",
            "frequency",
            "disease_id",
        ],
    )
    write_tsv(
        args.output_dir / "hpo_human_focal_gene_diseases.tsv",
        diseases,
        list(diseases[0])
        if diseases
        else [
            "ncbi_gene_id",
            "gene_symbol",
            "association_type",
            "disease_id",
            "source",
        ],
    )

    by_gene_terms: dict[str, set[str]] = defaultdict(set)
    by_gene_diseases: dict[str, set[str]] = defaultdict(set)
    by_gene_mendelian: dict[str, set[str]] = defaultdict(set)
    for row in phenotypes:
        by_gene_terms[row["gene_symbol"]].add(row["hpo_id"])
        by_gene_diseases[row["gene_symbol"]].add(row["disease_id"])
    for row in diseases:
        by_gene_diseases[row["gene_symbol"]].add(row["disease_id"])
        if row["association_type"] == "MENDELIAN":
            by_gene_mendelian[row["gene_symbol"]].add(row["disease_id"])

    summary: list[dict[str, object]] = []
    for gene in GENES:
        n_terms = len(by_gene_terms[gene])
        summary.append(
            {
                "gene": gene,
                "unique_hpo_terms": n_terms,
                "disease_ids": ";".join(sorted(by_gene_diseases[gene])),
                "mendelian_disease_ids": ";".join(sorted(by_gene_mendelian[gene])),
                "interpretation": (
                    "direct human disease/phenotype support is present in HPO"
                    if n_terms
                    else "no HPO gene–phenotype record in this release; this is absence of curated disease annotation, not evidence of no biological function"
                ),
            }
        )
    write_tsv(
        args.output_dir / "hpo_human_phenotype_evidence_summary.tsv",
        summary,
        list(summary[0]),
    )
    print(
        f"Recovered {len(phenotypes)} human phenotype rows and {len(diseases)} disease mappings for the seven genes"
    )


if __name__ == "__main__":
    main()
