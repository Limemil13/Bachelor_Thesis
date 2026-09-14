#!/usr/bin/env python3
"""Build a reproducible seven-gene mouse phenotype evidence matrix from MGI reports."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
SYMBOL_MAP = {gene.lower(): gene for gene in GENES}
GENE_PATTERN = re.compile(r"(?:^|[+/])\s*([A-Za-z0-9_.-]+)<")

CATEGORY_KEYWORDS = {
    "skeletal/cartilage/joint": (
        "skeletal",
        "skeleton",
        "bone",
        "cartilage",
        "chondro",
        "joint",
        "ossification",
        "osteo",
        "vertebr",
        "rib ",
        "femur",
        "tibia",
        "mandible",
        "craniofacial",
        "tooth",
        "teeth",
        "tendon",
        "ligament",
    ),
    "growth/body size": ("growth", "body size", "body weight", "dwarf", "length"),
    "eye/cornea": ("eye", "ocular", "cornea", "lens", "retina", "sclera"),
    "skin/connective tissue": (
        "skin",
        "derm",
        "connective tissue",
        "collagen",
        "wound",
        "fibrosis",
    ),
    "cardiovascular": (
        "cardiovascular",
        "heart",
        "cardiac",
        "vascular",
        "blood vessel",
        "aorta",
    ),
    "immune/inflammatory": (
        "immune",
        "inflamm",
        "leukocyte",
        "lymphocyte",
        "macrophage",
        "cytokine",
    ),
    "renal/metabolic": (
        "kidney",
        "renal",
        "metabolic",
        "glucose",
        "insulin",
        "adipose",
        "lipid",
    ),
    "hearing/ear": ("hearing", "auditory", "cochlea", "inner ear"),
    "reproductive/developmental": (
        "fertility",
        "reproductive",
        "embryo",
        "development",
        "gestation",
        "placenta",
    ),
}


def parse_obo(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    names: dict[str, str] = {}
    parents: dict[str, list[str]] = defaultdict(list)
    current: dict[str, object] | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current and current.get("id"):
                    term_id = str(current["id"])
                    names[term_id] = str(current.get("name", ""))
                    parents[term_id].extend(current.get("parents", []))
                current = {"parents": []}
            elif line.startswith("["):
                if current and current.get("id"):
                    term_id = str(current["id"])
                    names[term_id] = str(current.get("name", ""))
                    parents[term_id].extend(current.get("parents", []))
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("is_a: "):
                    current["parents"].append(line[6:].split()[0])
        if current and current.get("id"):
            term_id = str(current["id"])
            names[term_id] = str(current.get("name", ""))
            parents[term_id].extend(current.get("parents", []))
    return names, parents


def ancestor_names(
    term_id: str, names: dict[str, str], parents: dict[str, list[str]]
) -> list[str]:
    seen: set[str] = set()
    stack = [term_id]
    result: list[str] = []
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        if current in names:
            result.append(names[current])
        stack.extend(parents.get(current, []))
    return result


def categories_for(
    term_id: str, names: dict[str, str], parents: dict[str, list[str]]
) -> list[str]:
    context = " | ".join(ancestor_names(term_id, names, parents)).lower()
    categories = [
        category
        for category, keywords in CATEGORY_KEYWORDS.items()
        if any(keyword in context for keyword in keywords)
    ]
    return categories or ["other"]


def parse_mgi(
    path: Path, names: dict[str, str], parents: dict[str, list[str]]
) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < 8:
                continue
            (
                genotype,
                alleles,
                allele_ids,
                background,
                mp_id,
                pmid,
                marker_id,
                genotype_id,
            ) = fields[:8]
            symbols = {
                match.group(1).lower() for match in GENE_PATTERN.finditer(genotype)
            }
            focal = sorted(
                SYMBOL_MAP[symbol] for symbol in symbols if symbol in SYMBOL_MAP
            )
            if not focal:
                continue
            all_gene_symbols = sorted(
                {match.group(1) for match in GENE_PATTERN.finditer(genotype)}
            )
            evidence_scope = (
                "single_gene"
                if len({s.lower() for s in all_gene_symbols}) == 1
                else "compound_genotype"
            )
            categories = categories_for(mp_id, names, parents)
            for gene in focal:
                records.append(
                    {
                        "gene": gene,
                        "genotype": genotype,
                        "alleles": alleles,
                        "allele_ids": allele_ids,
                        "background": background,
                        "mp_id": mp_id,
                        "mp_term": names.get(mp_id, "unresolved MP term"),
                        "phenotype_categories": ";".join(categories),
                        "pmid": pmid,
                        "marker_mgi_id": marker_id,
                        "genotype_mgi_id": genotype_id,
                        "evidence_scope": evidence_scope,
                        "genes_in_genotype": ";".join(all_gene_symbols),
                    }
                )
    return records


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def aggregate(
    records: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    term_data: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in records:
        key = (row["gene"], row["mp_id"], row["evidence_scope"])
        if key not in term_data:
            term_data[key] = {
                "gene": row["gene"],
                "mp_id": row["mp_id"],
                "mp_term": row["mp_term"],
                "evidence_scope": row["evidence_scope"],
                "phenotype_categories": row["phenotype_categories"],
                "genotypes": set(),
                "pmids": set(),
                "record_count": 0,
            }
        item = term_data[key]
        item["genotypes"].add(row["genotype_mgi_id"])
        if row["pmid"]:
            item["pmids"].update(x for x in row["pmid"].split(",") if x)
        item["record_count"] += 1

    terms: list[dict[str, object]] = []
    for item in term_data.values():
        terms.append(
            {
                "gene": item["gene"],
                "mp_id": item["mp_id"],
                "mp_term": item["mp_term"],
                "evidence_scope": item["evidence_scope"],
                "phenotype_categories": item["phenotype_categories"],
                "genotype_count": len(item["genotypes"]),
                "record_count": item["record_count"],
                "pmids": ";".join(sorted(item["pmids"])),
                "publication_count": len(item["pmids"]),
            }
        )
    terms.sort(
        key=lambda row: (
            GENES.index(str(row["gene"])),
            str(row["evidence_scope"]),
            str(row["mp_id"]),
        )
    )

    category_rows: list[dict[str, object]] = []
    category_order = list(CATEGORY_KEYWORDS) + ["other"]
    for gene in GENES:
        for category in category_order:
            for scope in ("single_gene", "compound_genotype"):
                selected = [
                    row
                    for row in terms
                    if row["gene"] == gene
                    and row["evidence_scope"] == scope
                    and category in str(row["phenotype_categories"]).split(";")
                ]
                pmids = {
                    p for row in selected for p in str(row["pmids"]).split(";") if p
                }
                category_rows.append(
                    {
                        "gene": gene,
                        "category": category,
                        "evidence_scope": scope,
                        "unique_mp_terms": len(selected),
                        "genotype_count": len(
                            {
                                record["genotype_mgi_id"]
                                for record in records
                                if record["gene"] == gene
                                and record["evidence_scope"] == scope
                                and category
                                in record["phenotype_categories"].split(";")
                            }
                        ),
                        "publication_count": len(pmids),
                    }
                )

    evidence_rows: list[dict[str, object]] = []
    for gene in GENES:
        single = [
            row
            for row in terms
            if row["gene"] == gene and row["evidence_scope"] == "single_gene"
        ]
        compound = [
            row
            for row in terms
            if row["gene"] == gene and row["evidence_scope"] == "compound_genotype"
        ]
        single_pmids = {p for row in single for p in str(row["pmids"]).split(";") if p}
        skeletal = [
            row
            for row in single
            if "skeletal/cartilage/joint" in str(row["phenotype_categories"]).split(";")
        ]
        evidence_rows.append(
            {
                "gene": gene,
                "single_gene_unique_mp_terms": len(single),
                "single_gene_publications": len(single_pmids),
                "single_gene_skeletal_cartilage_joint_terms": len(skeletal),
                "compound_genotype_unique_mp_terms": len(compound),
                "interpretation": (
                    "direct mouse phenotype annotations include skeletal/cartilage/joint terms"
                    if skeletal
                    else "mouse phenotype annotations exist, but no skeletal/cartilage/joint term was recovered in this release"
                    if single
                    else "no single-gene mouse phenotype annotation recovered in this release"
                ),
            }
        )
    return terms, category_rows, evidence_rows


def plot_heatmap(category_rows: list[dict[str, object]], output: Path) -> None:
    categories = list(CATEGORY_KEYWORDS) + ["other"]
    matrix = np.zeros((len(GENES), len(categories)), dtype=int)
    for row in category_rows:
        if row["evidence_scope"] == "single_gene":
            matrix[
                GENES.index(str(row["gene"])), categories.index(str(row["category"]))
            ] = int(row["unique_mp_terms"])

    fig, ax = plt.subplots(figsize=(13, 4.8))
    image = ax.imshow(matrix, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(categories)), labels=categories, rotation=35, ha="right")
    ax.set_yticks(range(len(GENES)), labels=GENES)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j,
                i,
                str(matrix[i, j]),
                ha="center",
                va="center",
                fontsize=9,
                color="white" if matrix[i, j] > matrix.max() * 0.55 else "black",
            )
    ax.set_title("MGI single-gene phenotype annotations for the seven SLRPs")
    ax.set_xlabel("Broad ontology-derived phenotype category (categories may overlap)")
    ax.set_ylabel("Gene")
    fig.colorbar(image, ax=ax, label="Unique MP terms")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-pheno", required=True, type=Path)
    parser.add_argument("--ontology", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--figure-dir", required=True, type=Path)
    args = parser.parse_args()

    names, parents = parse_obo(args.ontology)
    records = parse_mgi(args.gene_pheno, names, parents)
    terms, categories, evidence = aggregate(records)

    write_tsv(
        args.output_dir / "mgi_focal_gene_phenotype_records.tsv",
        records,
        list(records[0]),
    )
    write_tsv(
        args.output_dir / "mgi_focal_gene_phenotype_terms.tsv", terms, list(terms[0])
    )
    write_tsv(
        args.output_dir / "mgi_phenotype_category_summary.tsv",
        categories,
        list(categories[0]),
    )
    write_tsv(
        args.output_dir / "mgi_phenotype_evidence_summary.tsv",
        evidence,
        list(evidence[0]),
    )
    plot_heatmap(categories, args.figure_dir / "mgi_single_gene_phenotype_heatmap.png")
    print(
        f"Recovered {len(records)} focal records and {len(terms)} unique gene/term/scope combinations"
    )


if __name__ == "__main__":
    main()
