#!/usr/bin/env python3
"""
reads the BGN/DCN tree, it assigns tips to each expected gene group then checks if
each group is seperated, the output lists the placements with conflicting labels
"""

from __future__ import annotations

import csv
from pathlib import Path

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[4]
TREE = (
    ROOT
    / "analyses/phylogenetics/combined_trees/seven_gene_tree/SLRP_selected_working_genes.treefile"
)
OUT = (
    ROOT
    / "analyses/phylogenetics/compact_panel/tables/class_I_BGN_DCN_tree_discrimination.tsv"
)


def parts(name: str) -> tuple[str, str, str]:
    # Decode canonical gene, species and accession fields from a tree tip.
    gene, species, accession = name.split("|", 2)
    return gene, species, accession


def main() -> None:
    tree = Phylo.read(TREE, "newick")
    tips = [tip for tip in tree.get_terminals() if parts(tip.name)[0] in {"BGN", "DCN"}]
    rows: list[dict[str, object]] = []
    for query in tips:
        gene, species, accession = parts(query.name)
        paralog = "DCN" if gene == "BGN" else "BGN"
        same = [tip for tip in tips if tip is not query and parts(tip.name)[0] == gene]
        other = [tip for tip in tips if parts(tip.name)[0] == paralog]
        nearest_same = min(same, key=lambda tip: tree.distance(query, tip))
        nearest_other = min(other, key=lambda tip: tree.distance(query, tip))
        same_distance = tree.distance(query, nearest_same)
        other_distance = tree.distance(query, nearest_other)
        margin = other_distance - same_distance
        rows.append(
            {
                "gene": gene,
                "species": species,
                "accession": accession,
                "nearest_same_gene_tip": nearest_same.name,
                "nearest_same_gene_distance": round(same_distance, 6),
                "nearest_paralog_tip": nearest_other.name,
                "nearest_paralog_distance": round(other_distance, 6),
                "paralog_minus_same_distance": round(margin, 6),
                "nearest_neighbor_support": (
                    "same_gene_closer" if margin > 0 else "paralog_closer_or_tied"
                ),
            }
        )

    rows.sort(key=lambda row: (str(row["gene"]), str(row["species"])))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT} ({len(rows)} class-I proteins)")
    for gene in ("BGN", "DCN"):
        subset = [row for row in rows if row["gene"] == gene]
        supported = sum(
            row["nearest_neighbor_support"] == "same_gene_closer" for row in subset
        )
        print(
            f"{gene}: {supported}/{len(subset)} have a closer same-gene than paralog tip"
        )


if __name__ == "__main__":
    main()
