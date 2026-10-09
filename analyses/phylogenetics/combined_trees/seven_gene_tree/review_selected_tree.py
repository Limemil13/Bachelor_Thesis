#!/usr/bin/env python3
"""reads the unrooted combined tree says if each focal gene is forming a cluster
or a split and reports nearests neighbours for flagged sequences"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import median

from Bio import Phylo

GENES = ("BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM")
EXPECTED_TIP_COUNT = 103


def gene_of(name: str) -> str:
    return name.split("|", 1)[0]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    # First verify that the tree contains the complete canonical panel. Then test
    # each gene as an unrooted split and inspect the nearest patristic neighbours of
    # candidates flagged elsewhere. This is a consistency check, not a new taxonomic
    # assignment based only on the closest tip.
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--protein-review", required=True, type=Path)
    parser.add_argument("--protspace-review", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    terminals = tree.get_terminals()
    by_name = {tip.name: tip for tip in terminals}
    if len(terminals) != EXPECTED_TIP_COUNT or len(by_name) != EXPECTED_TIP_COUNT:
        raise ValueError(
            f"Expected the final {EXPECTED_TIP_COUNT}-tip canonical compact tree"
        )

    gene_rows: list[dict[str, object]] = []
    all_names = set(by_name)
    for gene in GENES:
        tips = [tip for tip in terminals if gene_of(tip.name) == gene]
        target_names = {tip.name for tip in tips}
        ancestor = tree.common_ancestor(tips)
        descendants = ancestor.get_terminals()
        descendant_genes = sorted({gene_of(tip.name) for tip in descendants})
        matching_edge = None
        matching_side = ""
        largest_pure_group: set[str] = set()
        largest_pure_edge = None
        largest_pure_side = ""
        '''since its an unrooted tree, the gene could be on one or another side of an edge,
         so both sides are checked by looking into descendant and complement set
        '''
        for clade in tree.find_clades(order="preorder"):
            if clade is tree.root:
                continue
            side = {tip.name for tip in clade.get_terminals()}
            for side_name, candidate_side in (
                ("descendant", side),
                ("complement", all_names - side),
            ):
                if candidate_side <= target_names and len(candidate_side) > len(
                    largest_pure_group
                ):
                    largest_pure_group = candidate_side
                    largest_pure_edge = clade
                    largest_pure_side = side_name
            if side == target_names:
                matching_edge, matching_side = clade, "descendant"
                break
            if all_names - side == target_names:
                matching_edge, matching_side = clade, "complement"
                break
        within = [
            tree.distance(left, right)
            for i, left in enumerate(tips)
            for right in tips[i + 1 :]
        ]
        parent = next(
            (clade for clade in tree.find_clades() if ancestor in clade.clades), None
        )
        sister_genes: set[str] = set()
        if parent is not None:
            for child in parent.clades:
                if child is not ancestor:
                    sister_genes.update(
                        gene_of(tip.name) for tip in child.get_terminals()
                    )
        gene_rows.append(
            {
                "gene": gene,
                "tip_count": len(tips),
                "unrooted_monophyletic": "yes" if matching_edge is not None else "no",
                "supporting_split_side": matching_side,
                "supporting_split_support": (
                    matching_edge.confidence
                    if matching_edge is not None
                    and matching_edge.confidence is not None
                    else ""
                ),
                "supporting_split_label": matching_edge.name
                if matching_edge is not None
                else "",
                "largest_pure_gene_split_tip_count": len(largest_pure_group),
                "largest_pure_gene_split_side": largest_pure_side,
                "largest_pure_gene_split_support": (
                    largest_pure_edge.confidence
                    if largest_pure_edge is not None
                    and largest_pure_edge.confidence is not None
                    else ""
                ),
                "largest_pure_gene_split_label": (
                    largest_pure_edge.name if largest_pure_edge is not None else ""
                ),
                "tips_outside_largest_pure_gene_split": ",".join(
                    sorted(target_names - largest_pure_group)
                ),
                "mrca_descendant_tip_count": len(descendants),
                "mrca_descendant_genes": ",".join(descendant_genes),
                "mrca_support": ancestor.confidence
                if ancestor.confidence is not None
                else "",
                "mrca_label": ancestor.name or "",
                "sister_clade_genes": ",".join(sorted(sister_genes)),
                "median_within_gene_patristic_distance": round(median(within), 6),
                "maximum_within_gene_patristic_distance": round(max(within), 6),
            }
        )
    write_tsv(args.out_dir / "compact_tree_gene_clade_review.tsv", gene_rows)
    flagged: set[str] = set()
    with args.protein_review.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Protein QC status"] == "Needs inspection":
                flagged.add(
                    f"{row['Gene']}|{row['Species']}|{row['Protein accession']}"
                )
    with args.protspace_review.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["embedding_qc_status"] == "Needs inspection":
                species, gene, accession = row["species_accession"].split("|", 2)
                flagged.add(f"{gene}|{species}|{accession}")

    neighbor_rows: list[dict[str, object]] = []
    for query_name in sorted(flagged):
        query = by_name.get(query_name)
        if query is None:
            raise KeyError(f"Flagged protein is absent from tree: {query_name}")
        neighbors = sorted(
            (
                (tree.distance(query, tip), tip.name)
                for tip in terminals
                if tip is not query
            ),
            key=lambda item: (item[0], item[1]),
        )[:5]
        for rank, (distance, neighbor_name) in enumerate(neighbors, start=1):
            neighbor_rows.append(
                {
                    "query_tip": query_name,
                    "query_gene": gene_of(query_name),
                    "neighbor_rank": rank,
                    "neighbor_tip": neighbor_name,
                    "neighbor_gene": gene_of(neighbor_name),
                    "same_gene": "yes"
                    if gene_of(query_name) == gene_of(neighbor_name)
                    else "no",
                    "patristic_distance": round(distance, 6),
                }
            )
    write_tsv(args.out_dir / "compact_tree_flagged_tip_neighbors.tsv", neighbor_rows)

    print(f"Reviewed {len(terminals)} tips across {len(gene_rows)} genes")
    print(
        f"Unrooted-monophyletic genes: {sum(row['unrooted_monophyletic'] == 'yes' for row in gene_rows)}/{len(gene_rows)}"
    )
    print(f"Flagged tips with nearest-neighbor evidence: {len(flagged)}")


if __name__ == "__main__":
    main()
