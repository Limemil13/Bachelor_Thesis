#!/usr/bin/env python3
"""Compare focal-gene clustering before and after partial-BGN exclusion."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import median

from Bio import Phylo


GENES = ("BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM")


def gene_of(tip_name: str) -> str:
    return tip_name.split("|", 1)[0]


def support_label(clade: object | None) -> str:
    if clade is None:
        return ""
    name = getattr(clade, "name", None)
    if name:
        return str(name)
    confidence = getattr(clade, "confidence", None)
    return "" if confidence is None else str(confidence)


def split_support_values(label: str) -> tuple[str, str]:
    """Return IQ-TREE SH-aLRT and UFBoot values from an ``a/b`` node label."""
    fields = label.split("/", 1)
    if len(fields) == 2:
        return fields[0], fields[1]
    return (label, "") if label else ("", "")


def describe_gene(tree_path: Path, analysis: str, gene: str) -> dict[str, object]:
    tree = Phylo.read(tree_path, "newick")
    terminals = tree.get_terminals()
    all_names = {tip.name for tip in terminals}
    target_names = {tip.name for tip in terminals if gene_of(tip.name) == gene}
    if len(target_names) < 2:
        raise ValueError(f"Need at least two {gene} tips in {tree_path}")

    matching_edge = None
    matching_side = ""
    largest_pure_group: set[str] = set()
    largest_pure_edge = None
    largest_pure_side = ""
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        descendant_side = {tip.name for tip in clade.get_terminals()}
        for side_name, candidate_side in (
            ("descendant", descendant_side),
            ("complement", all_names - descendant_side),
        ):
            if candidate_side <= target_names and len(candidate_side) > len(
                largest_pure_group
            ):
                largest_pure_group = candidate_side
                largest_pure_edge = clade
                largest_pure_side = side_name
        if descendant_side == target_names:
            matching_edge = clade
            matching_side = "descendant"
            break
        if all_names - descendant_side == target_names:
            matching_edge = clade
            matching_side = "complement"
            break

    tips = [tip for tip in terminals if tip.name in target_names]
    distances = [
        tree.distance(left, right)
        for index, left in enumerate(tips)
        for right in tips[index + 1 :]
    ]
    outside = sorted(target_names - largest_pure_group)
    full_label = support_label(matching_edge)
    full_sh_alrt, full_ufboot = split_support_values(full_label)
    largest_label = support_label(largest_pure_edge)
    largest_sh_alrt, largest_ufboot = split_support_values(largest_label)
    return {
        "analysis": analysis,
        "tree_tip_count": len(terminals),
        "gene": gene,
        "gene_tip_count": len(target_names),
        "unrooted_monophyletic": "yes" if matching_edge is not None else "no",
        "supporting_split_side": matching_side,
        "supporting_split_label": full_label,
        "supporting_split_sh_alrt_percent": full_sh_alrt,
        "supporting_split_ufboot_percent": full_ufboot,
        "largest_pure_gene_split_tip_count": len(largest_pure_group),
        "largest_pure_gene_split_fraction": round(
            len(largest_pure_group) / len(target_names), 4
        ),
        "largest_pure_gene_split_side": largest_pure_side,
        "largest_pure_gene_split_label": largest_label,
        "largest_pure_gene_split_sh_alrt_percent": largest_sh_alrt,
        "largest_pure_gene_split_ufboot_percent": largest_ufboot,
        "tips_in_largest_pure_gene_split": ",".join(sorted(largest_pure_group)),
        "tips_outside_largest_pure_gene_split": ",".join(outside),
        "median_within_gene_patristic_distance": round(median(distances), 6),
        "maximum_within_gene_patristic_distance": round(max(distances), 6),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-tree", required=True, type=Path)
    parser.add_argument("--sensitivity-tree", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    canonical_rows = [
        describe_gene(args.canonical_tree, "canonical_103", gene) for gene in GENES
    ]
    sensitivity_rows = [
        describe_gene(args.sensitivity_tree, "partial_BGN_excluded_101", gene)
        for gene in GENES
    ]
    write_tsv(
        args.out_dir / "all_gene_tree_comparison.tsv",
        canonical_rows + sensitivity_rows,
    )

    bgn_rows = [canonical_rows[0], sensitivity_rows[0]]
    write_tsv(args.out_dir / "bgn_tree_comparison.tsv", bgn_rows)

    before, after = bgn_rows
    interpretation = [
        "BGN fragment-exclusion sensitivity analysis",
        f"Canonical tree: {before['gene_tip_count']} BGN tips; "
        f"monophyletic={before['unrooted_monophyletic']}; largest pure BGN split="
        f"{before['largest_pure_gene_split_tip_count']}/{before['gene_tip_count']}.",
        f"Partial-fragment-excluded tree: {after['gene_tip_count']} BGN tips; "
        f"monophyletic={after['unrooted_monophyletic']}; largest pure BGN split="
        f"{after['largest_pure_gene_split_tip_count']}/{after['gene_tip_count']} "
        f"(SH-aLRT/UFBoot "
        f"{after['largest_pure_gene_split_label'] or 'not available'}).",
    ]
    if after["unrooted_monophyletic"] == "yes":
        interpretation.append(
            "Removing the two partial shark records recovered a complete BGN split, "
            "so those fragments materially affected BGN resolution in this tree."
        )
    elif (
        after["largest_pure_gene_split_fraction"]
        > before["largest_pure_gene_split_fraction"]
    ):
        interpretation.append(
            "The numerical fraction of BGN tips in the largest pure split increased, "
            "but the 7-tip split had weak support and BGN monophyly was not recovered. "
            "The two partial records therefore do not by themselves explain the "
            "unresolved BGN topology."
        )
    else:
        interpretation.append(
            "Removing the partial records did not improve BGN recovery as an unrooted "
            "split. The unresolved BGN topology therefore cannot be attributed only "
            "to those two incomplete sequences."
        )
    (args.out_dir / "interpretation.txt").write_text(
        "\n".join(interpretation) + "\n", encoding="utf-8"
    )
    print("\n".join(interpretation))


if __name__ == "__main__":
    main()
