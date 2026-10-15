#!/usr/bin/env python3
"""Root a targeted tree and write compact clade and focus-record summaries."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from Bio import Phylo


def tip_group(name: str) -> str:
    return name.split("|", 1)[0]


def node_support(clade: object) -> str:
    name = getattr(clade, "name", None)
    confidence = getattr(clade, "confidence", None)
    if name not in {None, ""}:
        return str(name)
    return "" if confidence is None else str(confidence)


def split_support(value: str) -> tuple[str, str]:
    if "/" in value:
        return tuple(value.split("/", 1))  # type: ignore[return-value]
    return value, ""


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def root_tree(tree: object, outgroup_groups: set[str]) -> str:
    outgroup = [tip for tip in tree.get_terminals() if tip_group(tip.name) in outgroup_groups]
    if not outgroup:
        raise RuntimeError(f"No outgroup tips found for {sorted(outgroup_groups)}")
    mrca = tree.common_ancestor(outgroup)
    descendants = set(mrca.get_terminals())
    if descendants != set(outgroup):
        candidates: list[tuple[int, str, object]] = []
        for group in sorted(outgroup_groups):
            tips = [tip for tip in outgroup if tip_group(tip.name) == group]
            if not tips:
                continue
            group_mrca = tree.common_ancestor(tips)
            if set(group_mrca.get_terminals()) == set(tips):
                candidates.append((len(tips), group, group_mrca))
        if not candidates:
            raise RuntimeError("The combined outgroup and its individual groups are not monophyletic")
        _, group, mrca = max(candidates, key=lambda item: item[0])
        basis = f"monophyletic {group} group"
    else:
        basis = "combined " + ", ".join(sorted(outgroup_groups)) + " outgroup"
    branch = (mrca.branch_length or 0.0) / 2.0
    tree.root_with_outgroup(mrca, outgroup_branch_length=branch)
    return basis


def summarize_groups(tree: object) -> list[dict[str, str]]:
    groups = sorted({tip_group(tip.name) for tip in tree.get_terminals()})
    rows: list[dict[str, str]] = []
    for group in groups:
        tips = [tip for tip in tree.get_terminals() if tip_group(tip.name) == group]
        mrca = tree.common_ancestor(tips)
        descendants = mrca.get_terminals()
        same = [tip for tip in descendants if tip_group(tip.name) == group]
        support = node_support(mrca)
        alrt, ufboot = split_support(support)
        rows.append(
            {
                "group": group,
                "tips": str(len(tips)),
                "mrca_descendants": str(len(descendants)),
                "purity": f"{len(same) / len(descendants):.3f}",
                "monophyletic": "yes" if set(descendants) == set(tips) else "no",
                "support_raw": support,
                "sh_alrt": alrt,
                "ufboot": ufboot,
                "unexpected_groups": ",".join(
                    sorted({tip_group(tip.name) for tip in descendants if tip_group(tip.name) != group})
                ),
            }
        )
    return rows


def summarize_focus(tree: object, focus_groups: set[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    terminals = tree.get_terminals()
    for tip in terminals:
        if tip_group(tip.name) not in focus_groups:
            continue
        neighbours = sorted(
            ((tree.distance(tip, other), other) for other in terminals if other is not tip),
            key=lambda item: item[0],
        )
        for rank, (distance, neighbour) in enumerate(neighbours[:5], start=1):
            rows.append(
                {
                    "focus_tip": tip.name,
                    "rank": str(rank),
                    "neighbour_tip": neighbour.name,
                    "neighbour_group": tip_group(neighbour.name),
                    "patristic_distance": f"{distance:.8f}",
                }
            )
    return rows


def subset_row(tree: object, label: str, groups: set[str]) -> dict[str, str]:
    selected = [tip for tip in tree.get_terminals() if tip_group(tip.name) in groups]
    mrca = tree.common_ancestor(selected)
    descendants = mrca.get_terminals()
    unexpected = [tip for tip in descendants if tip not in selected]
    support = node_support(mrca)
    alrt, ufboot = split_support(support)
    return {
        "subset": label,
        "groups": ",".join(sorted(groups)),
        "selected_tips": str(len(selected)),
        "mrca_descendants": str(len(descendants)),
        "monophyletic": "yes" if not unexpected else "no",
        "support_raw": support,
        "sh_alrt": alrt,
        "ufboot": ufboot,
        "unexpected_groups": ",".join(sorted({tip_group(tip.name) for tip in unexpected})),
    }


def export_layout(tree: object, output_dir: Path) -> None:
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    y_position = {tip: float(index) for index, tip in enumerate(tree.get_terminals())}

    def assign_y(clade: object) -> float:
        if clade in y_position:
            return y_position[clade]
        values = [assign_y(child) for child in clade.clades]
        y_position[clade] = sum(values) / len(values)
        return y_position[clade]

    assign_y(tree.root)
    node_ids: dict[object, str] = {}
    node_rows: list[dict[str, str]] = []
    edge_rows: list[dict[str, str]] = []
    for index, clade in enumerate(tree.find_clades(order="preorder")):
        node_id = f"n{index}"
        node_ids[clade] = node_id
        terminal = clade.is_terminal()
        label = clade.name if terminal else ""
        parts = label.split("|") if label else []
        node_rows.append(
            {
                "node_id": node_id,
                "node_type": "tip" if terminal else "internal",
                "x": f"{depths[clade]:.9f}",
                "y": f"{y_position[clade]:.9f}",
                "label": label,
                "group": parts[0] if parts else "",
                "support": "" if terminal else node_support(clade),
            }
        )
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            edge_rows.append(
                {
                    "parent_id": node_ids[parent],
                    "child_id": node_ids[child],
                    "parent_x": f"{depths[parent]:.9f}",
                    "parent_y": f"{y_position[parent]:.9f}",
                    "child_x": f"{depths[child]:.9f}",
                    "child_y": f"{y_position[child]:.9f}",
                }
            )
    write_tsv(
        output_dir / "tree_nodes.tsv",
        ["node_id", "node_type", "x", "y", "label", "group", "support"],
        node_rows,
    )
    write_tsv(
        output_dir / "tree_edges.tsv",
        ["parent_id", "child_id", "parent_x", "parent_y", "child_x", "child_y"],
        edge_rows,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", choices=("class1", "class2"), required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.tree, args.output_dir / "unrooted.treefile")
    tree = Phylo.read(args.tree, "newick")

    if args.analysis == "class1":
        outgroups = {"FMOD", "PRELP"}
        focus = {"BGN", "DCN2", "SLRP1", "SLRP2", "SLRP3", "SLRP4"}
        subsets = (
            ("BGN_ASPN", {"BGN", "ASPN"}),
            ("DCN1_DCN2", {"DCN", "DCN2"}),
            ("CLASS1_TARGET", {"BGN", "ASPN", "DCN", "DCN2", "SLRP1", "SLRP2", "SLRP3", "SLRP4"}),
        )
    else:
        outgroups = {"KERA", "OMD", "PRELP"}
        focus = {"FMODB_CHECK", "LUM_CHECK", "SLRP6", "SLRP7", "SLRP8"}
        subsets = (
            ("FMOD_ALL", {"FMOD", "FMODA", "FMODB", "FMODB_CHECK"}),
            ("LUM_ALL", {"LUM", "LUM_CHECK"}),
            ("CLASS2A_TARGET", {"FMOD", "FMODA", "FMODB", "FMODB_CHECK", "LUM", "LUM_CHECK", "LUML", "SLRP6", "SLRP7", "SLRP8"}),
        )

    basis = root_tree(tree, outgroups)
    Phylo.write(tree, args.output_dir / "rooted.treefile", "newick")
    (args.output_dir / "rooting_note.txt").write_text(f"Rooted using {basis}.\n", encoding="utf-8")
    write_tsv(
        args.output_dir / "group_summary.tsv",
        ["group", "tips", "mrca_descendants", "purity", "monophyletic", "support_raw", "sh_alrt", "ufboot", "unexpected_groups"],
        summarize_groups(tree),
    )
    write_tsv(
        args.output_dir / "focus_nearest_neighbours.tsv",
        ["focus_tip", "rank", "neighbour_tip", "neighbour_group", "patristic_distance"],
        summarize_focus(tree, focus),
    )
    write_tsv(
        args.output_dir / "topology_summary.tsv",
        ["subset", "groups", "selected_tips", "mrca_descendants", "monophyletic", "support_raw", "sh_alrt", "ufboot", "unexpected_groups"],
        [subset_row(tree, label, groups) for label, groups in subsets],
    )
    export_layout(tree, args.output_dir)
    print(f"{args.analysis}: rooted tree and summaries written; {basis}")


if __name__ == "__main__":
    main()
