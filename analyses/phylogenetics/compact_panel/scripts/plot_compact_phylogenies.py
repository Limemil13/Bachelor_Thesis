#!/usr/bin/env python3
"""Render the valid compact canonical trees as static PNG/SVG figures.

The trees are unrooted maximum-likelihood topologies.  Rectangular rendering
requires a display root, but that orientation is arbitrary and must not be
interpreted as the direction of evolution.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[4]
GENES = ["BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN"]
GENE_COLORS = {
    "BGN": "#4C78A8",
    "DCN": "#9C6ADE",
    "FMOD": "#F58518",
    "PRELP": "#ECA82C",
    "EPYC": "#B279A2",
    "LUM": "#54A24B",
    "OGN": "#E45756",
}
TAXON_COLORS = {
    "human": "#B2182B",
    "mouse": "#B2182B",
    "rattus_norvegicus": "#B2182B",
    "canis_lupus_familiaris": "#B2182B",
    "bos_taurus": "#B2182B",
    "monodelphis_domestica": "#EF8A62",
    "chicken": "#D8A700",
    "anolis_carolinensis": "#66A61E",
    "frog": "#1B9E77",
    "zebrafish": "#2166AC",
    "lepisosteus_oculatus": "#2166AC",
    "latimeria_chalumnae": "#4393C3",
    "callorhinchus_milii": "#762A83",
    "scyliorhinus_canicula": "#762A83",
    "whale_shark": "#762A83",
    "lamprey": "#A6D854",
    "amphioxus": "#666666",
}
DISPLAY_SPECIES = {
    "rattus_norvegicus": "Rat",
    "canis_lupus_familiaris": "Dog",
    "bos_taurus": "Cow",
    "monodelphis_domestica": "Opossum",
    "anolis_carolinensis": "Anole",
    "lepisosteus_oculatus": "Spotted gar",
    "latimeria_chalumnae": "Coelacanth",
    "callorhinchus_milii": "Elephant shark",
    "scyliorhinus_canicula": "Catshark",
    "whale_shark": "Whale shark",
}


def tip_parts(name: str) -> tuple[str, str, str]:
    parts = name.split("|")
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    return name, "", ""


def tip_label(clade) -> str | None:
    if not clade.is_terminal():
        return None
    species, _gene, accession = tip_parts(clade.name)
    display = DISPLAY_SPECIES.get(species, species.replace("_", " ").title())
    marker = "  *" if accession == "XP_075930353.1" else ""
    return f"{display}  {accession}{marker}"


def support_label(clade) -> str | None:
    if clade.is_terminal():
        return None
    value = clade.name
    if value and "/" in str(value):
        return str(value)
    if clade.confidence is not None:
        return f"{clade.confidence:g}"
    return None


def model_from_report(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"Best-fit model according to BIC:\s*(\S+)", text)
    return match.group(1) if match else "model not parsed"


def save(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_compact_trees() -> None:
    fig, axes = plt.subplots(4, 2, figsize=(22, 30))
    for gene, ax in zip(GENES, axes.flat, strict=False):
        tree_path = (
            ROOT / f"analyses/phylogenetics/compact_panel/trees/{gene}/{gene}.treefile"
        )
        report_path = tree_path.with_suffix(".iqtree")
        tree = Phylo.read(tree_path, "newick")
        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            label_func=tip_label,
            branch_labels=support_label,
            show_confidence=False,
        )
        terminal_color_by_label = {}
        for terminal in tree.get_terminals():
            species, _gene, _accession = tip_parts(terminal.name)
            terminal_color_by_label[tip_label(terminal)] = TAXON_COLORS.get(
                species, "#333333"
            )
        for text in ax.texts:
            rendered = text.get_text().strip()
            if rendered in terminal_color_by_label:
                text.set_color(terminal_color_by_label[rendered])
                text.set_fontsize(8)
            elif "/" in rendered:
                text.set_fontsize(6)
                text.set_color("#555555")
        ax.set_title(
            f"{gene}: {len(tree.get_terminals())} proteins, {model_from_report(report_path)}",
            fontsize=13,
            color=GENE_COLORS[gene],
            fontweight="bold",
        )
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.set_xlabel("Substitutions per site", fontsize=9)
        ax.grid(axis="x", alpha=0.2)
    for ax in axes.flat[len(GENES) :]:
        ax.set_visible(False)
    fig.suptitle(
        "Compact canonical per-gene maximum-likelihood trees",
        fontsize=20,
        y=1.005,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.002,
        "Unrooted topologies; rectangular display roots are arbitrary. Internal labels are SH-aLRT/UFBoot support. * = deep-lineage assignment requiring caveat.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout()
    save(fig, ROOT / "analyses/phylogenetics/figures/compact_per_gene_trees")


def combined_label(clade) -> str | None:
    if not clade.is_terminal():
        return None
    parts = clade.name.split("|")
    gene, species, accession = parts[:3] if len(parts) >= 3 else ("", clade.name, "")
    display = DISPLAY_SPECIES.get(species, species.replace("_", " ").title())
    marker = "  *" if accession == "XP_075930353.1" else ""
    return f"{gene} | {display} | {accession}{marker}"


def plot_combined_tree() -> None:
    tree_path = (
        ROOT
        / "analyses/phylogenetics/combined_trees/six_gene_tree/SLRP_selected_working_genes.treefile"
    )
    report_path = tree_path.with_suffix(".iqtree")
    tree = Phylo.read(tree_path, "newick")
    fig, ax = plt.subplots(figsize=(18, 30))
    Phylo.draw(
        tree, axes=ax, do_show=False, label_func=combined_label, show_confidence=False
    )
    terminal_color_by_label = {}
    for terminal in tree.get_terminals():
        parts = terminal.name.split("|")
        gene = parts[0] if len(parts) >= 3 else ""
        terminal_color_by_label[combined_label(terminal)] = GENE_COLORS.get(
            gene, "#333333"
        )
    for text in ax.texts:
        rendered = text.get_text().strip()
        if rendered in terminal_color_by_label:
            text.set_color(terminal_color_by_label[rendered])
            text.set_fontsize(7)
        else:
            text.set_visible(False)
    ax.set_title(
        f"Seven-gene canonical SLRP tree — {len(tree.get_terminals())} proteins, {model_from_report(report_path)}",
        fontsize=18,
        fontweight="bold",
    )
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_xlabel("Substitutions per site")
    ax.grid(axis="x", alpha=0.2)
    fig.text(
        0.5,
        0.005,
        "Unrooted maximum-likelihood topology; display orientation is arbitrary. * = tentative lamprey LUM-like assignment.",
        ha="center",
        fontsize=10,
    )
    save(fig, ROOT / "analyses/phylogenetics/figures/six_gene_combined_tree")


def main() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plot_compact_trees()
    plot_combined_tree()
    print("Rendered compact per-gene and combined tree figures.")


if __name__ == "__main__":
    main()
