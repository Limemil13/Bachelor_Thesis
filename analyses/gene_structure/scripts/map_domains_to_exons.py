#!/usr/bin/env python3
"""Parse representative-protein Pfam hits and map them onto coding exons."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parents[3]
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN", "OMD")
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
GENE_COLORS = {
    "BGN": "#4C78A8",
    "DCN": "#9C6ADE",
    "FMOD": "#F58518",
    "PRELP": "#ECA82C",
    "EPYC": "#B279A2",
    "LUM": "#54A24B",
    "OGN": "#E45756",
    "OMD": "#8C8C8C",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_domtblout(path: Path, threshold: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split(maxsplit=22)
            if len(parts) < 23:
                continue
            i_evalue = float(parts[12])
            if i_evalue > threshold:
                continue
            query = parts[3]
            gene, species, protein_accession = query.split("|", 2)
            name, accession = parts[0], parts[1]
            description = parts[22]
            is_lrr = (
                "lrr" in name.lower()
                or "leucine rich repeat" in description.lower()
                or "leucine-rich repeat" in description.lower()
                or name.lower().startswith("slrp")
            )
            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "protein_accession": protein_accession,
                    "protein_length_aa": int(parts[5]),
                    "domain_name": name,
                    "domain_accession": accession,
                    "description": description,
                    "i_evalue": i_evalue,
                    "domain_score": float(parts[13]),
                    "ali_from_aa": int(parts[17]),
                    "ali_to_aa": int(parts[18]),
                    "env_from_aa": int(parts[19]),
                    "env_to_aa": int(parts[20]),
                    "is_lrr_related": "yes" if is_lrr else "no",
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domtblout", type=Path, required=True)
    parser.add_argument(
        "--blocks",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/extended/tables/coding_exon_splice_phase_5species.tsv",
    )
    parser.add_argument(
        "--qc",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/extended/tables/full_gene_structure_qc.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "analyses/gene_structure/extended",
    )
    parser.add_argument("--i-evalue", type=float, default=1e-3)
    args = parser.parse_args()

    hits = parse_domtblout(args.domtblout, args.i_evalue)
    blocks = read_tsv(args.blocks)
    qc = read_tsv(args.qc)
    blocks_by_key: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for block in blocks:
        blocks_by_key[(block["gene"], block["species"])].append(block)

    overlap_rows: list[dict[str, object]] = []
    for hit_index, hit in enumerate(hits, 1):
        domain_nt_start = (int(hit["ali_from_aa"]) - 1) * 3 + 1
        domain_nt_end = int(hit["ali_to_aa"]) * 3
        touched: list[int] = []
        for block in blocks_by_key[(str(hit["gene"]), str(hit["species"]))]:
            block_start = int(block["cumulative_cds_start_bp"])
            block_end = int(block["cumulative_cds_end_bp"])
            overlap_start = max(domain_nt_start, block_start)
            overlap_end = min(domain_nt_end, block_end)
            if overlap_start > overlap_end:
                continue
            touched.append(int(block["exon_index"]))
            overlap_rows.append(
                {
                    "domain_hit_id": hit_index,
                    **hit,
                    "cds_index": block["cds_index"],
                    "exon_index": block["exon_index"],
                    "overlap_start_aa": math.ceil(overlap_start / 3),
                    "overlap_end_aa": math.ceil(overlap_end / 3),
                    "overlap_length_aa": math.ceil(overlap_end / 3)
                    - math.ceil(overlap_start / 3)
                    + 1,
                }
            )
        hit["coding_exons_spanned"] = ",".join(map(str, sorted(set(touched))))
        hit["coding_exon_span_count"] = len(set(touched))

    tables = args.out_dir / "tables"
    hit_fields = [
        "gene",
        "species",
        "protein_accession",
        "protein_length_aa",
        "domain_name",
        "domain_accession",
        "description",
        "i_evalue",
        "domain_score",
        "ali_from_aa",
        "ali_to_aa",
        "env_from_aa",
        "env_to_aa",
        "is_lrr_related",
        "coding_exons_spanned",
        "coding_exon_span_count",
    ]
    write_tsv(tables / "representative_pfam_domain_hits.tsv", hits, hit_fields)
    write_tsv(
        tables / "representative_domain_exon_overlap.tsv",
        overlap_rows,
        list(overlap_rows[0]),
    )

    summary_rows: list[dict[str, object]] = []
    hit_df = pd.DataFrame(hits)
    qc_df = pd.DataFrame(qc)
    qc_df["translated_length_aa"] = pd.to_numeric(qc_df["translated_length_aa"])
    for gene in GENES:
        gene_hits = hit_df[hit_df["gene"].eq(gene)]
        counts = gene_hits.groupby("species").size().reindex(SPECIES, fill_value=0)
        lrr_counts = (
            gene_hits[gene_hits["is_lrr_related"].eq("yes")]
            .groupby("species")
            .size()
            .reindex(SPECIES, fill_value=0)
        )
        spanned = (
            gene_hits[gene_hits["is_lrr_related"].eq("yes")]
            .groupby("species")["coding_exon_span_count"]
            .max()
            .reindex(SPECIES, fill_value=0)
        )
        summary_rows.append(
            {
                "gene": gene,
                "species_count": qc_df[qc_df["gene"].eq(gene)]["species"].nunique(),
                "significant_pfam_hit_range": f"{counts.min()}-{counts.max()}",
                "lrr_related_hit_range": f"{lrr_counts.min()}-{lrr_counts.max()}",
                "max_coding_exons_spanned_by_one_lrr_hit_range": f"{spanned.min()}-{spanned.max()}",
                "all_species_have_lrr_support": "yes"
                if (lrr_counts > 0).all()
                else "no",
                "interpretation": (
                    "LRR-family architecture is supported in all five representative translations"
                    if (lrr_counts > 0).all()
                    else "one or more representative translations lack significant LRR-family support"
                ),
            }
        )
    write_tsv(
        tables / "representative_domain_exon_gene_summary.tsv",
        summary_rows,
        list(summary_rows[0]),
    )

    # Human overview: coding-exon positions provide the structural frame and
    # significant LRR-family hits show how the protein architecture crosses it.
    human_blocks = pd.DataFrame(blocks)
    human_blocks = human_blocks[human_blocks["species"].eq("human")]
    human_hits = hit_df[
        (hit_df["species"].eq("human")) & hit_df["is_lrr_related"].eq("yes")
    ]
    human_qc = qc_df[qc_df["species"].eq("human")].set_index("gene")
    fig, ax = plt.subplots(figsize=(14, 8))
    for row_index, gene in enumerate(reversed(GENES)):
        y = row_index
        gene_blocks = human_blocks[human_blocks["gene"].eq(gene)].sort_values(
            "cds_index"
        )
        length = float(human_qc.loc[gene, "translated_length_aa"])
        ax.hlines(y, 0, length, color="#C6C6C6", linewidth=12, zorder=1)
        for _, block in gene_blocks.iterrows():
            x0 = (float(block["cumulative_cds_start_bp"]) - 1) / 3
            x1 = float(block["cumulative_cds_end_bp"]) / 3
            alpha = 0.50 if int(block["cds_index"]) % 2 else 0.85
            ax.add_patch(
                Rectangle(
                    (x0, y - 0.12),
                    x1 - x0,
                    0.24,
                    color=GENE_COLORS[gene],
                    alpha=alpha,
                    zorder=2,
                )
            )
        levels: list[float] = []
        for _, hit in (
            human_hits[human_hits["gene"].eq(gene)]
            .sort_values("ali_from_aa")
            .iterrows()
        ):
            start, end = float(hit["ali_from_aa"]), float(hit["ali_to_aa"])
            level = 0
            while level < len(levels) and start <= levels[level] + 3:
                level += 1
            if level == len(levels):
                levels.append(end)
            else:
                levels[level] = end
            yy = y + 0.19 + level * 0.10
            ax.add_patch(
                Rectangle(
                    (start, yy),
                    end - start + 1,
                    0.065,
                    facecolor="#6F4E9C",
                    edgecolor="none",
                    zorder=3,
                )
            )
    ax.set_yticks(range(len(GENES)), list(reversed(GENES)))
    ax.set_xlim(0, max(human_qc["translated_length_aa"]) * 1.04)
    ax.set_ylim(-0.6, len(GENES) - 0.1)
    ax.set_xlabel("Position in human representative protein (amino acids)")
    ax.set_title(
        "LRR-family Pfam hits mapped onto representative coding exons",
        fontweight="bold",
    )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", alpha=0.15)
    ax.legend(
        handles=[
            Patch(color="#8F8F8F", label="Coding exons (alternating opacity)"),
            Patch(color="#6F4E9C", label="Significant LRR-family Pfam hit"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.11),
        ncol=2,
        frameon=False,
    )
    figures = args.out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    base = figures / "human_domain_exon_architecture"
    fig.savefig(
        base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Parsed {len(hits)} significant representative-protein Pfam hits")
    print(f"Mapped them to {len(overlap_rows)} domain-by-coding-exon overlaps")


if __name__ == "__main__":
    main()
