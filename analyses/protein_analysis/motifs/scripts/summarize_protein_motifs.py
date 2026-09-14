#!/usr/bin/env python3
"""Summarize MEME/STREME protein motif models and seven-way MAST scans."""

from __future__ import annotations

import argparse
import csv
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
HIGHLIGHTS = {
    "XP_048463897.1": "EPYC whale shark (manual tentative)",
    "XP_006631165.2": "OGN spotted gar (manual tentative)",
    "XP_056660002.1": "OGN opossum (manual retain)",
    "NP_001013588.1": "OGN zebrafish (manual retain)",
    "XP_038638277.1": "BGN catshark (partial-model review)",
    "XP_048475905.1": "BGN whale shark (partial-model review)",
    "XP_001363160.3": "DCN opossum (partial-model review)",
}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def meme_motifs(path: Path, scope: str) -> list[dict[str, object]]:
    root = ET.parse(path).getroot()
    container = root.find("motifs")
    if container is None:
        raise ValueError(f"No MEME motif container in {path}")
    rows: list[dict[str, object]] = []
    for motif in container.findall("motif"):
        rows.append(
            {
                "scope": scope,
                "motif": motif.attrib.get("alt", motif.attrib["id"]),
                "consensus": motif.attrib.get("name", ""),
                "width_aa": motif.attrib.get("width", ""),
                "sites": motif.attrib.get("sites", ""),
                "information_content_bits": motif.attrib.get("ic", ""),
                "p_value": motif.attrib.get("p_value", ""),
                "e_value": motif.attrib.get("e_value", ""),
            }
        )
    return rows


def streme_motifs(path: Path, gene: str) -> list[dict[str, object]]:
    root = ET.parse(path).getroot()
    rows: list[dict[str, object]] = []
    for motif in root.findall(".//motifs/motif"):
        motif_id = motif.attrib["id"]
        rows.append(
            {
                "gene": gene,
                "motif": motif.attrib.get("alt", motif_id),
                "consensus": motif_id.split("-", 1)[1] if "-" in motif_id else motif_id,
                "width_aa": motif.attrib.get("width", ""),
                "positive_sequences": motif.attrib.get("train_pos_count", ""),
                "negative_sequences": motif.attrib.get("train_neg_count", ""),
                "training_p_value": motif.attrib.get("train_pvalue", ""),
                "holdout_sequences": int(motif.attrib.get("test_pos_count", "0"))
                + int(motif.attrib.get("test_neg_count", "0")),
                "holdout_p_value": motif.attrib.get("test_pvalue", ""),
                "validation_status": "training-only; no holdout available",
            }
        )
    return rows


def safe_neg_log10(value: str, cap: float = 300.0) -> float:
    number = float(value)
    if number <= 0:
        return cap
    return min(-math.log10(number), cap)


def mast_scores(path: Path, model_gene: str) -> list[dict[str, object]]:
    root = ET.parse(path).getroot()
    rows: list[dict[str, object]] = []
    for sequence in root.findall(".//sequences/sequence"):
        score = sequence.find("score")
        if score is None:
            continue
        true_gene, species, accession = sequence.attrib["name"].split("|", 2)
        hits = sequence.findall(".//hit")
        hit_score = sum(safe_neg_log10(hit.attrib["pvalue"], cap=100.0) for hit in hits)
        rows.append(
            {
                "true_gene": true_gene,
                "species": species,
                "accession": accession,
                "model_gene": model_gene,
                "sequence_length_aa": sequence.attrib["length"],
                "mast_e_value": score.attrib["evalue"],
                "mast_neg_log10_e_value_capped_300": safe_neg_log10(
                    score.attrib["evalue"]
                ),
                "significant_motif_hits": len(hits),
                "motif_hit_score_sum_neg_log10_p": hit_score,
            }
        )
    return rows


def candidate_assignments(scores: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in scores:
        grouped[
            (str(row["true_gene"]), str(row["species"]), str(row["accession"]))
        ].append(row)
    result: list[dict[str, object]] = []
    for (true_gene, species, accession), rows in grouped.items():
        ordered = sorted(
            rows,
            key=lambda row: float(row["motif_hit_score_sum_neg_log10_p"]),
            reverse=True,
        )
        own = next(row for row in rows if row["model_gene"] == true_gene)
        best_other = max(
            (row for row in rows if row["model_gene"] != true_gene),
            key=lambda row: float(row["motif_hit_score_sum_neg_log10_p"]),
        )
        own_score = float(own["motif_hit_score_sum_neg_log10_p"])
        other_score = float(best_other["motif_hit_score_sum_neg_log10_p"])
        result.append(
            {
                "gene": true_gene,
                "species": species,
                "accession": accession,
                "best_motif_model": ordered[0]["model_gene"],
                "own_model_rank": 1
                + next(
                    i for i, row in enumerate(ordered) if row["model_gene"] == true_gene
                ),
                "own_model_mast_e_value": own["mast_e_value"],
                "own_model_significant_hits": own["significant_motif_hits"],
                "own_model_hit_score": own_score,
                "best_other_model": best_other["model_gene"],
                "best_other_hit_score": other_score,
                "own_minus_best_other_score": own_score - other_score,
                "motif_model_status": "own gene model best"
                if ordered[0]["model_gene"] == true_gene
                else "motif-model outlier",
                "manual_priority_note": HIGHLIGHTS.get(accession, ""),
            }
        )
    result.sort(key=lambda row: (GENES.index(str(row["gene"])), str(row["species"])))
    return result


def global_motif_coverage(mast_path: Path) -> list[dict[str, object]]:
    root = ET.parse(mast_path).getroot()
    motif_names = {
        str(i): motif.attrib.get("alt", motif.attrib["id"])
        for i, motif in enumerate(root.findall(".//motifs/motif"))
    }
    coverage: dict[tuple[str, str], set[str]] = defaultdict(set)
    hit_counts: dict[tuple[str, str], int] = defaultdict(int)
    for sequence in root.findall(".//sequences/sequence"):
        gene = sequence.attrib["name"].split("|", 1)[0]
        for hit in sequence.findall(".//hit"):
            motif = motif_names.get(hit.attrib["idx"], f"index_{hit.attrib['idx']}")
            coverage[(motif, gene)].add(sequence.attrib["name"])
            hit_counts[(motif, gene)] += 1
    rows: list[dict[str, object]] = []
    for motif in sorted(set(motif_names.values())):
        for gene in GENES:
            rows.append(
                {
                    "motif": motif,
                    "gene": gene,
                    "sequences_with_hit": len(coverage[(motif, gene)]),
                    "hit_count": hit_counts[(motif, gene)],
                }
            )
    return rows


def plot_model_matrix(scores: list[dict[str, object]], output: Path) -> None:
    matrix = np.zeros((len(GENES), len(GENES)))
    for i, true_gene in enumerate(GENES):
        for j, model_gene in enumerate(GENES):
            values = [
                float(row["motif_hit_score_sum_neg_log10_p"])
                for row in scores
                if row["true_gene"] == true_gene and row["model_gene"] == model_gene
            ]
            matrix[i, j] = float(np.median(values))
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    image = ax.imshow(matrix, cmap="viridis", aspect="auto")
    ax.set_xticks(
        range(len(GENES)),
        labels=[f"{gene} model" for gene in GENES],
        rotation=30,
        ha="right",
    )
    ax.set_yticks(range(len(GENES)), labels=[f"known {gene}" for gene in GENES])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j,
                i,
                f"{matrix[i, j]:.0f}",
                ha="center",
                va="center",
                fontsize=9,
                color="white" if matrix[i, j] > matrix.max() * 0.58 else "black",
            )
    ax.set_title("Gene-specific protein motif-model discrimination")
    fig.colorbar(image, ax=ax, label="Median sum of −log10 motif-hit p-values")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_candidate_margins(assignments: list[dict[str, object]], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 5.5))
    positions: list[float] = []
    labels: list[str] = []
    offset = 0
    colours = dict(
        zip(
            GENES,
            (
                "#4c78a8",
                "#f58518",
                "#54a24b",
                "#e45756",
                "#72b7b2",
                "#b279a2",
                "#8c564b",
            ),
            strict=False,
        )
    )
    for gene in GENES:
        rows = [row for row in assignments if row["gene"] == gene]
        xs = np.arange(offset, offset + len(rows))
        ys = [float(row["own_minus_best_other_score"]) for row in rows]
        ax.scatter(xs, ys, s=25, color=colours[gene], label=gene)
        for x, row, y in zip(xs, rows, ys, strict=False):
            if "tentative" in str(row["manual_priority_note"]):
                ax.annotate(
                    str(row["manual_priority_note"]).split(" (")[0],
                    (x, y),
                    xytext=(3, 5),
                    textcoords="offset points",
                    fontsize=7,
                    rotation=35,
                )
        positions.append(offset + (len(rows) - 1) / 2)
        labels.append(gene)
        offset += len(rows) + 2
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(positions, labels=labels)
    ax.set_ylabel("Own-gene motif score minus best other-gene score")
    ax.set_title("Candidate support from gene-specific protein motif models")
    ax.legend(ncol=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--raw-dir", type=Path)
    args = parser.parse_args()
    base = args.base
    raw = args.raw_dir or base / "raw_outputs"
    tables = base / "tables"
    figures = base / "figures"

    meme_rows = meme_motifs(raw / "meme_all_slrp" / "meme.xml", "all_seven_genes")
    for gene in GENES:
        meme_rows.extend(meme_motifs(raw / f"meme_{gene}" / "meme.xml", gene))
    write_tsv(tables / "protein_meme_motif_summary.tsv", meme_rows, list(meme_rows[0]))

    streme_rows: list[dict[str, object]] = []
    for gene in GENES:
        streme_rows.extend(
            streme_motifs(raw / f"streme_{gene}_vs_other_genes" / "streme.xml", gene)
        )
    write_tsv(
        tables / "protein_streme_gene_enriched_motifs.tsv",
        streme_rows,
        list(streme_rows[0]),
    )

    scores: list[dict[str, object]] = []
    for gene in GENES:
        scores.extend(
            mast_scores(raw / f"mast_{gene}_motifs_all_proteins" / "mast.xml", gene)
        )
    write_tsv(
        tables / "protein_candidate_motif_model_scores.tsv", scores, list(scores[0])
    )
    assignments = candidate_assignments(scores)
    write_tsv(
        tables / "protein_candidate_motif_model_assignments.tsv",
        assignments,
        list(assignments[0]),
    )

    coverage = global_motif_coverage(raw / "mast_all_slrp" / "mast.xml")
    write_tsv(
        tables / "protein_global_motif_gene_coverage.tsv", coverage, list(coverage[0])
    )

    plot_model_matrix(scores, figures / "protein_gene_motif_model_matrix.png")
    plot_candidate_margins(
        assignments, figures / "protein_candidate_motif_model_margins.png"
    )
    outliers = [
        row for row in assignments if row["motif_model_status"] != "own gene model best"
    ]
    print(
        f"Summarized {len(meme_rows)} MEME motifs, {len(streme_rows)} STREME motifs, and {len(scores)} candidate/model scores"
    )
    print(f"Candidates whose own-gene model was not best: {len(outliers)}")
    for row in outliers:
        print(
            f"  {row['gene']} {row['species']} {row['accession']} -> {row['best_motif_model']}"
        )


if __name__ == "__main__":
    main()
