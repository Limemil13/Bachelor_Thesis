#!/usr/bin/env python3
"""Summarize promoter extraction, STREME/Tomtom, AME, and FIMO outputs."""

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
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
FAMILY_ORDER = (
    "SOX5/6/9",
    "RUNX",
    "SMAD",
    "HIF",
    "GLI",
    "AP-1",
    "NF-kB",
    "WNT/TCF",
    "MEF2C",
    "STAT3",
    "CREB1",
    "ATF4",
    "FOXA2",
    "NFATC1",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def family_for(name: str) -> str:
    upper = name.upper()
    if "SOX5" in upper or "SOX6" in upper or "SOX9" in upper:
        return "SOX5/6/9"
    if "RUNX" in upper:
        return "RUNX"
    if "SMAD" in upper:
        return "SMAD"
    if "HIF1A" in upper:
        return "HIF"
    if "GLI" in upper:
        return "GLI"
    if "FOS" in upper or "JUN" in upper or "BATF" in upper:
        return "AP-1"
    if "NFKB" in upper or "RELA" in upper:
        return "NF-kB"
    if "TCF7L2" in upper or "LEF1" in upper:
        return "WNT/TCF"
    for family in ("MEF2C", "STAT3", "CREB1", "ATF4", "FOXA2", "NFATC1"):
        if family in upper:
            return family
    return "other"


def promoter_qc(metadata: list[dict[str, str]]) -> list[dict[str, object]]:
    return [
        {
            "gene": row["gene"],
            "species": row["species"],
            "window": row["window"],
            "transcript_id": row["transcript_id"],
            "strand": row["strand"],
            "length_bp": row["sequence_length_bp"],
            "expected_length_bp": row["expected_length_bp"],
            "clipped_at_contig_edge": row["clipped_at_contig_edge"],
            "gc_percent": row["gc_percent"],
            "ambiguous_percent": row["ambiguous_percent"],
        }
        for row in metadata
    ]


def parse_streme(
    path: Path, scope: str, tomtom_path: Path, jaspar_names: dict[str, str]
) -> list[dict[str, object]]:
    root = ET.parse(path).getroot()
    matches: dict[str, list[dict[str, str]]] = defaultdict(list)
    if tomtom_path.exists():
        for row in read_tsv(tomtom_path):
            matches[row["Query_ID"]].append(row)
    rows: list[dict[str, object]] = []
    for motif in root.findall(".//motif"):
        motif_id = motif.attrib["id"]
        candidates = sorted(
            matches.get(motif_id, []), key=lambda x: float(x["q-value"])
        )
        best = candidates[0] if candidates else None
        rows.append(
            {
                "scope": scope,
                "streme_motif_id": motif_id,
                "consensus": motif_id.split("-", 1)[1] if "-" in motif_id else motif_id,
                "width": motif.attrib.get("width", ""),
                "train_positive_sequences": motif.attrib.get("train_pos_count", ""),
                "train_negative_sequences": motif.attrib.get("train_neg_count", ""),
                "training_p_value": motif.attrib.get("train_pvalue", ""),
                "holdout_positive_sequences": motif.attrib.get("test_pos_count", ""),
                "holdout_negative_sequences": motif.attrib.get("test_neg_count", ""),
                "holdout_p_value": motif.attrib.get("test_pvalue", ""),
                "best_jaspar_motif_id": best["Target_ID"] if best else "",
                "best_jaspar_name": jaspar_names.get(best["Target_ID"], "")
                if best
                else "",
                "tomtom_q_value": best["q-value"] if best else "",
                "validation_status": "training-only; no holdout available"
                if motif.attrib.get("test_pos_count") == "0"
                else "holdout tested",
            }
        )
    return rows


def parse_jaspar_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MOTIF "):
                fields = line.split()
                names[fields[1]] = fields[2] if len(fields) > 2 else ""
    return names


def fimo_summary(
    fimo_rows: list[dict[str, str]], jaspar_names: dict[str, str]
) -> tuple[list[dict[str, object]], np.ndarray]:
    hits: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in fimo_rows:
        gene, species = row["sequence_name"].split("|", 2)[:2]
        row = dict(row)
        row["gene"] = gene
        row["species"] = species
        hits[(gene, row["motif_id"])].append(row)

    result: list[dict[str, object]] = []
    family_species: dict[tuple[str, str], set[str]] = defaultdict(set)
    for gene in GENES:
        for motif_id, name in jaspar_names.items():
            motif_hits = hits.get((gene, motif_id), [])
            species = sorted(
                {row["species"] for row in motif_hits}, key=lambda x: SPECIES.index(x)
            )
            family = family_for(name)
            family_species[(gene, family)].update(species)
            result.append(
                {
                    "gene": gene,
                    "motif_id": motif_id,
                    "motif_name": name,
                    "motif_family": family,
                    "species_with_p_le_1e-4": len(species),
                    "species": ";".join(species),
                    "hit_count": len(motif_hits),
                    "minimum_p_value": min(
                        (float(row["p-value"]) for row in motif_hits), default=""
                    ),
                    "minimum_reported_q_value": min(
                        (float(row["q-value"]) for row in motif_hits), default=""
                    ),
                    "interpretation": "exploratory uncorrected recurrence; not binding evidence",
                }
            )

    matrix = np.zeros((len(GENES), len(FAMILY_ORDER)), dtype=int)
    for i, gene in enumerate(GENES):
        for j, family in enumerate(FAMILY_ORDER):
            matrix[i, j] = len(family_species.get((gene, family), set()))
    return result, matrix


def plot_gc(metadata: list[dict[str, str]], output: Path) -> None:
    matrix = np.zeros((len(GENES), len(SPECIES)))
    by_key = {
        (row["gene"], row["species"]): float(row["gc_percent"])
        for row in metadata
        if row["window"] == "proximal"
    }
    for i, gene in enumerate(GENES):
        for j, species in enumerate(SPECIES):
            matrix[i, j] = by_key[(gene, species)]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    image = ax.imshow(matrix, cmap="YlGnBu", aspect="auto", vmin=25, vmax=75)
    ax.set_xticks(range(len(SPECIES)), labels=SPECIES, rotation=25, ha="right")
    ax.set_yticks(range(len(GENES)), labels=GENES)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.1f}", ha="center", va="center", fontsize=9)
    ax.set_title("GC content of 2 kb/+200 bp promoter windows")
    fig.colorbar(image, ax=ax, label="GC (%)")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fimo(matrix: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 5))
    image = ax.imshow(matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=5)
    ax.set_xticks(
        range(len(FAMILY_ORDER)), labels=FAMILY_ORDER, rotation=35, ha="right"
    )
    ax.set_yticks(range(len(GENES)), labels=GENES)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center", fontsize=9)
    ax.set_title("Exploratory recurrence of targeted TF motifs across five species")
    ax.set_xlabel("Motif family (FIMO p≤1e-4; uncorrected, not binding evidence)")
    fig.colorbar(image, ax=ax, label="Species with ≥1 motif hit")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ame(rows: list[dict[str, str]], output: Path) -> None:
    ordered = sorted(rows, key=lambda row: float(row["E-value"]))[:15]
    labels = [row["motif_alt_ID"] for row in ordered][::-1]
    values = [-math.log10(max(float(row["E-value"]), 1e-300)) for row in ordered][::-1]
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.barh(labels, values, color="#4c78a8")
    ax.axvline(-math.log10(0.05), color="#b22222", linestyle="--", label="E=0.05")
    ax.set_xlabel("−log10 AME E-value (motif-panel corrected)")
    ax.set_title("Targeted TF-motif enrichment versus dinucleotide-shuffled promoters")
    ax.legend(frameon=False)
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
    tables = base / "tables"
    figures = base / "figures"
    raw = args.raw_dir or base / "raw_outputs"

    metadata = read_tsv(base / "inputs" / "promoter_metadata.tsv")
    qc = promoter_qc(metadata)
    write_tsv(tables / "promoter_sequence_qc.tsv", qc, list(qc[0]))

    jaspar = parse_jaspar_names(
        base / "databases" / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
    )
    targeted = parse_jaspar_names(
        base / "databases" / "JASPAR2026_cartilage_growth_plate_targeted.meme"
    )

    denovo: list[dict[str, object]] = []
    denovo.extend(
        parse_streme(
            raw / "streme_all_vs_shuffled" / "streme.xml",
            "all_vs_shuffled",
            raw / "tomtom_all_vs_jaspar" / "tomtom.tsv",
            jaspar,
        )
    )
    for gene in GENES:
        denovo.extend(
            parse_streme(
                raw / f"streme_{gene}_vs_other_genes" / "streme.xml",
                gene,
                raw / f"tomtom_{gene}_vs_jaspar" / "tomtom.tsv",
                jaspar,
            )
        )
    write_tsv(tables / "promoter_denovo_motif_summary.tsv", denovo, list(denovo[0]))

    all_ame = read_tsv(raw / "ame_all_jaspar_vs_shuffled" / "ame.tsv")
    significant = [
        dict(row, significance="E-value < 0.05")
        for row in all_ame
        if float(row["E-value"]) < 0.05
    ]
    if significant:
        write_tsv(
            tables / "promoter_ame_all_jaspar_significant.tsv",
            significant,
            list(significant[0]),
        )
    else:
        write_tsv(
            tables / "promoter_ame_all_jaspar_significant.tsv",
            [],
            list(all_ame[0]) + ["significance"],
        )

    targeted_ame = read_tsv(raw / "ame_targeted_vs_shuffled" / "ame.tsv")
    targeted_rows: list[dict[str, object]] = []
    for row in targeted_ame:
        targeted_rows.append(
            dict(row, panel_significant="yes" if float(row["E-value"]) < 0.05 else "no")
        )
    write_tsv(
        tables / "promoter_ame_targeted_tf_summary.tsv",
        targeted_rows,
        list(targeted_rows[0]),
    )

    strict_fimo = read_tsv(raw / "fimo_targeted" / "fimo.tsv")
    exploratory_fimo = read_tsv(raw / "fimo_targeted_p1e4_exploratory" / "fimo.tsv")
    fimo_rows, matrix = fimo_summary(exploratory_fimo, targeted)
    write_tsv(
        tables / "promoter_targeted_tf_exploratory_recurrence.tsv",
        fimo_rows,
        list(fimo_rows[0]),
    )
    strict_summary = [
        {
            "scan": "targeted JASPAR FIMO q<=0.05",
            "motifs_tested": len(targeted),
            "promoters_tested": len(
                {row["header"] for row in metadata if row["window"] == "proximal"}
            ),
            "reported_hits": len(strict_fimo),
            "interpretation": "no motif instances survived global FIMO q-value correction"
            if not strict_fimo
            else "corrected motif instances detected",
        }
    ]
    write_tsv(
        tables / "promoter_targeted_tf_strict_scan_summary.tsv",
        strict_summary,
        list(strict_summary[0]),
    )

    plot_gc(metadata, figures / "promoter_proximal_gc_heatmap.png")
    plot_fimo(matrix, figures / "promoter_targeted_tf_species_recurrence.png")
    plot_ame(targeted_ame, figures / "promoter_targeted_tf_ame_enrichment.png")
    print(
        f"Summarized {len(metadata)} promoter windows, {len(denovo)} de novo motifs, and {len(exploratory_fimo)} exploratory FIMO hits"
    )


if __name__ == "__main__":
    main()
