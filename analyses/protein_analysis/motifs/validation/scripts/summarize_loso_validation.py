#!/usr/bin/env python3
"""Summarize independent leave-one-species-out protein motif validation."""

from __future__ import annotations

import argparse
import csv
import math
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def neg_log10(value: str, cap: float = 300.0) -> float:
    number = float(value)
    return cap if number <= 0 else min(-math.log10(number), cap)


def parse_mast(path: Path, model_gene: str, fold: str) -> list[dict[str, object]]:
    root = ET.parse(path).getroot()
    rows: list[dict[str, object]] = []
    for sequence in root.findall(".//sequences/sequence"):
        score = sequence.find("score")
        if score is None:
            continue
        true_gene, species, accession = sequence.attrib["name"].split("|", 2)
        hits = sequence.findall(".//hit")
        hit_score = sum(neg_log10(hit.attrib["pvalue"], cap=100.0) for hit in hits)
        rows.append(
            {
                "fold": fold,
                "true_gene": true_gene,
                "species": species,
                "accession": accession,
                "model_gene": model_gene,
                "mast_e_value": score.attrib["evalue"],
                "mast_neg_log10_e_value": neg_log10(score.attrib["evalue"]),
                "significant_motif_hits": len(hits),
                "motif_hit_score": hit_score,
            }
        )
    return rows


def assignments(
    scores: list[dict[str, object]], manifest: list[dict[str, str]]
) -> list[dict[str, object]]:
    exact = {
        (row["fold"], row["gene"], row["species"], row["accession"]): row[
            "exact_sequence_present_in_training"
        ]
        for row in manifest
    }
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in scores:
        grouped[
            (
                str(row["fold"]),
                str(row["true_gene"]),
                str(row["species"]),
                str(row["accession"]),
            )
        ].append(row)

    output: list[dict[str, object]] = []
    for key, rows in grouped.items():
        fold, true_gene, species, accession = key
        if len(rows) != len(GENES):
            raise ValueError(f"Expected seven model scores for {key}, found {len(rows)}")
        ranked_e = sorted(
            rows,
            key=lambda row: float(row["mast_neg_log10_e_value"]),
            reverse=True,
        )
        ranked_hits = sorted(
            rows,
            key=lambda row: float(row["motif_hit_score"]),
            reverse=True,
        )
        best, second = ranked_e[:2]
        output.append(
            {
                "fold": fold,
                "true_gene": true_gene,
                "species": species,
                "accession": accession,
                "predicted_gene": best["model_gene"],
                "correct": "yes" if best["model_gene"] == true_gene else "no",
                "best_mast_score": best["mast_neg_log10_e_value"],
                "second_gene": second["model_gene"],
                "mast_score_margin": float(best["mast_neg_log10_e_value"])
                - float(second["mast_neg_log10_e_value"]),
                "hit_score_prediction": ranked_hits[0]["model_gene"],
                "two_scoring_rules_agree": (
                    "yes" if ranked_hits[0]["model_gene"] == best["model_gene"] else "no"
                ),
                "exact_sequence_present_in_training": exact[key],
            }
        )
    output.sort(key=lambda row: (str(row["fold"]), str(row["true_gene"])))
    return output


def aggregate(
    predictions: list[dict[str, object]], key: str
) -> list[dict[str, object]]:
    values = GENES if key == "true_gene" else sorted({str(row[key]) for row in predictions})
    rows: list[dict[str, object]] = []
    for value in values:
        subset = [row for row in predictions if row[key] == value]
        correct = sum(row["correct"] == "yes" for row in subset)
        nonduplicate = [
            row for row in subset if row["exact_sequence_present_in_training"] == "no"
        ]
        nonduplicate_correct = sum(row["correct"] == "yes" for row in nonduplicate)
        rows.append(
            {
                key: value,
                "tested": len(subset),
                "correct": correct,
                "accuracy_percent": round(100 * correct / len(subset), 2) if subset else "",
                "tested_without_exact_training_duplicate": len(nonduplicate),
                "correct_without_exact_training_duplicate": nonduplicate_correct,
                "accuracy_without_exact_training_duplicate_percent": (
                    round(100 * nonduplicate_correct / len(nonduplicate), 2)
                    if nonduplicate
                    else ""
                ),
            }
        )
    return rows


def summarize_challenges(
    scores: list[dict[str, object]], excluded: list[dict[str, str]]
) -> list[dict[str, object]]:
    labels = {row["accession"]: row for row in excluded}
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in scores:
        grouped[str(row["accession"])].append(row)
    output: list[dict[str, object]] = []
    for accession, rows in grouped.items():
        ranked_e = sorted(
            rows,
            key=lambda row: float(row["mast_neg_log10_e_value"]),
            reverse=True,
        )
        ranked_hits = sorted(
            rows,
            key=lambda row: float(row["motif_hit_score"]),
            reverse=True,
        )
        best, second = ranked_e[:2]
        label = labels[accession]
        output.append(
            {
                "species": best["species"],
                "accession": accession,
                "original_target_label": label["original_gene"],
                "revised_reference_tree_label": label["revised_label"],
                "nearest_forced_target_model": best["model_gene"],
                "second_target_model": second["model_gene"],
                "mast_score_margin": float(best["mast_neg_log10_e_value"])
                - float(second["mast_neg_log10_e_value"]),
                "hit_score_nearest_model": ranked_hits[0]["model_gene"],
                "scoring_rules_agree": (
                    "yes" if ranked_hits[0]["model_gene"] == best["model_gene"] else "no"
                ),
                "interpretation": "nearest target model only; control is outside the seven target labels",
            }
        )
    return sorted(output, key=lambda row: str(row["accession"]))


def plot_confusion(predictions: list[dict[str, object]], output: Path) -> None:
    matrix = [[0 for _ in GENES] for _ in GENES]
    for row in predictions:
        matrix[GENES.index(str(row["true_gene"]))][
            GENES.index(str(row["predicted_gene"]))
        ] += 1

    # Write a compact SVG directly so that the summary step does not depend on
    # matplotlib, which is not installed in the lightweight motif environment.
    cell = 64
    left = 105
    top = 82
    width = left + cell * len(GENES) + 24
    height = top + cell * len(GENES) + 76
    maximum = max(max(row) for row in matrix) or 1
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.label{font-size:14px}.title{font-size:18px;font-weight:bold}.count{font-size:16px;font-weight:bold}</style>',
        f'<text class="title" x="{width / 2}" y="28" text-anchor="middle">Leave-one-species-out protein-pattern validation</text>',
        f'<text class="label" x="{left + cell * len(GENES) / 2}" y="{height - 12}" text-anchor="middle">Predicted gene</text>',
        f'<text class="label" x="18" y="{top + cell * len(GENES) / 2}" text-anchor="middle" transform="rotate(-90 18 {top + cell * len(GENES) / 2})">Curated gene</text>',
    ]
    for index, gene in enumerate(GENES):
        x = left + index * cell + cell / 2
        y = top - 12
        parts.append(f'<text class="label" x="{x}" y="{y}" text-anchor="middle">{gene}</text>')
        row_y = top + index * cell + cell / 2 + 5
        parts.append(f'<text class="label" x="{left - 12}" y="{row_y}" text-anchor="end">{gene}</text>')
    for row_index, row in enumerate(matrix):
        for column_index, value in enumerate(row):
            intensity = value / maximum
            blue = int(247 - 150 * intensity)
            red = int(247 - 210 * intensity)
            green = int(251 - 125 * intensity)
            fill = f"rgb({red},{green},{blue})"
            x = left + column_index * cell
            y = top + row_index * cell
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="#d9d9d9"/>'
            )
            text_fill = "white" if intensity > 0.55 else "#111"
            parts.append(
                f'<text class="count" x="{x + cell / 2}" y="{y + cell / 2 + 6}" text-anchor="middle" style="fill:{text_fill}">{value}</text>'
            )
    parts.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    args = parser.parse_args()
    inputs = args.analysis_dir / "inputs"
    results = args.analysis_dir / "results"
    manifest = read_tsv(inputs / "validation_manifest.tsv")

    scores: list[dict[str, object]] = []
    for fold_dir in sorted((results / "folds").iterdir()):
        if not fold_dir.is_dir():
            continue
        for gene in GENES:
            scores.extend(
                parse_mast(
                    fold_dir / gene / "mast" / "mast.xml",
                    model_gene=gene,
                    fold=fold_dir.name,
                )
            )
    predictions = assignments(scores, manifest)
    if len(predictions) != 101:
        raise SystemExit(f"Expected predictions for 101 target proteins, found {len(predictions)}")

    tables = results / "tables"
    write_tsv(tables / "loso_model_scores.tsv", scores, list(scores[0]))
    write_tsv(tables / "loso_predictions.tsv", predictions, list(predictions[0]))
    gene_rows = aggregate(predictions, "true_gene")
    species_rows = aggregate(predictions, "species")
    write_tsv(tables / "loso_gene_summary.tsv", gene_rows, list(gene_rows[0]))
    write_tsv(tables / "loso_species_summary.tsv", species_rows, list(species_rows[0]))

    confusion_rows = []
    for true_gene in GENES:
        counts = Counter(
            str(row["predicted_gene"])
            for row in predictions
            if row["true_gene"] == true_gene
        )
        confusion_rows.append(
            {"true_gene": true_gene, **{gene: counts[gene] for gene in GENES}}
        )
    write_tsv(
        tables / "loso_confusion_matrix.tsv",
        confusion_rows,
        ["true_gene", *GENES],
    )

    correct = sum(row["correct"] == "yes" for row in predictions)
    nonduplicate = [
        row for row in predictions if row["exact_sequence_present_in_training"] == "no"
    ]
    nonduplicate_correct = sum(row["correct"] == "yes" for row in nonduplicate)
    scoring_agreement = sum(row["two_scoring_rules_agree"] == "yes" for row in predictions)
    summary_rows = [
        {"metric": "tested_proteins", "value": len(predictions)},
        {"metric": "correct_predictions", "value": correct},
        {"metric": "overall_accuracy_percent", "value": round(100 * correct / len(predictions), 2)},
        {"metric": "proteins_without_exact_training_duplicate", "value": len(nonduplicate)},
        {"metric": "correct_without_exact_training_duplicate", "value": nonduplicate_correct},
        {
            "metric": "accuracy_without_exact_training_duplicate_percent",
            "value": round(100 * nonduplicate_correct / len(nonduplicate), 2),
        },
        {"metric": "mast_and_hit_score_predictions_agree", "value": scoring_agreement},
    ]
    write_tsv(tables / "loso_summary.tsv", summary_rows, ["metric", "value"])

    challenge_scores: list[dict[str, object]] = []
    for gene in GENES:
        challenge_scores.extend(
            parse_mast(
                results / "final_models" / gene / "mast" / "mast.xml",
                model_gene=gene,
                fold="challenge_controls",
            )
        )
    excluded = read_tsv(inputs / "excluded_from_target_models.tsv")
    challenge_rows = summarize_challenges(challenge_scores, excluded)
    write_tsv(
        tables / "challenge_control_model_scores.tsv",
        challenge_scores,
        list(challenge_scores[0]),
    )
    write_tsv(
        tables / "challenge_control_nearest_models.tsv",
        challenge_rows,
        list(challenge_rows[0]),
    )
    plot_confusion(predictions, results / "figures" / "loso_confusion_matrix.svg")
    print(
        f"LOSO validation: {correct}/{len(predictions)} correct; "
        f"excluding exact training duplicates: {nonduplicate_correct}/{len(nonduplicate)} correct"
    )


if __name__ == "__main__":
    main()
