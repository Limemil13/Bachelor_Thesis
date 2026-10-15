#!/usr/bin/env python3
"""Summarize same-gene promoter alignments across representative vertebrates."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from itertools import combinations
from pathlib import Path


GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
WINDOWS = ("core", "proximal")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line.upper())
    return {key: "".join(value) for key, value in records.items()}


def species_name(header: str) -> str:
    parts = header.split("|")
    if len(parts) < 2:
        raise ValueError(f"Unexpected FASTA header: {header}")
    return parts[1]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pairwise_metrics(first: str, second: str) -> tuple[int, int, float, float]:
    compared = 0
    matches = 0
    for a, b in zip(first, second):
        if a == "-" or b == "-":
            continue
        compared += 1
        matches += a == b
    shorter = min(len(first.replace("-", "")), len(second.replace("-", "")))
    identity = 100.0 * matches / compared if compared else 0.0
    coverage = 100.0 * compared / shorter if shorter else 0.0
    return compared, matches, identity, coverage


def alignment_metrics(records: dict[str, str]) -> dict[str, object]:
    sequences = list(records.values())
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Aligned FASTA contains unequal sequence lengths")
    sequence_count = len(sequences)
    minimum_present = math.ceil(0.8 * sequence_count)
    fully_covered = 0
    fully_identical = 0
    recurrent = 0
    longest_recurrent_run = 0
    current_run = 0
    for column in zip(*sequences):
        bases = [base for base in column if base != "-"]
        if len(bases) == sequence_count:
            fully_covered += 1
            if len(set(bases)) == 1:
                fully_identical += 1
        if len(bases) >= minimum_present:
            most_common = max(bases.count(base) for base in set(bases))
            is_recurrent = most_common / len(bases) >= 0.8
        else:
            is_recurrent = False
        if is_recurrent:
            recurrent += 1
            current_run += 1
            longest_recurrent_run = max(longest_recurrent_run, current_run)
        else:
            current_run = 0
    shortest = min(len(sequence.replace("-", "")) for sequence in sequences)
    return {
        "sequence_count": sequence_count,
        "species": ";".join(sorted(species_name(header) for header in records)),
        "alignment_length_bp": lengths.pop(),
        "shortest_promoter_bp": shortest,
        "fully_covered_columns": fully_covered,
        "fully_identical_columns": fully_identical,
        "recurrent_columns": recurrent,
        "longest_recurrent_run_bp": longest_recurrent_run,
        "fully_identical_percent_of_shortest": round(100 * fully_identical / shortest, 2),
        "recurrent_percent_of_shortest": round(100 * recurrent / shortest, 2),
    }


def heatmap_svg(rows: list[dict[str, object]], output: Path) -> None:
    values = {(str(row["gene"]), str(row["window"])): float(row["median_pairwise_identity_percent"]) for row in rows}
    cell_w = 150
    cell_h = 54
    left = 90
    top = 72
    width = left + cell_w * len(WINDOWS) + 30
    height = top + cell_h * len(GENES) + 58
    all_values = list(values.values())
    minimum = min(all_values)
    maximum = max(all_values)
    span = maximum - minimum or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.label{font-size:14px}.title{font-size:18px;font-weight:bold}.value{font-size:15px;font-weight:bold}</style>',
        f'<text class="title" x="{width / 2}" y="28" text-anchor="middle">Cross-species promoter sequence identity</text>',
    ]
    for column, window in enumerate(WINDOWS):
        x = left + column * cell_w + cell_w / 2
        parts.append(f'<text class="label" x="{x}" y="{top - 14}" text-anchor="middle">{window}</text>')
    for row_index, gene in enumerate(GENES):
        y = top + row_index * cell_h + cell_h / 2 + 5
        parts.append(f'<text class="label" x="{left - 12}" y="{y}" text-anchor="end">{gene}</text>')
        for column, window in enumerate(WINDOWS):
            value = values[(gene, window)]
            fraction = (value - minimum) / span
            red = int(239 - 185 * fraction)
            green = int(247 - 95 * fraction)
            blue = int(255 - 40 * fraction)
            x = left + column * cell_w
            y0 = top + row_index * cell_h
            parts.append(f'<rect x="{x}" y="{y0}" width="{cell_w}" height="{cell_h}" fill="rgb({red},{green},{blue})" stroke="#d9d9d9"/>')
            text_fill = "white" if fraction > 0.65 else "#111"
            parts.append(f'<text class="value" x="{x + cell_w / 2}" y="{y0 + cell_h / 2 + 5}" text-anchor="middle" style="fill:{text_fill}">{value:.1f}%</text>')
    parts.append(f'<text class="label" x="{width / 2}" y="{height - 12}" text-anchor="middle">Median identity across all species pairs</text>')
    parts.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    args = parser.parse_args()
    results = args.analysis_dir / "results"
    pair_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for window in WINDOWS:
        for gene in GENES:
            records = read_fasta(results / "alignments" / window / f"{gene}.aln.fna")
            pair_identities: list[float] = []
            pair_coverages: list[float] = []
            for (header_a, sequence_a), (header_b, sequence_b) in combinations(records.items(), 2):
                compared, matches, identity, coverage = pairwise_metrics(sequence_a, sequence_b)
                pair_identities.append(identity)
                pair_coverages.append(coverage)
                pair_rows.append(
                    {
                        "gene": gene,
                        "window": window,
                        "species_a": species_name(header_a),
                        "species_b": species_name(header_b),
                        "compared_bp": compared,
                        "matching_bp": matches,
                        "identity_percent": round(identity, 2),
                        "shorter_sequence_coverage_percent": round(coverage, 2),
                    }
                )
            summary = alignment_metrics(records)
            summary_rows.append(
                {
                    "gene": gene,
                    "window": window,
                    **summary,
                    "pairwise_comparisons": len(pair_identities),
                    "median_pairwise_identity_percent": round(statistics.median(pair_identities), 2),
                    "minimum_pairwise_identity_percent": round(min(pair_identities), 2),
                    "maximum_pairwise_identity_percent": round(max(pair_identities), 2),
                    "median_pairwise_coverage_percent": round(statistics.median(pair_coverages), 2),
                }
            )
    tables = results / "tables"
    write_tsv(tables / "promoter_pairwise_similarity.tsv", pair_rows, list(pair_rows[0]))
    write_tsv(tables / "promoter_alignment_summary.tsv", summary_rows, list(summary_rows[0]))
    heatmap_svg(summary_rows, results / "figures" / "promoter_cross_species_identity.svg")
    print(f"Summarized {len(summary_rows)} alignments and {len(pair_rows)} pairwise comparisons")


if __name__ == "__main__":
    main()
