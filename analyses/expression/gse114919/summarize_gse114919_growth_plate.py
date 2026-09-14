#!/usr/bin/env python3
"""Summarize the published GSE114919 mouse/rat growth-plate matrices.

Values are kept on the authors' published normalized-value scale.  The script
does not combine those values with local Salmon TPMs or compare absolute mouse
and rat magnitudes.  Cross-species interpretation uses detection and ranks
calculated separately within each species/age/bone/zone condition.
"""

from __future__ import annotations

import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
GENES = ("ASPN", "BGN", "DCN", "EPYC", "FMOD", "LUM", "OGN", "OMD", "PRELP")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def rat_name_map() -> dict[int, str]:
    rows = read_tsv(TABLES / "gse114919_rat_sample_name_map.tsv")
    return {int(row["No."]): row["Sample "].strip() for row in rows}


def parse_mouse(sample: str) -> dict[str, str | int]:
    match = re.fullmatch(r"([14])w(T|Ph|P)_(HZ|PZ)(\d+)", sample)
    if not match:
        raise ValueError(f"Unrecognized mouse sample header: {sample}")
    age, bone_code, zone, replicate = match.groups()
    return {
        "published_sample_name": sample,
        "age_weeks": int(age),
        "bone": "tibia" if bone_code == "T" else "phalanx",
        "zone": zone,
        "replicate": int(replicate),
    }


def parse_rat(sample: str, names: dict[int, str]) -> dict[str, str | int]:
    number_match = re.search(r"_S(\d+)$", sample)
    if not number_match:
        raise ValueError(f"Unrecognized rat matrix header: {sample}")
    number = int(number_match.group(1))
    published_name = names[number]
    match = re.fullmatch(r"(T|Ph)([14])wk (PZ|HZ)(\d+)", published_name)
    if not match:
        raise ValueError(f"Unrecognized rat sample-map name: {published_name}")
    bone_code, age, zone, replicate = match.groups()
    return {
        "published_sample_name": published_name,
        "age_weeks": int(age),
        "bone": "tibia" if bone_code == "T" else "phalanx",
        "zone": zone,
        "replicate": int(replicate),
    }


def mean_sd_cv(values: list[float]) -> tuple[float, float, float | str]:
    mean = statistics.fmean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    cv = sd / mean if mean else ""
    return mean, sd, cv


def main() -> None:
    rat_names = rat_name_map()
    long_rows: list[dict] = []
    sources = {
        "mouse": TABLES / "gse114919_mouse_slrp_normalized_counts.tsv",
        "rat": TABLES / "gse114919_rat_slrp_normalized_counts.tsv",
    }
    for species, source in sources.items():
        source_rows = read_tsv(source)
        expected_samples = 29 if species == "mouse" else 30
        assert len(source_rows) == len(GENES) * expected_samples
        assert {x["gene"] for x in source_rows} == set(GENES)
        for row in source_rows:
            metadata = (
                parse_mouse(row["sample"])
                if species == "mouse"
                else parse_rat(row["sample"], rat_names)
            )
            long_rows.append(
                {
                    **row,
                    **metadata,
                    "normalized_count": float(row["normalized_count"]),
                    "detected": "yes" if float(row["normalized_count"]) > 0 else "no",
                    "source_accession": "GSE114919",
                    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114919",
                }
            )

    long_fields = [
        "species",
        "gene",
        "matrix_gene_id",
        "ensembl_gene_id",
        "ncbi_gene_id",
        "source_row",
        "sample",
        "published_sample_name",
        "bone",
        "age_weeks",
        "zone",
        "replicate",
        "normalized_count",
        "detected",
        "source_accession",
        "source_url",
    ]
    write_tsv(TABLES / "gse114919_slrp_growth_plate_long.tsv", long_rows, long_fields)

    groups: dict[tuple, list[float]] = defaultdict(list)
    for row in long_rows:
        key = (row["species"], row["bone"], row["age_weeks"], row["zone"], row["gene"])
        groups[key].append(float(row["normalized_count"]))

    summary_rows: list[dict] = []
    for (species, bone, age, zone, gene), values in sorted(groups.items()):
        mean, sd, cv = mean_sd_cv(values)
        summary_rows.append(
            {
                "species": species,
                "bone": bone,
                "age_weeks": age,
                "zone": zone,
                "gene": gene,
                "replicate_count": len(values),
                "detected_replicates": sum(x > 0 for x in values),
                "detection_fraction": sum(x > 0 for x in values) / len(values),
                "mean_published_normalized_value": round(mean, 6),
                "sd_published_normalized_value": round(sd, 6),
                "cv": round(cv, 6) if cv != "" else "",
                "minimum": round(min(values), 6),
                "maximum": round(max(values), 6),
            }
        )

    by_condition: dict[tuple, list[dict]] = defaultdict(list)
    for row in summary_rows:
        by_condition[
            (row["species"], row["bone"], row["age_weeks"], row["zone"])
        ].append(row)
    for condition_rows in by_condition.values():
        ordered = sorted(
            condition_rows,
            key=lambda x: (-float(x["mean_published_normalized_value"]), x["gene"]),
        )
        for rank, row in enumerate(ordered, start=1):
            row["within_condition_rank_of_9"] = rank

    summary_fields = [
        "species",
        "bone",
        "age_weeks",
        "zone",
        "gene",
        "replicate_count",
        "detected_replicates",
        "detection_fraction",
        "mean_published_normalized_value",
        "sd_published_normalized_value",
        "cv",
        "minimum",
        "maximum",
        "within_condition_rank_of_9",
    ]
    write_tsv(
        TABLES / "gse114919_slrp_condition_summary.tsv", summary_rows, summary_fields
    )

    tibia_rows = [x for x in summary_rows if x["bone"] == "tibia"]
    by_species_gene: dict[tuple, list[dict]] = defaultdict(list)
    for row in tibia_rows:
        by_species_gene[(row["species"], row["gene"])].append(row)
    gene_rows: list[dict] = []
    for (species, gene), rows in sorted(by_species_gene.items()):
        gene_rows.append(
            {
                "species": species,
                "gene": gene,
                "tibia_conditions": len(rows),
                "conditions_detected_in_all_replicates": sum(
                    x["detected_replicates"] == x["replicate_count"] for x in rows
                ),
                "mean_within_condition_rank": round(
                    statistics.fmean(
                        float(x["within_condition_rank_of_9"]) for x in rows
                    ),
                    3,
                ),
                "best_within_condition_rank": min(
                    int(x["within_condition_rank_of_9"]) for x in rows
                ),
                "worst_within_condition_rank": max(
                    int(x["within_condition_rank_of_9"]) for x in rows
                ),
                "mean_of_condition_means_descriptive_only": round(
                    statistics.fmean(
                        float(x["mean_published_normalized_value"]) for x in rows
                    ),
                    6,
                ),
                "interpretation_scale": "within-species ranks; do not compare absolute values between species",
            }
        )
    gene_fields = list(gene_rows[0])
    write_tsv(
        TABLES / "gse114919_slrp_tibia_cross_condition_summary.tsv",
        gene_rows,
        gene_fields,
    )

    # Auditable within-species, within-zone age contrasts for tibial growth plate.
    lookup = {
        (x["species"], x["gene"], x["age_weeks"], x["zone"]): x for x in tibia_rows
    }
    contrast_rows: list[dict] = []
    for species in ("mouse", "rat"):
        for gene in GENES:
            for zone in ("PZ", "HZ"):
                younger = lookup[(species, gene, 1, zone)]
                older = lookup[(species, gene, 4, zone)]
                contrast_rows.append(
                    {
                        "species": species,
                        "gene": gene,
                        "zone": zone,
                        "one_week_mean": younger["mean_published_normalized_value"],
                        "four_week_mean": older["mean_published_normalized_value"],
                        "four_minus_one_week_delta": round(
                            float(older["mean_published_normalized_value"])
                            - float(younger["mean_published_normalized_value"]),
                            6,
                        ),
                        "scale_note": "difference on authors' published normalized-value scale; not a fold change",
                    }
                )
    write_tsv(
        TABLES / "gse114919_slrp_tibia_age_contrasts.tsv",
        contrast_rows,
        list(contrast_rows[0]),
    )

    print(f"Long rows: {len(long_rows)} (mouse 9x29; rat 9x30)")
    print(f"Condition summaries: {len(summary_rows)}")
    print(f"Tibia species-gene summaries: {len(gene_rows)}")
    print(
        "Mouse source workbook lacks the published 1w phalanx PZ replicate 3 column; tibia groups remain n=5."
    )


if __name__ == "__main__":
    main()
