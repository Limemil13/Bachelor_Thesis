#!/usr/bin/env python3
"""Check cartilage and hypertrophic-zone markers in GSE114919.

The deposited mouse and rat matrices contain laser-capture-microdissected
proliferative-zone (PZ) and hypertrophic-zone (HZ) samples.  This script checks
those labels descriptively with five familiar markers.  SOX9, COL2A1 and ACAN
are used as cartilage-lineage/matrix controls; COL10A1 and MMP13 are expected to
be higher in hypertrophic cartilage.  Values remain on the authors' published
normalized scale, and no absolute mouse-versus-rat comparison is made.
"""

from __future__ import annotations

import csv
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
SOURCE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114919"

MARKERS = {
    "SOX9": {
        "role": "cartilage_lineage_control",
        "expected_zone": "none",
    },
    "COL2A1": {
        "role": "cartilage_matrix_control",
        "expected_zone": "none",
    },
    "ACAN": {
        "role": "cartilage_matrix_control",
        "expected_zone": "none",
    },
    "COL10A1": {
        "role": "hypertrophic_marker",
        "expected_zone": "HZ",
    },
    "MMP13": {
        "role": "hypertrophic_marker",
        "expected_zone": "HZ",
    },
}


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write a table with a stable, explicit column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_rat_name_map() -> dict[int, str]:
    """Link rat matrix suffixes to the sample names deposited with the study."""
    path = TABLES / "gse114919_rat_sample_name_map.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            int(row["No."]): row["Sample "].strip()
            for row in csv.DictReader(handle, delimiter="\t")
        }


def parse_mouse_sample(sample: str) -> dict[str, str | int]:
    """Read age, bone, zone and replicate from a mouse column name."""
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


def parse_rat_sample(sample: str, name_map: dict[int, str]) -> dict[str, str | int]:
    """Decode a rat matrix column through the deposited companion name map."""
    suffix = re.search(r"_S(\d+)$", sample)
    if not suffix:
        raise ValueError(f"Unrecognized rat sample header: {sample}")
    published_name = name_map[int(suffix.group(1))]
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


def checked_values(row: pd.Series, sample_columns: list[str]) -> pd.Series:
    """Convert one matrix row to numbers without silently filling missing cells."""
    values = pd.to_numeric(row[sample_columns], errors="coerce")
    if values.isna().any():
        missing = ", ".join(values.index[values.isna()].tolist())
        raise ValueError(f"Missing or nonnumeric values in columns: {missing}")
    return values.astype(float)


def extract_mouse() -> tuple[list[dict], list[dict]]:
    """Extract exact marker-symbol rows from the mouse workbook."""
    matrix = pd.read_excel(RAW / "GSE114919_Mouse_normalizedcounts.xlsx")
    symbol_column = matrix.columns[0]
    sample_columns = list(matrix.columns[1:])
    symbols = matrix[symbol_column].astype(str).str.upper()
    long_rows: list[dict] = []
    availability: list[dict] = []

    for gene, metadata in MARKERS.items():
        selected = matrix[symbols == gene]
        if len(selected) != 1:
            raise ValueError(f"Expected one mouse row for {gene}; found {len(selected)}")
        source_row = int(selected.index[0]) + 2
        availability.append(
            {
                "species": "mouse",
                "gene": gene,
                **metadata,
                "available_in_source_matrix": "yes",
                "source_identifier": gene,
                "source_row": source_row,
                "note": "one exact gene-symbol row in the deposited matrix",
            }
        )
        for sample, value in checked_values(selected.iloc[0], sample_columns).items():
            long_rows.append(
                {
                    "species": "mouse",
                    "gene": gene,
                    **metadata,
                    "source_identifier": gene,
                    "source_row": source_row,
                    "sample": sample,
                    **parse_mouse_sample(sample),
                    "normalized_value": float(value),
                    "detected": "yes" if value > 0 else "no",
                    "source_accession": "GSE114919",
                    "source_url": SOURCE_URL,
                }
            )
    return long_rows, availability


def extract_rat() -> tuple[list[dict], list[dict]]:
    """Extract exact marker-symbol rows from the rat workbook."""
    matrix = pd.read_excel(RAW / "GSE114919_Rat_normalizedcounts.xlsx")
    accession_column, symbol_column = matrix.columns[:2]
    sample_columns = list(matrix.columns[3:])
    symbols = matrix[symbol_column].fillna("").astype(str).str.upper()
    name_map = read_rat_name_map()
    long_rows: list[dict] = []
    availability: list[dict] = []

    for gene, metadata in MARKERS.items():
        selected = matrix[symbols == gene]
        if selected.empty:
            availability.append(
                {
                    "species": "rat",
                    "gene": gene,
                    **metadata,
                    "available_in_source_matrix": "no",
                    "source_identifier": "",
                    "source_row": "",
                    "note": "no exact gene-symbol row; unavailable, not treated as zero",
                }
            )
            continue
        if len(selected) != 1:
            raise ValueError(f"Expected at most one rat row for {gene}; found {len(selected)}")
        row = selected.iloc[0]
        source_row = int(selected.index[0]) + 2
        source_identifier = str(row[accession_column])
        availability.append(
            {
                "species": "rat",
                "gene": gene,
                **metadata,
                "available_in_source_matrix": "yes",
                "source_identifier": source_identifier,
                "source_row": source_row,
                "note": "one exact gene-symbol row in the deposited matrix",
            }
        )
        for sample, value in checked_values(row, sample_columns).items():
            long_rows.append(
                {
                    "species": "rat",
                    "gene": gene,
                    **metadata,
                    "source_identifier": source_identifier,
                    "source_row": source_row,
                    "sample": sample,
                    **parse_rat_sample(sample, name_map),
                    "normalized_value": float(value),
                    "detected": "yes" if value > 0 else "no",
                    "source_accession": "GSE114919",
                    "source_url": SOURCE_URL,
                }
            )
    return long_rows, availability


def summarize_conditions(long_rows: list[dict]) -> list[dict]:
    """Compare HZ and PZ means only within one species, bone and age."""
    grouped: dict[tuple, list[float]] = defaultdict(list)
    for row in long_rows:
        key = (
            row["species"],
            row["bone"],
            row["age_weeks"],
            row["gene"],
            row["zone"],
        )
        grouped[key].append(float(row["normalized_value"]))

    conditions = sorted(
        {(r["species"], r["bone"], r["age_weeks"], r["gene"]) for r in long_rows}
    )
    rows: list[dict] = []
    for species, bone, age, gene in conditions:
        pz = grouped[(species, bone, age, gene, "PZ")]
        hz = grouped[(species, bone, age, gene, "HZ")]
        pz_mean = statistics.fmean(pz)
        hz_mean = statistics.fmean(hz)
        expected = MARKERS[gene]["expected_zone"]
        difference = hz_mean - pz_mean
        rows.append(
            {
                "species": species,
                "bone": bone,
                "age_weeks": age,
                "gene": gene,
                **MARKERS[gene],
                "pz_replicates": len(pz),
                "hz_replicates": len(hz),
                "pz_detected_fraction": round(sum(x > 0 for x in pz) / len(pz), 6),
                "hz_detected_fraction": round(sum(x > 0 for x in hz) / len(hz), 6),
                "pz_mean_published_value": round(pz_mean, 6),
                "hz_mean_published_value": round(hz_mean, 6),
                "hz_minus_pz_published_value": round(difference, 6),
                "observed_higher_zone": "HZ" if difference > 0 else "PZ" if difference < 0 else "equal",
                "matches_expected_direction": (
                    "not_applicable" if expected == "none" else "yes" if difference > 0 else "no"
                ),
            }
        )
    return rows


def summarize_markers(condition_rows: list[dict], availability: list[dict]) -> list[dict]:
    """Condense condition-level checks without pooling mouse and rat values."""
    available_lookup = {(r["species"], r["gene"]): r for r in availability}
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in condition_rows:
        grouped[(row["species"], row["gene"])].append(row)

    rows: list[dict] = []
    for species in ("mouse", "rat"):
        for gene, metadata in MARKERS.items():
            source = available_lookup[(species, gene)]
            conditions = grouped.get((species, gene), [])
            differences = [float(r["hz_minus_pz_published_value"]) for r in conditions]
            expected = metadata["expected_zone"]
            rows.append(
                {
                    "species": species,
                    "gene": gene,
                    **metadata,
                    "available_in_source_matrix": source["available_in_source_matrix"],
                    "condition_count": len(conditions),
                    "conditions_higher_in_hz": sum(x > 0 for x in differences),
                    "conditions_higher_in_pz": sum(x < 0 for x in differences),
                    "conditions_matching_expected_direction": (
                        "not_applicable"
                        if expected == "none"
                        else sum(x > 0 for x in differences)
                    ),
                    "mean_hz_minus_pz_descriptive_only": (
                        round(statistics.fmean(differences), 6) if differences else ""
                    ),
                    "source_note": source["note"],
                }
            )
    return rows


def plot_heatmap(condition_rows: list[dict]) -> None:
    """Show within-condition HZ-minus-PZ differences for each species."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    genes = list(MARKERS)
    species_order = ("mouse", "rat")
    conditions_by_species = {
        species: sorted(
            {(r["bone"], int(r["age_weeks"])) for r in condition_rows if r["species"] == species},
            key=lambda x: (x[0] != "phalanx", x[1]),
        )
        for species in species_order
    }
    lookup = {
        (r["species"], r["bone"], int(r["age_weeks"]), r["gene"]): float(
            r["hz_minus_pz_published_value"]
        )
        for r in condition_rows
    }
    finite = list(lookup.values())
    limit = max(abs(min(finite)), abs(max(finite))) if finite else 1.0
    limit = max(limit, 0.1)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2), sharey=True)
    image = None
    for axis, species in zip(axes, species_order):
        conditions = conditions_by_species[species]
        matrix = np.full((len(genes), len(conditions)), np.nan)
        for row_index, gene in enumerate(genes):
            for column_index, (bone, age) in enumerate(conditions):
                matrix[row_index, column_index] = lookup.get(
                    (species, bone, age, gene), np.nan
                )
        image = axis.imshow(
            matrix,
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            aspect="auto",
        )
        axis.set_title(species.capitalize())
        axis.set_xticks(range(len(conditions)))
        axis.set_xticklabels(
            [f"{bone}\n{age} week" for bone, age in conditions], fontsize=9
        )
        axis.set_yticks(range(len(genes)))
        axis.set_yticklabels(genes)
        for row_index in range(len(genes)):
            for column_index in range(len(conditions)):
                value = matrix[row_index, column_index]
                label = "NA" if math.isnan(value) else f"{value:+.2f}"
                color = "#555555" if math.isnan(value) else (
                    "white" if abs(value) > limit * 0.55 else "black"
                )
                axis.text(
                    column_index,
                    row_index,
                    label,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=color,
                )
        axis.spines[:].set_visible(False)

    assert image is not None
    colorbar = fig.colorbar(image, ax=axes, fraction=0.035, pad=0.04)
    colorbar.set_label("HZ mean minus PZ mean (published normalized scale)")
    fig.suptitle("Descriptive zone-marker check for GSE114919", fontsize=14, y=1.01)
    fig.text(
        0.5,
        -0.02,
        "Positive values indicate higher expression in HZ. "
        "SOX9, COL2A1 and ACAN are cartilage controls, not zone-specific tests.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 0.93, 1))
    for suffix in ("png", "svg"):
        fig.savefig(
            FIGURES / f"gse114919_zone_marker_qc.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    mouse_rows, mouse_availability = extract_mouse()
    rat_rows, rat_availability = extract_rat()
    long_rows = mouse_rows + rat_rows
    availability = mouse_availability + rat_availability
    condition_rows = summarize_conditions(long_rows)
    marker_rows = summarize_markers(condition_rows, availability)

    write_tsv(
        TABLES / "gse114919_zone_marker_availability.tsv",
        availability,
        [
            "species",
            "gene",
            "role",
            "expected_zone",
            "available_in_source_matrix",
            "source_identifier",
            "source_row",
            "note",
        ],
    )
    write_tsv(
        TABLES / "gse114919_zone_marker_long.tsv",
        long_rows,
        [
            "species",
            "gene",
            "role",
            "expected_zone",
            "source_identifier",
            "source_row",
            "sample",
            "published_sample_name",
            "bone",
            "age_weeks",
            "zone",
            "replicate",
            "normalized_value",
            "detected",
            "source_accession",
            "source_url",
        ],
    )
    write_tsv(
        TABLES / "gse114919_zone_marker_condition_summary.tsv",
        condition_rows,
        list(condition_rows[0]),
    )
    write_tsv(
        TABLES / "gse114919_zone_marker_overall_summary.tsv",
        marker_rows,
        list(marker_rows[0]),
    )
    plot_heatmap(condition_rows)

    print(f"Wrote {len(long_rows)} marker measurements")
    print(f"Wrote {len(condition_rows)} within-condition HZ/PZ comparisons")
    for row in marker_rows:
        if row["role"] == "hypertrophic_marker":
            print(
                row["species"],
                row["gene"],
                row["available_in_source_matrix"],
                row["conditions_matching_expected_direction"],
                "/",
                row["condition_count"],
            )


if __name__ == "__main__":
    main()
