#!/usr/bin/env python3
"""Run a separate full-family SLRP screen on the GSE114919 matrices.

This sensitivity analysis expands the original nine-gene comparison to the 18
canonical SLRP genes.  It does not change the seven-gene comparative panel.
Ranks are calculated separately within each species, age, bone and growth-plate
zone because the published mouse and rat values must not be compared directly.

The deposited rat matrix does not contain rows for ECM2, CHADL, PODN or PODNL1.
Those entries are recorded as unavailable rather than as zero expression.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

SOURCE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114919"

# Record SLRP class and the gene's pre-existing thesis role.
GENES = (
    ("ASPN", "I", "background"),
    ("BGN", "I", "final_panel"),
    ("DCN", "I", "final_panel"),
    ("ECM2", "I", "not_preselected"),
    ("FMOD", "II", "final_panel"),
    ("KERA", "II", "not_preselected"),
    ("LUM", "II", "final_panel"),
    ("OMD", "II", "context"),
    ("PRELP", "II", "final_panel"),
    ("EPYC", "III", "final_panel"),
    ("OGN", "III", "final_panel"),
    ("OPTC", "III", "not_preselected"),
    ("CHAD", "IV", "not_preselected"),
    ("CHADL", "IV", "not_preselected"),
    ("NYX", "IV", "not_preselected"),
    ("TSKU", "IV", "not_preselected"),
    ("PODN", "V", "not_preselected"),
    ("PODNL1", "V", "not_preselected"),
)

GENE_METADATA = {
    gene: {"slrp_class": slrp_class, "panel_role": panel_role}
    for gene, slrp_class, panel_role in GENES
}

# Record current mouse Gene IDs as an audit field; extraction uses source symbols.
MOUSE_GENE_IDS = {
    "ASPN": "66695",
    "BGN": "12111",
    "DCN": "13179",
    "ECM2": "407800",
    "FMOD": "14264",
    "KERA": "16545",
    "LUM": "17022",
    "OMD": "27047",
    "PRELP": "116847",
    "EPYC": "13516",
    "OGN": "18295",
    "OPTC": "116846",
    "CHAD": "12643",
    "CHADL": "214685",
    "NYX": "18116",
    "TSKU": "244152",
    "PODN": "242608",
    "PODNL1": "244550",
}


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write one tab-separated analysis table with a fixed column order."""
    # Write a stable TSV schema.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def rat_name_map() -> dict[int, str]:
    """Map numeric rat matrix suffixes to the published sample names."""
    # Recover rat condition names from the deposited companion map.
    path = TABLES / "gse114919_rat_sample_name_map.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {int(row["No."]): row["Sample "].strip() for row in rows}


def parse_mouse_sample(sample: str) -> dict[str, str | int]:
    """Decode age, bone, zone and replicate from a mouse matrix column."""
    # Reject mouse sample headers that do not match the deposited format.
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


def parse_rat_sample(sample: str, names: dict[int, str]) -> dict[str, str | int]:
    """Use the companion name map to decode a rat matrix column."""
    number_match = re.search(r"_S(\d+)$", sample)
    if not number_match:
        raise ValueError(f"Unrecognized rat sample header: {sample}")
    published_name = names[int(number_match.group(1))]
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


def numeric_values(row: pd.Series, sample_columns: list[str]) -> pd.Series:
    """Convert one selected gene row to numbers and reject missing cells."""
    # Reject missing/nonnumeric source cells rather than replacing them with zero.
    values = pd.to_numeric(row[sample_columns], errors="coerce")
    if values.isna().any():
        missing = ", ".join(values.index[values.isna()].tolist())
        raise ValueError(f"Nonnumeric or missing source values in: {missing}")
    return values.astype(float)


def extract_mouse() -> tuple[list[dict], list[dict]]:
    """Extract all 18 mouse SLRP symbol rows from the deposited workbook."""
    path = RAW / "GSE114919_Mouse_normalizedcounts.xlsx"
    matrix = pd.read_excel(path)
    gene_column = matrix.columns[0]
    sample_columns = list(matrix.columns[1:])
    gene_symbols = matrix[gene_column].astype(str).str.upper()

    long_rows: list[dict] = []
    availability: list[dict] = []
    for gene, slrp_class, panel_role in GENES:
        # Require one exact symbol row per mouse gene.
        selected = matrix[gene_symbols == gene]
        if len(selected) != 1:
            raise ValueError(
                f"Expected one mouse row for {gene}; found {len(selected)}"
            )
        source_index = int(selected.index[0])
        # Store the visible Excel row for source checking.
        source_row = source_index + 2
        values = numeric_values(selected.iloc[0], sample_columns)
        availability.append(
            {
                "species": "mouse",
                "gene": gene,
                "slrp_class": slrp_class,
                "panel_role": panel_role,
                "available_in_source_matrix": "yes",
                "source_identifier": gene,
                "ncbi_gene_id": MOUSE_GENE_IDS[gene],
                "source_row": source_row,
                "note": "one exact gene-symbol row in the deposited matrix",
            }
        )
        for sample, value in values.items():
            long_rows.append(
                {
                    "species": "mouse",
                    "gene": gene,
                    "slrp_class": slrp_class,
                    "panel_role": panel_role,
                    "source_identifier": gene,
                    "ncbi_gene_id": MOUSE_GENE_IDS[gene],
                    "source_row": source_row,
                    "sample": sample,
                    **parse_mouse_sample(sample),
                    "normalized_value": value,
                    "detected": "yes" if value > 0 else "no",
                    "source_accession": "GSE114919",
                    "source_url": SOURCE_URL,
                }
            )
    return long_rows, availability


def extract_rat() -> tuple[list[dict], list[dict]]:
    """Extract every canonical SLRP represented in the deposited rat matrix."""
    path = RAW / "GSE114919_Rat_normalizedcounts.xlsx"
    matrix = pd.read_excel(path)
    accession_column, symbol_column, _description_column = matrix.columns[:3]
    sample_columns = list(matrix.columns[3:])
    symbols = matrix[symbol_column].fillna("").astype(str).str.upper()
    names = rat_name_map()

    long_rows: list[dict] = []
    availability: list[dict] = []
    for gene, slrp_class, panel_role in GENES:
        selected = matrix[symbols == gene]
        if selected.empty:
            # Record a missing rat row as unavailable, not as zero.
            availability.append(
                {
                    "species": "rat",
                    "gene": gene,
                    "slrp_class": slrp_class,
                    "panel_role": panel_role,
                    "available_in_source_matrix": "no",
                    "source_identifier": "",
                    "ncbi_gene_id": "",
                    "source_row": "",
                    "note": "no matching gene-symbol row; unavailable, not treated as zero",
                }
            )
            continue
        if len(selected) != 1:
            raise ValueError(
                f"Expected at most one rat row for {gene}; found {len(selected)}"
            )

        source_index = int(selected.index[0])
        source_row = source_index + 2
        source_identifier = str(selected.iloc[0][accession_column])
        values = numeric_values(selected.iloc[0], sample_columns)
        availability.append(
            {
                "species": "rat",
                "gene": gene,
                "slrp_class": slrp_class,
                "panel_role": panel_role,
                "available_in_source_matrix": "yes",
                "source_identifier": source_identifier,
                "ncbi_gene_id": "",
                "source_row": source_row,
                "note": "one exact gene-symbol row in the deposited matrix",
            }
        )
        for sample, value in values.items():
            long_rows.append(
                {
                    "species": "rat",
                    "gene": gene,
                    "slrp_class": slrp_class,
                    "panel_role": panel_role,
                    "source_identifier": source_identifier,
                    "ncbi_gene_id": "",
                    "source_row": source_row,
                    "sample": sample,
                    **parse_rat_sample(sample, names),
                    "normalized_value": value,
                    "detected": "yes" if value > 0 else "no",
                    "source_accession": "GSE114919",
                    "source_url": SOURCE_URL,
                }
            )
    return long_rows, availability


def build_condition_summary(long_data: pd.DataFrame) -> pd.DataFrame:
    """Summarize replicates and rank genes inside each biological condition."""
    # Summarize replicates within species, bone, age, and zone.
    group_columns = ["species", "bone", "age_weeks", "zone", "gene"]
    summary = (
        long_data.groupby(group_columns, as_index=False)
        .agg(
            slrp_class=("slrp_class", "first"),
            panel_role=("panel_role", "first"),
            replicate_count=("normalized_value", "size"),
            detected_replicates=("normalized_value", lambda x: int((x > 0).sum())),
            mean_published_normalized_value=("normalized_value", "mean"),
            sd_published_normalized_value=("normalized_value", "std"),
            minimum=("normalized_value", "min"),
            maximum=("normalized_value", "max"),
        )
        .sort_values(group_columns)
        .reset_index(drop=True)
    )
    summary["detection_fraction"] = (
        summary["detected_replicates"] / summary["replicate_count"]
    )

    condition_columns = ["species", "bone", "age_weeks", "zone"]
    summary["genes_available_in_condition"] = summary.groupby(condition_columns)[
        "gene"
    ].transform("size")
    # Rank available SLRPs separately within each condition.
    summary["within_condition_rank"] = (
        summary.groupby(condition_columns)["mean_published_normalized_value"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    denominator = summary["genes_available_in_condition"] - 1
    summary["within_condition_rank_percentile"] = np.where(
        denominator > 0,
        (summary["within_condition_rank"] - 1) / denominator,
        0.0,
    )
    return summary


def build_tibia_summary(condition_summary: pd.DataFrame) -> pd.DataFrame:
    """Combine the four tibial age/zone conditions within each species."""
    # Combine tibial conditions using within-condition ranks.
    tibia = condition_summary[condition_summary["bone"] == "tibia"].copy()
    result = (
        tibia.groupby(["species", "gene"], as_index=False)
        .agg(
            slrp_class=("slrp_class", "first"),
            panel_role=("panel_role", "first"),
            tibia_conditions=("gene", "size"),
            conditions_detected_in_all_replicates=(
                "detection_fraction",
                lambda x: int((x == 1).sum()),
            ),
            mean_within_condition_rank=("within_condition_rank", "mean"),
            best_within_condition_rank=("within_condition_rank", "min"),
            worst_within_condition_rank=("within_condition_rank", "max"),
            mean_rank_percentile=("within_condition_rank_percentile", "mean"),
            rank_denominator=("genes_available_in_condition", "first"),
            mean_of_condition_means_descriptive_only=(
                "mean_published_normalized_value",
                "mean",
            ),
        )
        .sort_values(["species", "mean_rank_percentile", "gene"])
        .reset_index(drop=True)
    )
    return result


def build_overall_summary(
    tibia_summary: pd.DataFrame, availability: pd.DataFrame
) -> pd.DataFrame:
    """Place mouse and rat rank summaries side by side for easy review."""
    rows: list[dict] = []
    for gene, slrp_class, panel_role in GENES:
        row: dict[str, object] = {
            "gene": gene,
            "slrp_class": slrp_class,
            "panel_role": panel_role,
        }
        percentiles: list[float] = []
        for species in ("mouse", "rat"):
            available = availability[
                (availability["species"] == species) & (availability["gene"] == gene)
            ].iloc[0]
            present = available["available_in_source_matrix"] == "yes"
            row[f"{species}_available_in_source_matrix"] = "yes" if present else "no"
            match = tibia_summary[
                (tibia_summary["species"] == species) & (tibia_summary["gene"] == gene)
            ]
            if match.empty:
                row[f"{species}_mean_rank"] = ""
                row[f"{species}_rank_denominator"] = ""
                row[f"{species}_mean_rank_percentile"] = ""
            else:
                values = match.iloc[0]
                percentile = float(values["mean_rank_percentile"])
                row[f"{species}_mean_rank"] = round(
                    float(values["mean_within_condition_rank"]), 3
                )
                row[f"{species}_rank_denominator"] = int(values["rank_denominator"])
                row[f"{species}_mean_rank_percentile"] = round(percentile, 6)
                percentiles.append(percentile)
        row["mean_rank_percentile_available_species"] = round(
            float(np.mean(percentiles)), 6
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["mean_rank_percentile_available_species", "gene"]
    )


def plot_tibia_rank_heatmap(
    condition_summary: pd.DataFrame, overall_summary: pd.DataFrame
) -> None:
    """Plot all available tibial ranks; missing source rows remain visibly blank."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tibia = condition_summary[condition_summary["bone"] == "tibia"].copy()
    order = overall_summary["gene"].tolist()
    conditions = [
        ("mouse", 1, "PZ"),
        ("mouse", 1, "HZ"),
        ("mouse", 4, "PZ"),
        ("mouse", 4, "HZ"),
        ("rat", 1, "PZ"),
        ("rat", 1, "HZ"),
        ("rat", 4, "PZ"),
        ("rat", 4, "HZ"),
    ]
    labels = [
        "Mouse 1w PZ",
        "Mouse 1w HZ",
        "Mouse 4w PZ",
        "Mouse 4w HZ",
        "Rat 1w PZ",
        "Rat 1w HZ",
        "Rat 4w PZ",
        "Rat 4w HZ",
    ]

    values = np.full((len(order), len(conditions)), np.nan)
    annotations = np.full((len(order), len(conditions)), "NA", dtype=object)
    for row_index, gene in enumerate(order):
        for column_index, (species, age, zone) in enumerate(conditions):
            match = tibia[
                (tibia["gene"] == gene)
                & (tibia["species"] == species)
                & (tibia["age_weeks"] == age)
                & (tibia["zone"] == zone)
            ]
            if match.empty:
                continue
            record = match.iloc[0]
            values[row_index, column_index] = float(
                record["within_condition_rank_percentile"]
            )
            annotations[row_index, column_index] = (
                f"{int(record['within_condition_rank'])}/"
                f"{int(record['genes_available_in_condition'])}"
            )

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#D9D9D9")
    fig, axis = plt.subplots(figsize=(11.5, 8.2))
    image = axis.imshow(values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    axis.set_xticks(range(len(labels)), labels, rotation=35, ha="right")
    axis.set_yticks(range(len(order)), order)
    axis.axvline(3.5, color="white", linewidth=2.5)

    for label, gene in zip(axis.get_yticklabels(), order, strict=True):
        if GENE_METADATA[gene]["panel_role"] == "final_panel":
            label.set_fontweight("bold")

    for row_index in range(len(order)):
        for column_index in range(len(conditions)):
            value = values[row_index, column_index]
            text_color = (
                "#FFFFFF" if not np.isnan(value) and value < 0.48 else "#111111"
            )
            axis.text(
                column_index,
                row_index,
                annotations[row_index, column_index],
                ha="center",
                va="center",
                fontsize=8,
                color=text_color,
            )

    axis.set_title(
        "Full-family SLRP expression ranks in GSE114919 tibial growth plate",
        loc="left",
        fontsize=14,
        pad=16,
    )
    axis.set_xlabel("Condition (PZ, proliferative zone; HZ, hypertrophic zone)")
    axis.set_ylabel("SLRP gene (bold genes form the existing seven-gene panel)")
    colorbar = fig.colorbar(image, ax=axis, fraction=0.035, pad=0.025)
    colorbar.set_label(
        "Rank percentile within the available family set (lower is higher)"
    )
    fig.text(
        0.01,
        0.01,
        "Cells show rank/number of SLRPs available in that species. NA means the gene row was absent from the deposited rat matrix, not zero expression.",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))

    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        fig.savefig(
            FIGURES / f"gse114919_full_slrp_family_tibia_rank_heatmap.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    mouse_rows, mouse_availability = extract_mouse()
    rat_rows, rat_availability = extract_rat()
    long_rows = mouse_rows + rat_rows
    availability_rows = mouse_availability + rat_availability

    long_fields = [
        "species",
        "gene",
        "slrp_class",
        "panel_role",
        "source_identifier",
        "ncbi_gene_id",
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
    ]
    write_tsv(TABLES / "gse114919_full_slrp_family_long.tsv", long_rows, long_fields)

    availability_fields = [
        "species",
        "gene",
        "slrp_class",
        "panel_role",
        "available_in_source_matrix",
        "source_identifier",
        "ncbi_gene_id",
        "source_row",
        "note",
    ]
    write_tsv(
        TABLES / "gse114919_full_slrp_family_source_availability.tsv",
        availability_rows,
        availability_fields,
    )

    long_data = pd.DataFrame(long_rows)
    availability = pd.DataFrame(availability_rows)
    condition_summary = build_condition_summary(long_data)
    condition_summary.to_csv(
        TABLES / "gse114919_full_slrp_family_condition_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    tibia_summary = build_tibia_summary(condition_summary)
    tibia_summary.to_csv(
        TABLES / "gse114919_full_slrp_family_tibia_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    overall_summary = build_overall_summary(tibia_summary, availability)
    overall_summary.to_csv(
        TABLES / "gse114919_full_slrp_family_overall_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    try:
        plot_tibia_rank_heatmap(condition_summary, overall_summary)
    except ModuleNotFoundError:
        print(
            "Matplotlib is unavailable in this environment; tables were written "
            "and the heatmap can be generated with the companion plot script."
        )

    missing_rat = availability[
        (availability["species"] == "rat")
        & (availability["available_in_source_matrix"] == "no")
    ]["gene"].tolist()
    print(f"Mouse genes extracted: {len(mouse_availability)} of {len(GENES)}")
    print(
        "Rat genes extracted: "
        f"{len(rat_availability) - len(missing_rat)} of {len(GENES)}"
    )
    print(f"Rat genes unavailable in the source matrix: {', '.join(missing_rat)}")
    print("\nTibial rank summary (lower percentile means higher expression):")
    print(
        overall_summary[
            [
                "gene",
                "panel_role",
                "mouse_mean_rank",
                "rat_mean_rank",
                "mean_rank_percentile_available_species",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
