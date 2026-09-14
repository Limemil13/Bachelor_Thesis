"""Validate and summarize the canonical ProtSpace/ProtT5 result.

Outlier calls use cosine distances in the original 1024-dimensional ProtT5
embedding, not visual distances in PCA/UMAP/t-SNE. The 2-D projections are
generated only for exploratory figures.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE = Path(__file__).resolve().parents[2]
INPUT_ANNOTATIONS = (
    BASE / "protspace" / "input" / "slrp_annotations_canonical_with_controls.csv"
)
RUN_DIR = (
    BASE
    / "protspace"
    / "output"
    / "slrp_protspace_seven_gene_curated_with_controls_20260903"
)
BUNDLE = RUN_DIR / "data.parquetbundle"
EMBEDDINGS = RUN_DIR / "tmp" / "prot_t5.h5"
TABLE_DIR = BASE / "protspace" / "tables"
FIGURE_DIR = BASE / "protspace" / "figures"
DELIMITER = b"---PARQUET_DELIMITER---"


def read_bundle_tables(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    parts = path.read_bytes().split(DELIMITER)
    if len(parts) not in (3, 4):
        raise ValueError(f"Expected 3 or 4 ProtSpace bundle parts, found {len(parts)}")
    frames = [pq.read_table(io.BytesIO(part)).to_pandas() for part in parts[:3]]
    return frames[0], frames[1], frames[2]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def robust_z(values: np.ndarray) -> np.ndarray:
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros_like(values)
    return 0.67448975 * (values - median) / mad


def silhouette_from_distance(distance: np.ndarray, labels: list[str]) -> np.ndarray:
    result = np.zeros(len(labels), dtype=float)
    for i, label in enumerate(labels):
        same = [j for j, other in enumerate(labels) if other == label and j != i]
        a = float(np.mean(distance[i, same]))
        other_means = []
        for other_label in sorted(set(labels) - {label}):
            members = [j for j, item in enumerate(labels) if item == other_label]
            other_means.append(float(np.mean(distance[i, members])))
        b = min(other_means)
        result[i] = (b - a) / max(a, b) if max(a, b) else 0.0
    return result


def main() -> None:
    annotations = pd.read_csv(INPUT_ANNOTATIONS, dtype=str).fillna("")
    expected_total = 107  # 103 canonical SLRPs plus four retained controls.
    if (
        len(annotations) != expected_total
        or annotations["identifier"].nunique() != expected_total
    ):
        raise ValueError(
            f"Expected {expected_total} unique canonical ProtSpace annotations"
        )

    bundle_annotations, projection_metadata, projection_data = read_bundle_tables(
        BUNDLE
    )
    canonical_ids = set(annotations["identifier"])
    if len(projection_metadata) != 3 or not canonical_ids.issubset(
        set(bundle_annotations["protein_id"])
    ):
        raise ValueError(
            "Canonical IDs are not a subset of the validated ProtSpace bundle"
        )
    # The retained bundle contains one additional diagnostic embedding for the
    # excluded ASPN-like chicken BGN-locus protein. Quantitative distances are
    # recomputed below using only the 107 retained IDs; excluded projection rows are
    # also removed before export/plotting.
    projection_data = projection_data[
        projection_data["identifier"].isin(canonical_ids)
    ].copy()
    if len(projection_data) != expected_total * 3:
        raise ValueError("Unexpected canonical ProtSpace projection subset dimensions")

    identifiers = annotations["identifier"].tolist()
    with h5py.File(EMBEDDINGS, "r") as handle:
        if not set(identifiers).issubset(set(handle.keys())):
            raise ValueError(
                "Canonical embedding IDs are missing from the validated bundle"
            )
        vectors = np.vstack(
            [handle[identifier][:] for identifier in identifiers]
        ).astype(float)
    if vectors.shape != (expected_total, 1024) or not np.isfinite(vectors).all():
        raise ValueError(f"Unexpected embedding matrix: {vectors.shape}")

    normalized = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    distance = np.clip(1.0 - normalized @ normalized.T, 0.0, 2.0)

    slrp_mask = annotations["analysis_group"].eq("SLRP_candidate").to_numpy()
    slrp_indices = np.flatnonzero(slrp_mask)
    control_indices = np.flatnonzero(~slrp_mask)
    slrp_labels = annotations.loc[slrp_mask, "gene"].tolist()
    slrp_distance = distance[np.ix_(slrp_indices, slrp_indices)]
    silhouettes = silhouette_from_distance(slrp_distance, slrp_labels)

    preliminary: list[dict[str, object]] = []
    median_within_by_gene: dict[str, list[tuple[int, float]]] = {}
    genes = sorted(set(slrp_labels))
    for local_i, global_i in enumerate(slr_indices := slrp_indices.tolist()):
        row = annotations.iloc[global_i]
        gene = row["gene"]
        same = [
            idx
            for idx in slr_indices
            if annotations.iloc[idx]["gene"] == gene and idx != global_i
        ]
        other = [idx for idx in slr_indices if annotations.iloc[idx]["gene"] != gene]
        nearest_same = min(same, key=lambda idx: distance[global_i, idx])
        nearest_other = min(other, key=lambda idx: distance[global_i, idx])
        nearest_control = min(control_indices, key=lambda idx: distance[global_i, idx])

        centroid_distances: dict[str, float] = {}
        for candidate_gene in genes:
            members = [
                idx
                for idx in slr_indices
                if annotations.iloc[idx]["gene"] == candidate_gene and idx != global_i
            ]
            centroid = normalized[members].mean(axis=0)
            centroid /= np.linalg.norm(centroid)
            centroid_distances[candidate_gene] = float(
                1.0 - normalized[global_i] @ centroid
            )
        predicted_gene = min(centroid_distances, key=centroid_distances.get)
        median_within = float(np.median(distance[global_i, same]))
        median_within_by_gene.setdefault(gene, []).append((global_i, median_within))

        preliminary.append(
            {
                "identifier": row["identifier"],
                "gene": gene,
                "species_accession": row["old_header"].split()[0],
                "sequence_correction": row["sequence_correction"],
                "nearest_same_gene": annotations.iloc[nearest_same]["identifier"],
                "nearest_same_gene_cosine_distance": round(
                    float(distance[global_i, nearest_same]), 6
                ),
                "nearest_other_gene": annotations.iloc[nearest_other]["identifier"],
                "nearest_other_gene_label": annotations.iloc[nearest_other]["gene"],
                "nearest_other_gene_cosine_distance": round(
                    float(distance[global_i, nearest_other]), 6
                ),
                "nearest_gene_separation_margin": round(
                    float(
                        distance[global_i, nearest_other]
                        - distance[global_i, nearest_same]
                    ),
                    6,
                ),
                "nearest_control": annotations.iloc[nearest_control]["identifier"],
                "nearest_control_cosine_distance": round(
                    float(distance[global_i, nearest_control]), 6
                ),
                "median_within_gene_cosine_distance": round(median_within, 6),
                "leave_one_out_centroid_gene": predicted_gene,
                "leave_one_out_centroid_correct": "yes"
                if predicted_gene == gene
                else "no",
                "expected_gene_centroid_distance": round(centroid_distances[gene], 6),
                "best_other_gene_centroid_distance": round(
                    min(
                        value
                        for key, value in centroid_distances.items()
                        if key != gene
                    ),
                    6,
                ),
                "silhouette_cosine": round(float(silhouettes[local_i]), 6),
            }
        )

    robust_scores: dict[int, float] = {}
    for values in median_within_by_gene.values():
        indices = [index for index, _ in values]
        scores = robust_z(np.array([value for _, value in values]))
        robust_scores.update(dict(zip(indices, scores, strict=True)))

    rows: list[dict[str, object]] = []
    for preliminary_row in preliminary:
        global_i = identifiers.index(str(preliminary_row["identifier"]))
        score = float(robust_scores[global_i])
        concerns: list[str] = []
        if preliminary_row["leave_one_out_centroid_correct"] == "no":
            concerns.append("nearest_centroid_is_another_gene")
        if score > 3.5:
            concerns.append("high_within_gene_embedding_outlier_score")
        elif score > 2.5:
            concerns.append("moderate_within_gene_embedding_outlier_score")
        if float(preliminary_row["nearest_gene_separation_margin"]) < 0:
            concerns.append("nearest_sequence_is_another_gene")
        # A negative per-gene silhouette is descriptive, not by itself a QC
        # failure: EPYC and OGN are both class III SLRPs and may be close in a
        # protein-language-model space even when their orthology is correct.

        if "nearest_centroid_is_another_gene" in concerns or score > 3.5:
            status = "Needs inspection"
        elif concerns:
            status = "Watch"
        else:
            status = "Supported"
        rows.append(
            {
                **preliminary_row,
                "within_gene_robust_z": round(score, 4),
                "embedding_qc_status": status,
                "embedding_qc_concerns": ";".join(concerns),
            }
        )

    write_tsv(TABLE_DIR / "protspace_canonical_embedding_qc.tsv", rows)
    write_tsv(
        TABLE_DIR / "protspace_canonical_embedding_review.tsv",
        [row for row in rows if row["embedding_qc_status"] != "Supported"],
    )

    gene_rows: list[dict[str, object]] = []
    for gene in genes:
        gene_rows_subset = [row for row in rows if row["gene"] == gene]
        gene_rows.append(
            {
                "gene": gene,
                "protein_count": len(gene_rows_subset),
                "centroid_assignment_correct_count": sum(
                    row["leave_one_out_centroid_correct"] == "yes"
                    for row in gene_rows_subset
                ),
                "mean_silhouette_cosine": round(
                    float(
                        np.mean(
                            [
                                float(row["silhouette_cosine"])
                                for row in gene_rows_subset
                            ]
                        )
                    ),
                    4,
                ),
                "mean_median_within_gene_cosine_distance": round(
                    float(
                        np.mean(
                            [
                                float(row["median_within_gene_cosine_distance"])
                                for row in gene_rows_subset
                            ]
                        )
                    ),
                    4,
                ),
                "supported_count": sum(
                    row["embedding_qc_status"] == "Supported"
                    for row in gene_rows_subset
                ),
                "watch_count": sum(
                    row["embedding_qc_status"] == "Watch" for row in gene_rows_subset
                ),
                "needs_inspection_count": sum(
                    row["embedding_qc_status"] == "Needs inspection"
                    for row in gene_rows_subset
                ),
            }
        )
    write_tsv(TABLE_DIR / "protspace_canonical_gene_summary.tsv", gene_rows)

    projection_names = projection_metadata["projection_name"].tolist()
    projection_export = projection_data.merge(
        annotations[["identifier", "gene", "analysis_group"]],
        on="identifier",
        how="left",
    )
    projection_export["embedding_qc_status"] = (
        projection_export["identifier"]
        .map({str(row["identifier"]): str(row["embedding_qc_status"]) for row in rows})
        .fillna("Control")
    )
    projection_export.to_csv(
        TABLE_DIR / "protspace_canonical_projections.tsv", sep="\t", index=False
    )

    print(
        f"Validated canonical subset: {expected_total} proteins from a "
        f"{len(bundle_annotations)}-protein embedding bundle, {len(projection_names)} projections"
    )
    print(f"Embedding matrix: {vectors.shape}")
    print(
        f"Supported: {sum(row['embedding_qc_status'] == 'Supported' for row in rows)}"
    )
    print(f"Watch: {sum(row['embedding_qc_status'] == 'Watch' for row in rows)}")
    print(
        f"Needs inspection: {sum(row['embedding_qc_status'] == 'Needs inspection' for row in rows)}"
    )


if __name__ == "__main__":
    main()
