#!/usr/bin/env python3
"""Plot replicate-level GSE114919 growth-plate expression values."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
INPUT = BASE / "tables/gse114919_slrp_growth_plate_long.tsv"
OUTPUT = BASE / "figures"
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN", "ASPN", "OMD")
SPECIES_COLORS = {"mouse": "#4C78A8", "rat": "#E45756"}
ZONE_MARKERS = {"PZ": "o", "HZ": "^"}


def condition_label(row: pd.Series) -> str:
    bone = "P" if row["bone"] == "phalanx" else "T"
    return f"{row['species'][0].upper()}{int(row['age_weeks'])}{bone}-{row['zone']}"


def main() -> None:
    data = pd.read_csv(INPUT, sep="\t")
    data = data[data["gene"].isin(GENES)].copy()
    data["condition"] = data.apply(condition_label, axis=1)

    conditions = []
    for species in ("mouse", "rat"):
        for bone, age in (("phalanx", 1), ("tibia", 1), ("tibia", 4)):
            for zone in ("PZ", "HZ"):
                subset = data[
                    (data["species"] == species)
                    & (data["bone"] == bone)
                    & (data["age_weeks"] == age)
                    & (data["zone"] == zone)
                ]
                if not subset.empty:
                    conditions.append(subset.iloc[0]["condition"])

    rng = np.random.default_rng(20260922)
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 10.5), sharex=True)
    for axis, gene in zip(axes.flat, GENES, strict=True):
        gene_data = data[data["gene"] == gene]
        for position, condition in enumerate(conditions):
            values = gene_data.loc[
                gene_data["condition"] == condition, "normalized_count"
            ].to_numpy()
            if not len(values):
                continue
            species = "mouse" if condition.startswith("M") else "rat"
            zone = condition.rsplit("-", 1)[1]
            jitter = rng.uniform(-0.12, 0.12, len(values))
            axis.scatter(
                position + jitter,
                values,
                s=23,
                marker=ZONE_MARKERS[zone],
                color=SPECIES_COLORS[species],
                alpha=0.78,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
            mean = values.mean()
            axis.plot(
                [position - 0.18, position + 0.18],
                [mean, mean],
                color="black",
                linewidth=1.1,
                zorder=4,
            )

        axis.set_title(gene, loc="left", fontsize=11, fontweight="bold")
        axis.grid(axis="y", alpha=0.18)
        axis.axvline(5.5, color="#777777", linewidth=0.8, linestyle="--")
        axis.set_xticks(range(len(conditions)), conditions, rotation=60, ha="right")
        axis.tick_params(axis="x", labelsize=7.5)
        axis.tick_params(axis="y", labelsize=8)

    for axis in axes[:, 0]:
        axis.set_ylabel("Published normalized-value scale")

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="white",
            label=species.capitalize(),
        )
        for species, color in SPECIES_COLORS.items()
    ]
    legend_handles.extend(
        [
            plt.Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="none",
                color="#555555",
                label="Proliferative zone" if zone == "PZ" else "Hypertrophic zone",
            )
            for zone, marker in ZONE_MARKERS.items()
        ]
    )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.985),
    )
    fig.suptitle(
        "Replicate-level SLRP expression in the GSE114919 growth-plate dataset",
        fontsize=14,
        y=1.015,
    )
    fig.text(
        0.5,
        0.006,
        "Condition code: species (M/R), age in weeks, bone (P = phalanx; T = tibia), zone. "
        "Points are biological replicates; black bars are condition means.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.95))

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        fig.savefig(
            OUTPUT / f"gse114919_slrp_replicates.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"Wrote replicate-level figure to {OUTPUT}")


if __name__ == "__main__":
    main()
