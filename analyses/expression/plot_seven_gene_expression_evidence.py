#!/usr/bin/env python3
"""Plot the non-pooled GSE305415 and growth-plate expression evidence.

Panel A uses GSE305415 WT P0 Prx1-lineage TPM (mean +/- SD). Panel B uses only
within-species ranks from the independent zone-resolved GSE114919 tibial data.
The two scales are deliberately not combined.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
TABLES = BASE / "tables"
GSE = BASE / "gse114919" / "tables"
FIGURES = BASE / "figures"

GENES = ["BGN", "FMOD", "PRELP", "EPYC", "LUM", "OGN", "DCN"]
MAIN_COLORS = {
    "BGN": "#4c78a8",
    "FMOD": "#f58518",
    "PRELP": "#eeca3b",
    "EPYC": "#b279a2",
    "LUM": "#54a24b",
    "OGN": "#e45756",
    "DCN": "#9d9d9d",
}


def main() -> None:
    # Panels A and B come from separate studies and stay on separate scales. The
    # script deliberately does not calculate a combined expression score.
    local = pd.read_csv(
        TABLES / "mouse_expression_descriptive_metadata_corrected.tsv", sep="\t"
    )
    local["gene"] = local["gene_symbol"].str.upper()
    local = local.set_index("gene").loc[GENES]

    cross = pd.read_csv(
        GSE / "gse114919_slrp_tibia_cross_condition_summary.tsv", sep="\t"
    )
    # Reshape the long growth-plate table into gene-by-species ranks.
    ranks = (
        cross[cross["gene"].isin(GENES)]
        .pivot(index="gene", columns="species", values="mean_within_condition_rank")
        .loc[GENES, ["mouse", "rat"]]
    )

    # The wider first panel leaves room for bars and error bars; the second is a
    # compact heat map where rank 1 denotes the highest expression.
    fig, axes = plt.subplots(
        1, 2, figsize=(11.5, 5.7), gridspec_kw={"width_ratios": [1.55, 1]}
    )

    x = np.arange(len(GENES))
    means = local["p0_chondroprogenitor_mean_TPM"].to_numpy()
    sds = local["p0_chondroprogenitor_sd_TPM"].to_numpy()
    axes[0].bar(
        x,
        means,
        yerr=sds,
        capsize=3,
        color=[MAIN_COLORS[g] for g in GENES],
        edgecolor="black",
        linewidth=0.6,
    )
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, GENES, rotation=35, ha="right")
    axes[0].set_ylabel("TPM, log scale (mean +/- SD; n = 3)")
    axes[0].set_title("A  GSE305415 WT P0 Prx1-lineage cells", loc="left")
    axes[0].grid(axis="y", which="both", alpha=0.2)

    # Reverse the colour map so lower (better) numerical ranks stand out.
    image = axes[1].imshow(
        ranks.to_numpy(), cmap="viridis_r", vmin=1, vmax=9, aspect="auto"
    )
    axes[1].set_xticks([0, 1], ["Mouse", "Rat"])
    axes[1].set_yticks(np.arange(len(GENES)), GENES)
    axes[1].set_title("B  Tibial growth-plate mean rank", loc="left")
    for row in range(len(GENES)):
        for col in range(2):
            value = ranks.iloc[row, col]
            axes[1].text(
                col,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                color="white" if value <= 5 else "black",
                fontsize=9,
                fontweight="bold",
            )
    colorbar = fig.colorbar(image, ax=axes[1], fraction=0.05, pad=0.04)
    colorbar.set_label("Within-condition rank (1 = highest)")
    axes[1].text(
        0.5,
        -0.14,
        "Mean across age x zone conditions; values are not pooled across species",
        transform=axes[1].transAxes,
        fontsize=8.5,
        ha="center",
    )

    fig.suptitle(
        "Expression evidence for the seven-gene SLRP panel",
        fontsize=13,
        y=1.01,
    )
    fig.tight_layout()
    # Save vector and raster versions for manuscript use and quick viewing.
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        fig.savefig(
            FIGURES / f"seven_gene_expression_evidence.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"Wrote seven-gene expression figure to {FIGURES}")


if __name__ == "__main__":
    main()
