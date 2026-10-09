#!/usr/bin/env python3
"""Plot the separate full-family GSE114919 tibial-rank sensitivity check."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from summarize_gse114919_full_slrp_family import plot_tibia_rank_heatmap

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"


def main() -> None:
    # Plot the saved full-family summary tables.
    condition_summary = pd.read_csv(
        TABLES / "gse114919_full_slrp_family_condition_summary.tsv", sep="\t"
    )
    overall_summary = pd.read_csv(
        TABLES / "gse114919_full_slrp_family_overall_summary.tsv", sep="\t"
    )
    plot_tibia_rank_heatmap(condition_summary, overall_summary)
    print("Wrote the full-family tibial-rank heatmap to the figures directory")


if __name__ == "__main__":
    main()
