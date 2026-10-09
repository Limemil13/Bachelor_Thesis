"""description of expression of the three wild type samples
reads salmon quant.sf files and metadata that was corrected
sums transcript TPM by gene, calculates mean, SD, coef var and rank
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
METADATA = BASE / "sample_metadata_corrected.tsv"
TABLE_DIR = BASE / "tables"
FIGURE_DIR = BASE / "figures"

GENES_TO_PLOT = [
    "Col2a1",
    "Acan",
    "Sox9",
    "Comp",
    "Matn1",
    "Col9a1",
    "Col9a2",
    "Col11a1",
    "Bgn",
    "Fmod",
    "Epyc",
    "Ogn",
    "Omd",
    "Prelp",
    "Lum",
    "Aspn",
    "Dcn",
    "Chad",
    "Chadl",
    "Tsku",
    "Podn",
    "Podnl1",
    "Kera",
    "Ecm2",
]


def main() -> None:
    #source directory has for each sample one modified salmon gene TPM table
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-root",
        required=True,
        type=Path,
        help="RNA-seq workspace containing results/*_gene_tpm.tsv",
    )
    args = parser.parse_args()


    #read id and TPM val from every sample
    metadata = pd.read_csv(METADATA, sep="\t")
    tables: list[pd.DataFrame] = []


    for row in metadata.itertuples(index=False):
        local_name = str(row.local_quantification).replace("_quant", "")
        source = args.source_root / "results" / f"{local_name}_gene_tpm.tsv"
        if not source.exists():
            raise FileNotFoundError(source)
        table = pd.read_csv(source, sep="\t")[["gene_symbol", "total_TPM"]]
        tables.append(table.rename(columns={"total_TPM": row.analysis_column}))


    combined = tables[0]
    for table in tables[1:]:
        combined = combined.merge(
            table, on="gene_symbol", how="outer", validate="one_to_one"
        )


    focal = ["P0_WT1_paired", "P0_WT2_paired", "P0_WT3_single"]
    combined["p0_chondroprogenitor_mean_TPM"] = combined[focal].mean(axis=1)
    combined["p0_chondroprogenitor_sd_TPM"] = combined[focal].std(axis=1, ddof=1)
    combined["p0_chondroprogenitor_CV"] = np.where(
        combined["p0_chondroprogenitor_mean_TPM"] > 0,
        combined["p0_chondroprogenitor_sd_TPM"]
        / combined["p0_chondroprogenitor_mean_TPM"],
        np.nan,
    )
    combined["p0_detected_in_all_three"] = (combined[focal] > 0).all(axis=1)


    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    table_out = TABLE_DIR / "mouse_expression_descriptive_metadata_corrected.tsv"
    combined.to_csv(table_out, sep="\t", index=False)

    #need to leave pseudocount out
    plot = combined[combined["gene_symbol"].isin(GENES_TO_PLOT)].copy()
    plot["gene_symbol"] = pd.Categorical(
        plot["gene_symbol"], categories=GENES_TO_PLOT, ordered=True
    )
    plot = plot.sort_values("gene_symbol")
    columns = metadata["analysis_column"].tolist()
    matrix = plot.set_index("gene_symbol")[columns]

    labels = [
        "P0 WT1\npaired",
        "P0 WT2\npaired",
        "P0 WT3\nsingle",
        "Skin\n6 mo",
        "Heart\n30 mo",
        "Kidney\n30 mo",
        "Skin\n15 mo",
    ]
    for values, name, colorbar in (
        (matrix, "mouse_expression_descriptive_raw_TPM.png", "TPM"),
        (
            np.log10(matrix + 1),
            "mouse_expression_descriptive_log10_TPM.png",
            "log10(TPM + 1)",
        ),
    ):
        fig, ax = plt.subplots(figsize=(11, 10))
        image = ax.imshow(values.values, aspect="auto")
        ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(values.index)), labels=values.index)
        ax.set_xlabel("Sample (age and layout shown where relevant)")
        ax.set_ylabel("Gene")
        ax.set_title(
            "Descriptive mouse TPM panel\n"
            "P0 chondroprogenitors and age-/study-confounded tissue comparators"
        )
        fig.colorbar(image, ax=ax, label=colorbar)
        fig.tight_layout()
        fig.savefig(FIGURE_DIR / name, dpi=300)
        plt.close(fig)

    print(f"Wrote {table_out} ({len(combined)} genes)")
    print(f"Wrote figures to {FIGURE_DIR}")
    print("No cross-study tissue-enrichment ratio was calculated.")


if __name__ == "__main__":
    main()
