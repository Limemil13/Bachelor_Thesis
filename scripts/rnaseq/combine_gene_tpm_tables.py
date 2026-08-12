import pandas as pd
from pathlib import Path

files = {
    "growth_plate_WT1": "results/WT1_gene_tpm.tsv",
    "growth_plate_WT2": "results/WT2_gene_tpm.tsv",
    "growth_plate_WT3": "results/WT3_gene_tpm.tsv",
    "skin_1": "results/skin_1_gene_tpm.tsv",
    "heart_1": "results/heart_1_gene_tpm.tsv",
    "kidney_1": "results/kidney_1_gene_tpm.tsv",
    "brain_1": "results/brain_1_gene_tpm.tsv",
}

tables = []

for sample, file in files.items():
    df = pd.read_csv(file, sep="\t")
    df = df[["gene_symbol", "total_TPM"]]
    df = df.rename(columns={"total_TPM": sample})
    tables.append(df)

combined = tables[0]
for df in tables[1:]:
    combined = combined.merge(df, on="gene_symbol", how="outer")

gp_cols = ["growth_plate_WT1", "growth_plate_WT2", "growth_plate_WT3"]
combined["growth_plate_mean"] = combined[gp_cols].mean(axis=1)

other_cols = ["skin_1", "heart_1", "kidney_1", "brain_1"]
combined["other_tissue_mean"] = combined[other_cols].mean(axis=1)

combined["growth_plate_vs_other_ratio"] = (
    combined["growth_plate_mean"] / (combined["other_tissue_mean"] + 0.01)
)

combined.to_csv("results/combined_mouse_tissue_expression.tsv", sep="\t", index=False)

print(combined)
