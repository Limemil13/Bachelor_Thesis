import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Input file
input_file = "results/combined_mouse_tissue_expression.tsv"

# Output folder
output_raw = "figures/expression_heatmap_raw_TPM.png"
output_log = "figures/expression_heatmap_log10_TPM.png"

# Read table
df = pd.read_csv(input_file, sep="\t")

# Columns to plot
plot_cols = {
    "growth_plate_mean": "Growth Plate",
    "skin_1": "Skin",
    "heart_1": "Heart",
    "kidney_1": "Kidney",
    "brain_1": "Brain"
}

# Genes to show
genes_to_plot = [
    # cartilage / growth plate controls
    "Col2a1", "Acan", "Sox9", "Comp", "Matn1", "Col9a1", "Col9a2", "Col11a1",

    # SLRP candidates
    "Bgn", "Fmod", "Epyc", "Ogn", "Omd", "Prelp", "Lum", "Aspn", "Dcn",
    "Chad", "Chadl", "Tsku", "Podn", "Podnl1", "Kera", "Ecm2"
]

# Filter and order genes
df = df[df["gene_symbol"].isin(genes_to_plot)].copy()
df["gene_symbol"] = pd.Categorical(df["gene_symbol"], categories=genes_to_plot, ordered=True)
df = df.sort_values("gene_symbol")

# Build matrix
matrix = df[["gene_symbol"] + list(plot_cols.keys())].set_index("gene_symbol")
matrix = matrix.rename(columns=plot_cols)


def plot_heatmap(data, output_file, title, colorbar_label, number_format):
    fig, ax = plt.subplots(figsize=(10, 12))

    im = ax.imshow(data.values, aspect="auto")

    # Axis labels
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(data.columns)
    ax.set_yticklabels(data.index)

    # Rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Add values inside cells
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(
                j,
                i,
                format(data.values[i, j], number_format),
                ha="center",
                va="center",
                fontsize=7,
                color="white"
            )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(colorbar_label)

    ax.set_title(title)
    ax.set_xlabel("Tissue")
    ax.set_ylabel("Gene")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


# Raw TPM heatmap
plot_heatmap(
    matrix,
    output_raw,
    "Expression of Cartilage Markers and SLRP Genes Across Mouse Tissues",
    "TPM (Transcripts Per Million)",
    ".1f"
)

# Log10 heatmap
matrix_log = np.log10(matrix + 1)

plot_heatmap(
    matrix_log,
    output_log,
    "Expression of Cartilage Markers and SLRP Genes Across Mouse Tissues\nlog10(TPM + 1)",
    "log10(TPM + 1)",
    ".2f"
)

print("Saved:")
print(output_raw)
print(output_log)
