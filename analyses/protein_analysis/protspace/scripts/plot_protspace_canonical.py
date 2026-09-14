"""Plot canonical ProtSpace projections from the validated projection table."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
TABLE = BASE / "protspace" / "tables" / "protspace_canonical_projections.tsv"
OUT_DIR = BASE / "protspace" / "figures"


def main() -> None:
    data = pd.read_csv(TABLE, sep="\t")
    genes = ["BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM"]
    colors = {
        "BGN": "#1f77b4",
        "DCN": "#8c564b",
        "EPYC": "#ff7f0e",
        "FMOD": "#2ca02c",
        "OGN": "#d62728",
        "PRELP": "#9467bd",
        "LUM": "#e7298a",
        "Control": "#555555",
    }
    projection_names = data["projection_name"].drop_duplicates().tolist()
    fig, axes = plt.subplots(
        1, len(projection_names), figsize=(16, 5), constrained_layout=True
    )
    for ax, projection_name in zip(axes, projection_names, strict=False):
        current = data[data["projection_name"] == projection_name]
        for label in genes + ["Control"]:
            if label == "Control":
                subset = current[current["analysis_group"] != "SLRP_candidate"]
                marker = "X"
            else:
                subset = current[current["gene"] == label]
                marker = "o"
            ax.scatter(
                subset["x"],
                subset["y"],
                s=34,
                alpha=0.85,
                label=label,
                color=colors[label],
                marker=marker,
                edgecolor="white",
                linewidth=0.4,
            )
        labels = current[
            (current["analysis_group"] != "SLRP_candidate")
            | (current["embedding_qc_status"] == "Needs inspection")
        ]
        for point in labels.itertuples(index=False):
            ax.annotate(
                point.identifier,
                (point.x, point.y),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
            )
        ax.set_title(projection_name)
        ax.set_xlabel("dimension 1")
        ax.set_ylabel("dimension 2")
        ax.grid(alpha=0.15)
    axes[-1].legend(fontsize=8, loc="best")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "protspace_canonical_projections.png", dpi=220)
    fig.savefig(OUT_DIR / "protspace_canonical_projections.svg")
    plt.close(fig)
    print(f"Wrote canonical ProtSpace projections for {len(projection_names)} methods")


if __name__ == "__main__":
    main()
