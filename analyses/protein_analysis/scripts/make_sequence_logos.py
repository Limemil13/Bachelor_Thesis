#!/usr/bin/env python3
"""Create canonical amino-acid sequence logos without external logo tools.

Letter height is amino-acid frequency multiplied by information content in bits.
The complete stack is additionally scaled by non-gap occupancy, preventing a
column represented by only one or two proteins from appearing fully conserved.
Coordinates are alignment columns, not human-protein residue numbers.
"""

from __future__ import annotations

from collections import Counter
from math import log, log2
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

ROOT = Path(__file__).resolve().parents[3]
ALIGNMENTS = ROOT / "analyses" / "protein_analysis" / "alignments" / "canonical"
OUT = ROOT / "analyses" / "protein_analysis" / "figures" / "sequence_logos"
GENES = ("BGN", "DCN", "EPYC", "FMOD", "LUM", "OGN", "PRELP")
AA = tuple("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA)
MAX_BITS = log2(20)
WINDOW = 50

COLORS = {
    **{aa: "#4C78A8" for aa in "AILMFWVY"},  # hydrophobic/aromatic
    **{aa: "#E45756" for aa in "KR"},  # positively charged
    **{aa: "#B279A2" for aa in "DE"},  # negatively charged
    **{aa: "#54A24B" for aa in "NQST"},  # polar
    "C": "#F58518",
    "G": "#F2CF5B",
    "P": "#9D755D",
    "H": "#72B7B2",
}

FONT = FontProperties(family="DejaVu Sans", weight="bold")
LETTER_PATHS = {aa: TextPath((0, 0), aa, size=1, prop=FONT) for aa in AA}


def read_fasta(path: Path) -> list[str]:
    sequences: list[str] = []
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    sequences.append("".join(chunks).upper())
                    chunks = []
            else:
                chunks.append(line)
    if chunks:
        sequences.append("".join(chunks).upper())
    if not sequences or len({len(seq) for seq in sequences}) != 1:
        raise ValueError(f"Not a rectangular alignment: {path}")
    return sequences


def letter_heights(column: list[str], n_sequences: int) -> list[tuple[str, float]]:
    observed = [residue for residue in column if residue in AA_SET]
    n = len(observed)
    if n == 0:
        return []
    counts = Counter(observed)
    probabilities = {aa: count / n for aa, count in counts.items()}
    entropy = -sum(p * log2(p) for p in probabilities.values())
    small_sample_correction = 19 / (2 * log(2) * n)
    information = max(0.0, MAX_BITS - entropy - small_sample_correction)
    occupancy = n / n_sequences
    return sorted(
        ((aa, p * information * occupancy) for aa, p in probabilities.items()),
        key=lambda item: item[1],
    )


def draw_letter(ax, aa: str, x: float, y: float, height: float) -> None:
    if height < 0.012:
        return
    path = LETTER_PATHS[aa]
    box = path.get_extents()
    transform = (
        Affine2D()
        .translate(-box.xmin, -box.ymin)
        .scale(0.86 / box.width, height / box.height)
        .translate(x + 0.07, y)
    )
    ax.add_patch(
        PathPatch(path, transform=transform + ax.transData, color=COLORS[aa], lw=0)
    )


def make_logo(gene: str) -> None:
    sequences = read_fasta(ALIGNMENTS / f"{gene}_canonical_aligned.faa")
    length = len(sequences[0])
    n_rows = (length + WINDOW - 1) // WINDOW
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 2.2 * n_rows), squeeze=False)

    for row, start in enumerate(range(0, length, WINDOW)):
        ax = axes[row, 0]
        end = min(start + WINDOW, length)
        for local_x, column_index in enumerate(range(start, end)):
            column = [seq[column_index] for seq in sequences]
            baseline = 0.0
            for aa, height in letter_heights(column, len(sequences)):
                draw_letter(ax, aa, local_x, baseline, height)
                baseline += height

        width = end - start
        ax.set_xlim(0, width)
        ax.set_ylim(0, MAX_BITS)
        tick_positions = list(range(0, width, 5))
        ax.set_xticks(
            [position + 0.5 for position in tick_positions],
            [str(start + position + 1) for position in tick_positions],
            fontsize=8,
        )
        ax.set_yticks((0, 2, 4))
        ax.set_ylabel("bits", fontsize=9)
        ax.grid(axis="y", alpha=0.18, linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"Alignment columns {start + 1}–{end}", fontsize=9, loc="left")
        if row == n_rows - 1:
            ax.set_xlabel("MAFFT alignment column (not human residue number)")

    fig.suptitle(
        f"{gene} amino-acid sequence logo — {len(sequences)} canonical proteins",
        fontsize=14,
        y=1.002,
    )
    fig.text(
        0.995,
        0.002,
        "Stack height = information content × non-gap occupancy",
        ha="right",
        va="bottom",
        fontsize=8,
    )
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        fig.savefig(
            OUT / f"{gene}_sequence_logo.{suffix}", dpi=300, bbox_inches="tight"
        )
    plt.close(fig)
    print(f"{gene}: {len(sequences)} sequences, {length} columns")


def main() -> None:
    for gene in GENES:
        make_logo(gene)
    print(f"Wrote sequence logos to {OUT}")


if __name__ == "__main__":
    main()
