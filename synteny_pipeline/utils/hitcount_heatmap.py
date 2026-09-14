from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

DEFAULT_SPECIES_GROUPS_JSON = Path("data/assets/species_groups.json")
DEFAULT_SPECIES_ICONS_DIR = Path("data/assets")


def parse_args():
    p = argparse.ArgumentParser(
        description="Draw pairwise heatmap from pairwise_summary.tsv using hits_found or hits_passing."
    )
    p.add_argument(
        "--pairwise-summary",
        required=True,
        help="Path to *_pairwise_summary.tsv",
    )
    p.add_argument(
        "--output-png",
        required=True,
        help="Output PNG path",
    )
    p.add_argument(
        "--value-column",
        choices=["hits_found", "hits_passing", "total_neighbors", "score_fraction"],
        default="hits_passing",
        help="Which column to plot",
    )
    p.add_argument(
        "--exclude-species",
        action="append",
        default=[],
        help="Species to exclude. Can be used multiple times.",
    )
    p.add_argument(
        "--species-groups-json",
        default=str(DEFAULT_SPECIES_GROUPS_JSON),
        help="JSON file with species groups",
    )
    p.add_argument(
        "--species-icons-dir",
        default=str(DEFAULT_SPECIES_ICONS_DIR),
        help="Directory containing species icons",
    )
    p.add_argument(
        "--sort-species",
        choices=["alphabetical", "groups"],
        default="groups",
    )
    p.add_argument(
        "--title",
        default=None,
        help="Custom plot title",
    )
    return p.parse_args()


def load_species_groups(path: str | None) -> dict[str, list[str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def order_species_by_groups(
    species_list: list[str], species_groups: dict[str, list[str]]
) -> list[str]:
    present = set(species_list)
    ordered = []
    used = set()

    for _, group_species in species_groups.items():
        for sp in group_species:
            if sp in present and sp not in used:
                ordered.append(sp)
                used.add(sp)

    remaining = [sp for sp in species_list if sp not in used]
    ordered.extend(sorted(remaining))
    return ordered


def get_group_spans(
    species_order: list[str], species_groups: dict[str, list[str]]
) -> list[tuple[str, int, int]]:
    spans = []
    species_pos = {sp: i for i, sp in enumerate(species_order)}

    for group_name, group_species in species_groups.items():
        present_positions = [
            species_pos[sp] for sp in group_species if sp in species_pos
        ]
        if not present_positions:
            continue
        start = min(present_positions)
        end = max(present_positions) + 1
        spans.append((group_name, start, end))

    spans.sort(key=lambda x: x[1])
    return spans


def group_color(group_name: str) -> str:
    g = group_name.lower()

    if "semi-aquatic" in g or "water" in g or "aquatic" in g:
        return "#9ecae1"
    if "amphib" in g:
        return "#6baed6"
    if "fish" in g:
        return "#3182bd"
    if "rept" in g:
        return "#a1d99b"
    if "bird" in g:
        return "#fdd0a2"
    if "non-placental" in g:
        return "#d9b3ff"
    if "placental" in g or "mammal" in g:
        return "#f4cccc"

    return "#d9d9d9"


def get_species_icon_path(species: str, icons_dir: Path | None) -> Path | None:
    if icons_dir is None:
        return None

    candidates = [
        icons_dir / f"{species}.png",
        icons_dir / f"{species}.PNG",
        icons_dir / f"{species}.jpg",
        icons_dir / f"{species}.JPG",
        icons_dir / f"{species}.jpeg",
        icons_dir / f"{species}.JPEG",
        icons_dir / f"{species}.webp",
        icons_dir / f"{species}.WEBP",
        icons_dir / f"{species.capitalize()}.png",
        icons_dir / f"{species.capitalize()}.PNG",
    ]

    for p in candidates:
        if p.exists():
            return p
    return None


def draw_icon_preserve_aspect(
    ax,
    img,
    box_x0: float,
    box_x1: float,
    box_y0: float,
    box_y1: float,
    zorder: int = 5,
):
    img_h, img_w = img.shape[:2]
    if img_h <= 0 or img_w <= 0:
        return

    box_w = box_x1 - box_x0
    box_h = box_y1 - box_y0
    img_ratio = img_w / img_h
    box_ratio = box_w / box_h

    if img_ratio >= box_ratio:
        draw_w = box_w
        draw_h = box_w / img_ratio
    else:
        draw_h = box_h
        draw_w = box_h * img_ratio

    cx = (box_x0 + box_x1) / 2
    cy = (box_y0 + box_y1) / 2

    x0 = cx - draw_w / 2
    x1 = cx + draw_w / 2
    y0 = cy - draw_h / 2
    y1 = cy + draw_h / 2

    ax.imshow(img, extent=(x0, x1, y0, y1), aspect="auto", zorder=zorder)


def main():
    args = parse_args()

    df = pd.read_csv(args.pairwise_summary, sep="\t")

    exclude_species = set(args.exclude_species)

    if exclude_species:
        df = df[
            (~df["reference_species"].isin(exclude_species))
            & (~df["target_species"].isin(exclude_species))
        ].copy()

    species_groups = load_species_groups(args.species_groups_json)
    icons_dir = Path(args.species_icons_dir) if args.species_icons_dir else None

    pivot = df.pivot(
        index="reference_species",
        columns="target_species",
        values=args.value_column,
    )

    row_species = list(pivot.index)
    col_species = list(pivot.columns)
    all_species = sorted(set(row_species) | set(col_species))

    if args.sort_species == "groups" and species_groups:
        ordered = order_species_by_groups(all_species, species_groups)
    else:
        ordered = sorted(all_species)

    row_order = [sp for sp in ordered if sp in pivot.index]
    col_order = [sp for sp in ordered if sp in pivot.columns]

    pivot = pivot.loc[row_order, col_order]
    mat = pivot.values

    n_rows = len(row_order)
    n_cols = len(col_order)

    show_groups = bool(species_groups) and args.sort_species == "groups"
    show_icons = icons_dir is not None

    left_margin = 2.2
    if show_icons:
        left_margin += 0.9
    if show_groups:
        left_margin += 1.4

    fig_w = max(8, left_margin + n_cols * 0.45 + 2.0)
    fig_h = max(6, 1.5 + n_rows * 0.45)

    cmap = plt.cm.get_cmap("viridis").copy()
    cmap.set_bad("#e6e6e6")

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=180)

    # red color scale: light = low conservation, dark red = high conservation
    cmap = matplotlib.colormaps["Reds"].copy()
    cmap.set_bad(color="#f2f2f2")  # missing values shown in light grey

    im = ax.imshow(
        mat, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0, interpolation="nearest"
    )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_order, rotation=90, fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_order, fontsize=8)

    ax.set_xlabel("Target species", fontsize=10)
    ax.set_ylabel("Reference species", fontsize=10)
    ax.set_title(args.title, fontsize=12, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Fraction of conserved neighbors", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    if show_groups:
        spans = get_group_spans(row_order, species_groups)
        bar_x0 = -left_margin + 0.2
        bar_w = 0.25
        label_x = bar_x0 - 0.1

        for group_name, start_idx, end_idx in spans:
            color = group_color(group_name)
            y0 = start_idx - 0.5
            h = end_idx - start_idx
            ax.add_patch(
                Rectangle(
                    (bar_x0, y0),
                    bar_w,
                    h,
                    facecolor=color,
                    edgecolor="none",
                    clip_on=False,
                )
            )
            ax.text(
                label_x,
                start_idx + h / 2 - 0.5,
                group_name.replace("_", " "),
                ha="right",
                va="center",
                fontsize=8,
                fontweight="bold",
                clip_on=False,
            )
            if end_idx < n_rows:
                ax.hlines(end_idx - 0.5, -0.5, n_cols - 0.5, colors=color, linewidth=2)

    if show_icons:
        icon_x0 = -0.95
        icon_x1 = -0.15
        for i, sp in enumerate(row_order):
            icon_path = get_species_icon_path(sp, icons_dir)
            if icon_path is not None:
                try:
                    img = mpimg.imread(icon_path)
                    draw_icon_preserve_aspect(
                        ax, img, icon_x0, icon_x1, i - 0.42, i + 0.42, zorder=5
                    )
                except Exception:
                    ax.text(
                        icon_x0,
                        i,
                        sp,
                        ha="right",
                        va="center",
                        fontsize=7,
                        clip_on=False,
                    )
            else:
                ax.text(
                    icon_x0, i, sp, ha="right", va="center", fontsize=7, clip_on=False
                )

    plt.tight_layout()
    plt.savefig(args.output_png, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote heatmap: {args.output_png}")
    print(f"Plotted value column: {args.value_column}")


if __name__ == "__main__":
    main()
