from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle

DEFAULT_MIN_PIDENT = 30.0
DEFAULT_MAX_EVALUE = 1e-5
DEFAULT_MIN_LEN = 60
DEFAULT_MIN_COVERAGE = 0.05

DEFAULT_EXCLUDE_GENES = {
    "CCDST",
    "CRCT1",
}

SFTP_GENES = {
    "HRNR",
    "FLG",
    "FLG2",
    "CRNN",
    "RPTN",
    "TCHH",
    "TCHHL1",
}

REDUCED_GENE_SET = [
    "CELF3",
    "RIIAD1",
    "MRPL9",
    "OAZ3",
    "TDRKH",
    "LINGO4",
    "RORC",
    "C2CD4D",
    "THEM5",
    "THEM4",
    "S100A10",
    "S100A11",
    "HRNR",
    "FLG",
    "FLG2",
    "CRNN",
    "RPTN",
    "TCHH",
    "TCHHL1",
    "LCE3E",
    "LCE3D",
    "LCE3C",
    "LCE3B",
    "LCE3A",
    "LCE2D",
    "LCE2C",
    "KPLCE",
]

DEFAULT_SPECIES_GROUPS_JSON = Path("data/assets/species_groups.json")
DEFAULT_SPECIES_ICONS_DIR = Path("data/assets")


def parse_args():
    p = argparse.ArgumentParser(
        description="Build all-vs-all consensus tables from mapping.json files."
    )
    p.add_argument("--run-folder", action="append", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--output-prefix", required=True)

    p.add_argument("--min-pident", type=float, default=DEFAULT_MIN_PIDENT)
    p.add_argument("--max-evalue", type=float, default=DEFAULT_MAX_EVALUE)
    p.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN)
    p.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)

    p.add_argument(
        "--require-reciprocal-for-weak",
        action="store_true",
        help="Require reciprocal pass for weak hits if reciprocal results exist.",
    )
    p.add_argument("--weak-pident-threshold", type=float, default=45.0)

    p.add_argument("--exclude-loc", action="store_true")
    p.add_argument("--exclude-gene", action="append", default=[])

    p.add_argument(
        "--collapse-sftp",
        action="store_true",
        help="Collapse SFTP family into one SFTP_cluster column in consensus tables.",
    )

    p.add_argument(
        "--sort-species",
        choices=["score", "alphabetical", "groups"],
        default="groups",
    )
    p.add_argument(
        "--species-groups-json",
        default=str(DEFAULT_SPECIES_GROUPS_JSON),
        help="JSON file defining species group order.",
    )
    p.add_argument(
        "--species-icons-dir",
        default=str(DEFAULT_SPECIES_ICONS_DIR),
        help="Directory containing PNG/JPG/WebP icons named like homo_sapiens.png",
    )

    p.add_argument("--pairwise-title", default=None)
    p.add_argument("--full-consensus-title", default=None)
    p.add_argument("--reduced-consensus-title", default=None)

    return p.parse_args()


def load_species_groups(path: str | None) -> dict[str, list[str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def load_reciprocal_lookup(mapping_json: Path) -> dict[str, dict]:
    rec_json = mapping_json.parent / "reciprocal" / "reciprocal.json"
    if not rec_json.exists():
        return {}
    try:
        rows = json.loads(rec_json.read_text(encoding="utf-8"))
    except Exception:
        return {}

    lookup = {}
    for row in rows:
        neighbor = row.get("neighbor")
        if neighbor:
            lookup[neighbor] = row
    return lookup


def should_exclude_gene(
    gene: str | None, exclude_loc: bool, exclude_genes: set[str]
) -> bool:
    if gene is None:
        return True
    g = gene.upper()
    if exclude_loc and g.startswith("LOC"):
        return True
    if g in exclude_genes:
        return True
    return False


def collect_mapping_files(run_folders: Iterable[Path], focal_gene: str) -> list[Path]:
    mapping_files: list[Path] = []
    for run_folder in run_folders:
        if not run_folder.exists():
            continue
        mapping_files.extend(
            run_folder.glob(f"neighbor_mapping/*/{focal_gene}/mapping.json")
        )
    return sorted(mapping_files)


def classify_significant(
    entry: dict,
    reciprocal_lookup: dict[str, dict],
    min_pident: float,
    max_evalue: float,
    min_len: int,
    min_coverage: float,
    require_reciprocal_for_weak: bool,
    weak_pident_threshold: float,
) -> tuple[str, str, float | None]:
    if not entry.get("hit_found", False):
        return "-", "no_hit", None

    pident = entry.get("pident")
    evalue = entry.get("evalue")
    length = entry.get("length")

    if pident is None or evalue is None or length is None:
        return "x", "missing_metrics", None

    try:
        pident = float(pident)
        evalue = float(evalue)
        length = int(length)
    except Exception:
        return "x", "bad_metrics", None

    if not (pident >= min_pident and evalue <= max_evalue and length >= min_len):
        return "~", "weak_thresholds", None

    coverage = None
    query_protein_length = entry.get("query_protein_length")
    if query_protein_length is not None:
        try:
            query_protein_length = int(query_protein_length)
            if query_protein_length > 0:
                coverage = length / query_protein_length
        except Exception:
            coverage = None

    if coverage is not None and coverage < min_coverage:
        return "~", "low_coverage", coverage

    if require_reciprocal_for_weak and pident < weak_pident_threshold:
        rec = reciprocal_lookup.get(entry.get("neighbor_name", ""))
        if rec is None:
            return "~", "no_reciprocal_info", coverage
        if not bool(rec.get("reciprocal_pass", False)):
            return "~", "reciprocal_failed", coverage

    return "✓", "significant", coverage


def collapse_sftp_consensus(df: pd.DataFrame) -> pd.DataFrame:
    sftp_cols = [c for c in df.columns if c in SFTP_GENES]
    if not sftp_cols:
        return df

    base_cols = [c for c in df.columns if c not in SFTP_GENES]

    def summarize_row(row):
        vals = [row[c] for c in sftp_cols]
        if any(v == "✓" for v in vals):
            return "✓"
        if any(v == "~" for v in vals):
            return "~"
        return "-"

    out = df.copy()
    out["SFTP_cluster"] = out.apply(summarize_row, axis=1)

    insert_at = min(df.columns.get_loc(c) for c in sftp_cols)
    ordered = base_cols.copy()
    ordered.insert(insert_at, "SFTP_cluster")
    return out[ordered]


def symbol_to_numeric(symbol: str) -> int:
    if symbol == "✓":
        return 2
    if symbol == "~":
        return 1
    return 0


def consensus_symbol_from_fractions(strong_fraction: float, any_fraction: float) -> str:
    if strong_fraction >= 0.60:
        return "✓"
    if any_fraction >= 0.20:
        return "~"
    return "-"


def build_pairwise_tables(
    mapping_files: list[Path],
    min_pident: float,
    max_evalue: float,
    min_len: int,
    min_coverage: float,
    require_reciprocal_for_weak: bool,
    weak_pident_threshold: float,
    exclude_loc: bool,
    exclude_genes: set[str],
):
    long_rows = []
    pairwise_summary_rows = []
    all_genes = set()
    all_targets = set()
    all_refs = set()

    for mapping_file in mapping_files:
        with mapping_file.open("r", encoding="utf-8") as f:
            data = json.load(f)

        ref = data.get("ref_assembly", "")
        tgt = data.get("target_assembly", "")
        reciprocal_lookup = load_reciprocal_lookup(mapping_file)

        all_refs.add(ref)
        all_targets.add(tgt)

        hits_found = 0
        hits_passing = 0
        total_neighbors = 0

        for entry in data.get("hits", []):
            gene = entry.get("neighbor_name")
            if should_exclude_gene(gene, exclude_loc, exclude_genes):
                continue

            all_genes.add(gene)
            total_neighbors += 1

            symbol, reason, coverage = classify_significant(
                entry=entry,
                reciprocal_lookup=reciprocal_lookup,
                min_pident=min_pident,
                max_evalue=max_evalue,
                min_len=min_len,
                min_coverage=min_coverage,
                require_reciprocal_for_weak=require_reciprocal_for_weak,
                weak_pident_threshold=weak_pident_threshold,
            )

            if entry.get("hit_found"):
                hits_found += 1
            if symbol == "✓":
                hits_passing += 1

            long_rows.append(
                {
                    "reference_species": ref,
                    "target_species": tgt,
                    "neighbor_name": gene,
                    "symbol": symbol,
                    "reason": reason,
                    "coverage": coverage,
                    "hit_found": entry.get("hit_found"),
                    "method": entry.get("method"),
                    "pident": entry.get("pident"),
                    "evalue": entry.get("evalue"),
                    "length": entry.get("length"),
                    "mapping_json": str(mapping_file),
                }
            )

        score_fraction = (hits_passing / total_neighbors) if total_neighbors else 0.0
        pairwise_summary_rows.append(
            {
                "reference_species": ref,
                "target_species": tgt,
                "hits_found": hits_found,
                "hits_passing": hits_passing,
                "total_neighbors": total_neighbors,
                "score_fraction": score_fraction,
                "mapping_json": str(mapping_file),
            }
        )

    long_df = pd.DataFrame(long_rows)
    pairwise_summary_df = pd.DataFrame(pairwise_summary_rows)
    return (
        long_df,
        pairwise_summary_df,
        sorted(all_genes),
        sorted(all_targets),
        sorted(all_refs),
    )


def build_consensus_tables(
    long_df: pd.DataFrame,
    all_targets: list[str],
    genes_to_use: list[str],
):
    consensus_rows = []
    stats_rows = []

    for tgt in all_targets:
        subset_tgt = long_df[long_df["target_species"] == tgt]
        refs_for_tgt = sorted(
            subset_tgt["reference_species"].dropna().unique().tolist()
        )
        n_refs = len(refs_for_tgt)

        row = {"target_species": tgt}
        stats_row = {"target_species": tgt, "n_references_used": n_refs}

        for gene in genes_to_use:
            sg = subset_tgt[subset_tgt["neighbor_name"] == gene]

            if sg.empty or n_refs == 0:
                strong_fraction = 0.0
                weak_fraction = 0.0
                any_fraction = 0.0
                symbol = "-"
            else:
                strong_fraction = (sg["symbol"] == "✓").sum() / n_refs
                weak_fraction = (sg["symbol"] == "~").sum() / n_refs
                any_fraction = (
                    (sg["symbol"] == "✓") | (sg["symbol"] == "~")
                ).sum() / n_refs
                symbol = consensus_symbol_from_fractions(strong_fraction, any_fraction)

            row[gene] = symbol
            stats_row[f"{gene}__strong_fraction"] = strong_fraction
            stats_row[f"{gene}__weak_fraction"] = weak_fraction
            stats_row[f"{gene}__any_fraction"] = any_fraction

        consensus_rows.append(row)
        stats_rows.append(stats_row)

    consensus_df = pd.DataFrame(consensus_rows).set_index("target_species")
    stats_df = pd.DataFrame(stats_rows).set_index("target_species")
    return consensus_df, stats_df


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


def complete_pairwise_summary(
    summary_df: pd.DataFrame,
    all_species: list[str],
) -> pd.DataFrame:
    rows = []
    existing = {
        (r["reference_species"], r["target_species"]): r
        for _, r in summary_df.iterrows()
    }

    for ref in all_species:
        for tgt in all_species:
            key = (ref, tgt)
            if key in existing:
                rec = dict(existing[key])
                rec["missing_mapping"] = False
            else:
                rec = {
                    "reference_species": ref,
                    "target_species": tgt,
                    "hits_found": None,
                    "hits_passing": None,
                    "total_neighbors": None,
                    "score_fraction": None,
                    "mapping_json": "",
                    "missing_mapping": True,
                }

            if ref == tgt and pd.isna(rec["score_fraction"]):
                rec["score_fraction"] = 1.0
                rec["missing_mapping"] = False

            rows.append(rec)

    return pd.DataFrame(rows)


def draw_pairwise_score_heatmap(
    summary_df: pd.DataFrame,
    out_png: Path,
    title: str,
    species_groups: dict[str, list[str]] | None = None,
    species_icons_dir: Path | None = None,
    sort_species: str = "groups",
):
    if summary_df.empty:
        raise ValueError("summary_df is empty, cannot draw heatmap.")

    pivot = summary_df.pivot(
        index="reference_species",
        columns="target_species",
        values="score_fraction",
    ).fillna(0.0)

    row_order = list(pivot.index)
    col_order = list(pivot.columns)

    if sort_species == "groups" and species_groups:
        ordered = order_species_by_groups(
            sorted(set(row_order) | set(col_order)), species_groups
        )
        row_order = [sp for sp in ordered if sp in pivot.index]
        col_order = [sp for sp in ordered if sp in pivot.columns]
    elif sort_species == "alphabetical":
        row_order = sorted(row_order)
        col_order = sorted(col_order)

    mat = pivot.loc[row_order, col_order].values

    n_rows = len(row_order)
    n_cols = len(col_order)
    show_groups = species_groups is not None and len(species_groups) > 0

    fig_w = max(11, 5.5 + n_cols * 0.42)
    fig_h = max(7, 2.0 + n_rows * 0.42)

    if show_groups:
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
        gs = fig.add_gridspec(1, 3, width_ratios=[1.9, 12, 0.45], wspace=0.05)
        ax_group = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_group)
        cax = fig.add_subplot(gs[0, 2])
    else:
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
        gs = fig.add_gridspec(1, 2, width_ratios=[12, 0.45], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])
        ax_group = None

    cmap = matplotlib.colormaps["Reds"].copy()
    im = ax.imshow(
        mat,
        aspect="auto",
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
    )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_order, rotation=90, fontsize=8)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_order, fontsize=8)

    ax.set_xlabel("Target species", fontsize=10)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=12, fontweight="bold")

    ax.tick_params(axis="y", pad=4)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Fraction of conserved neighbors", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    if show_groups and ax_group is not None:
        spans = get_group_spans(row_order, species_groups)

        ax_group.set_xlim(0, 1)
        ax_group.set_ylim(n_rows - 0.5, -0.5)
        ax_group.axis("off")

        for group_name, start_idx, end_idx in spans:
            color = group_color(group_name)
            y0 = start_idx - 0.5
            h = end_idx - start_idx

            # slightly more left
            ax_group.add_patch(
                Rectangle(
                    (0.45, y0),  # was 0.55
                    0.10,
                    h,
                    facecolor=color,
                    edgecolor="none",
                    clip_on=False,
                    zorder=2,
                )
            )

            ax_group.text(
                0.40,  # was 0.50
                start_idx + h / 2 - 0.5,
                group_name.replace("_", " "),
                ha="right",
                va="center",
                fontsize=7,
                fontweight="bold",
                clip_on=False,
            )

            # separator line across heatmap
            if end_idx < n_rows:
                ax.hlines(
                    end_idx - 0.5,
                    -0.5,
                    n_cols - 0.5,
                    colors=color,
                    linewidth=1.2,
                )

    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.22, top=0.90)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def draw_consensus_heatmap(
    df: pd.DataFrame,
    out_png: Path,
    title: str,
    species_groups: dict[str, list[str]] | None = None,
    species_icons_dir: Path | None = None,
    sort_species: str = "groups",
):
    if df.empty:
        raise ValueError("Consensus dataframe is empty.")

    species_order = list(df.index)
    gene_order = list(df.columns)

    if sort_species == "groups" and species_groups:
        species_order = [
            sp
            for sp in order_species_by_groups(species_order, species_groups)
            if sp in df.index
        ]
    elif sort_species == "alphabetical":
        species_order = sorted(species_order)
    else:

        def row_score(row: pd.Series) -> int:
            return sum(symbol_to_numeric(v) for v in row)

        species_order = sorted(
            species_order, key=lambda sp: row_score(df.loc[sp]), reverse=True
        )

    numeric = (
        df.loc[species_order, gene_order]
        .replace(
            {
                "-": 0,
                "~": 1,
                "✓": 2,
            }
        )
        .astype(float)
    )

    n_rows = len(species_order)
    n_cols = len(gene_order)

    show_groups = species_groups is not None and len(species_groups) > 0
    show_icons = species_icons_dir is not None

    left_extra = 0.0
    if show_icons:
        left_extra += 1.3
    if show_groups:
        left_extra += 2.2

    fig_w = min(max(10, 4.5 + left_extra + n_cols * 0.48), 20)
    fig_h = min(max(8, 2.5 + n_rows * 0.36), 18)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=180)

    cmap = ListedColormap(["#f2f2f2", "#ffe699", "#a9d18e"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    ax.imshow(
        numeric.values, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest"
    )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(gene_order, rotation=90, fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(species_order, fontsize=7)

    ax.set_title(title, fontsize=11, fontweight="bold")

    ax.set_xticks([x - 0.5 for x in range(1, n_cols)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_rows)], minor=True)
    ax.grid(which="minor", color="#999999", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    if show_groups:
        spans = get_group_spans(species_order, species_groups)
        bar_x0 = -2.2
        bar_w = 0.22
        label_x = bar_x0 - 0.08

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
                fontsize=7,
                fontweight="bold",
                clip_on=False,
            )
            if end_idx < n_rows:
                ax.hlines(
                    end_idx - 0.5, -0.5, n_cols - 0.5, colors=color, linewidth=1.8
                )

    if show_icons:
        icon_x0 = -1.15
        icon_x1 = -0.20
        for i, sp in enumerate(species_order):
            icon_path = get_species_icon_path(sp, species_icons_dir)
            if icon_path is not None:
                try:
                    img = mpimg.imread(icon_path)
                    draw_icon_preserve_aspect(
                        ax, img, icon_x0, icon_x1, i - 0.38, i + 0.38, zorder=5
                    )
                except Exception:
                    pass

    legend_x = n_cols + 0.35
    ax.add_patch(
        Rectangle(
            (legend_x, -0.15),
            0.25,
            0.25,
            facecolor="#a9d18e",
            edgecolor="gray",
            linewidth=0.5,
            clip_on=False,
        )
    )
    ax.text(
        legend_x + 0.33, -0.02, "consistent hit", va="center", fontsize=7, clip_on=False
    )

    ax.add_patch(
        Rectangle(
            (legend_x, 0.45),
            0.25,
            0.25,
            facecolor="#ffe699",
            edgecolor="gray",
            linewidth=0.5,
            clip_on=False,
        )
    )
    ax.text(
        legend_x + 0.33,
        0.58,
        "inconsistent / weak",
        va="center",
        fontsize=7,
        clip_on=False,
    )

    ax.add_patch(
        Rectangle(
            (legend_x, 1.05),
            0.25,
            0.25,
            facecolor="#f2f2f2",
            edgecolor="gray",
            linewidth=0.5,
            clip_on=False,
        )
    )
    ax.text(
        legend_x + 0.33, 1.18, "mostly absent", va="center", fontsize=7, clip_on=False
    )

    fig.subplots_adjust(left=0.22, right=0.90, bottom=0.22, top=0.90)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main():
    args = parse_args()

    run_folders = [Path(x) for x in args.run_folder]
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    exclude_genes = set(DEFAULT_EXCLUDE_GENES)
    exclude_genes.update(g.upper() for g in args.exclude_gene)

    species_groups = load_species_groups(args.species_groups_json)
    species_icons_dir = None

    mapping_files = collect_mapping_files(run_folders, args.gene)
    if not mapping_files:
        raise ValueError("No mapping.json files found in the provided run folder(s).")

    long_df, pairwise_summary_df, all_genes, all_targets, all_refs = (
        build_pairwise_tables(
            mapping_files=mapping_files,
            min_pident=args.min_pident,
            max_evalue=args.max_evalue,
            min_len=args.min_len,
            min_coverage=args.min_coverage,
            require_reciprocal_for_weak=args.require_reciprocal_for_weak,
            weak_pident_threshold=args.weak_pident_threshold,
            exclude_loc=args.exclude_loc,
            exclude_genes=exclude_genes,
        )
    )

    if long_df.empty:
        raise ValueError("No usable hits remained after filtering.")

    all_species = sorted(set(all_refs) | set(all_targets))
    pairwise_summary_df = complete_pairwise_summary(pairwise_summary_df, all_species)

    reduced_genes = [g for g in REDUCED_GENE_SET if g in all_genes]
    full_genes = [g for g in all_genes]

    full_consensus_df, full_stats_df = build_consensus_tables(
        long_df=long_df,
        all_targets=all_targets,
        genes_to_use=full_genes,
    )

    reduced_consensus_df, reduced_stats_df = build_consensus_tables(
        long_df=long_df,
        all_targets=all_targets,
        genes_to_use=reduced_genes,
    )

    if args.collapse_sftp:
        full_consensus_df = collapse_sftp_consensus(full_consensus_df)
        reduced_consensus_df = collapse_sftp_consensus(reduced_consensus_df)

    validation_df = pairwise_summary_df.copy()
    validation_df["suspicious_low_neighbor_count"] = (
        validation_df["total_neighbors"].fillna(0) < 10
    )
    validation_df["suspicious_zero_score"] = (
        validation_df["score_fraction"].fillna(-1) == 0.0
    )

    out_long = output_prefix.with_name(output_prefix.name + "_long").with_suffix(".tsv")
    out_pairwise = output_prefix.with_name(
        output_prefix.name + "_pairwise_summary"
    ).with_suffix(".tsv")
    out_validation = output_prefix.with_name(
        output_prefix.name + "_matrix_validation"
    ).with_suffix(".tsv")

    out_full_consensus = output_prefix.with_name(
        output_prefix.name + "_consensus_full"
    ).with_suffix(".tsv")
    out_full_stats = output_prefix.with_name(
        output_prefix.name + "_consensus_full_stats"
    ).with_suffix(".tsv")

    out_reduced_consensus = output_prefix.with_name(
        output_prefix.name + "_consensus_reduced"
    ).with_suffix(".tsv")
    out_reduced_stats = output_prefix.with_name(
        output_prefix.name + "_consensus_reduced_stats"
    ).with_suffix(".tsv")

    out_pairwise_heatmap = output_prefix.with_name(
        output_prefix.name + "_pairwise_heatmap"
    ).with_suffix(".png")
    out_full_heatmap = output_prefix.with_name(
        output_prefix.name + "_consensus_full"
    ).with_suffix(".png")
    out_reduced_heatmap = output_prefix.with_name(
        output_prefix.name + "_consensus_reduced"
    ).with_suffix(".png")

    long_df.to_csv(out_long, sep="\t", index=False)
    pairwise_summary_df.to_csv(out_pairwise, sep="\t", index=False)
    validation_df.to_csv(out_validation, sep="\t", index=False)

    full_consensus_df.to_csv(out_full_consensus, sep="\t")
    full_stats_df.to_csv(out_full_stats, sep="\t")

    reduced_consensus_df.to_csv(out_reduced_consensus, sep="\t")
    reduced_stats_df.to_csv(out_reduced_stats, sep="\t")

    pairwise_title = (
        args.pairwise_title or f"{args.gene} pairwise all-vs-all score heatmap"
    )
    full_title = args.full_consensus_title or f"{args.gene} consensus table (all genes)"
    reduced_title = (
        args.reduced_consensus_title or f"{args.gene} consensus table (reduced genes)"
    )

    draw_pairwise_score_heatmap(
        summary_df=pairwise_summary_df,
        out_png=out_pairwise_heatmap,
        title=pairwise_title,
        species_groups=species_groups if args.sort_species == "groups" else None,
        species_icons_dir=species_icons_dir,
        sort_species=args.sort_species,
    )

    draw_consensus_heatmap(
        df=full_consensus_df,
        out_png=out_full_heatmap,
        title=full_title,
        species_groups=species_groups if args.sort_species == "groups" else None,
        species_icons_dir=species_icons_dir,
        sort_species=args.sort_species,
    )

    draw_consensus_heatmap(
        df=reduced_consensus_df,
        out_png=out_reduced_heatmap,
        title=reduced_title,
        species_groups=species_groups if args.sort_species == "groups" else None,
        species_icons_dir=species_icons_dir,
        sort_species=args.sort_species,
    )

    print(f"Wrote long TSV: {out_long}")
    print(f"Wrote pairwise summary TSV: {out_pairwise}")
    print(f"Wrote matrix validation TSV: {out_validation}")
    print(f"Wrote full consensus TSV: {out_full_consensus}")
    print(f"Wrote full consensus stats TSV: {out_full_stats}")
    print(f"Wrote reduced consensus TSV: {out_reduced_consensus}")
    print(f"Wrote reduced consensus stats TSV: {out_reduced_stats}")
    print(f"Wrote pairwise heatmap PNG: {out_pairwise_heatmap}")
    print(f"Wrote full consensus PNG: {out_full_heatmap}")
    print(f"Wrote reduced consensus PNG: {out_reduced_heatmap}")


if __name__ == "__main__":
    main()
