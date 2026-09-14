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
from matplotlib.patches import Rectangle

DEFAULT_MIN_PIDENT = 30.0
DEFAULT_MAX_EVALUE = 1e-5
DEFAULT_MIN_LEN = 60
DEFAULT_MIN_COVERAGE = 0.05

DEFAULT_EXCLUDE_GENES = {"CCDST", "CRCT1"}

SFTP_GENES = {
    "RPTN",
    "FLG",
    "FLG2",
    "CRNN",
    "TCHH",
    "TCHHL1",
    "HRNR",
}

QUERIES_DIR = Path("data/queries")
DEFAULT_SPECIES_GROUPS_JSON = Path("data/assets/species_groups.json")
DEFAULT_SPECIES_ICONS_DIR = Path("data/assets")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a filtered synteny table from mapping.json files."
    )
    parser.add_argument(
        "--run-folder",
        action="append",
        required=True,
        help="Run folder to scan. Can be used multiple times.",
    )
    parser.add_argument(
        "--gene",
        required=True,
        help="Focal gene folder, e.g. HRNR or TCHH",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix, e.g. data/plots/hrnr_strict",
    )
    parser.add_argument(
        "--min-pident",
        type=float,
        default=DEFAULT_MIN_PIDENT,
        help=f"Minimum percent identity (default: {DEFAULT_MIN_PIDENT})",
    )
    parser.add_argument(
        "--max-evalue",
        type=float,
        default=DEFAULT_MAX_EVALUE,
        help=f"Maximum e-value (default: {DEFAULT_MAX_EVALUE})",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=DEFAULT_MIN_LEN,
        help=f"Minimum alignment length (default: {DEFAULT_MIN_LEN})",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE,
        help=f"Minimum coverage threshold (default: {DEFAULT_MIN_COVERAGE})",
    )
    parser.add_argument(
        "--require-reciprocal-for-weak",
        action="store_true",
        help="Require reciprocal support for weak hits if reciprocal results exist",
    )
    parser.add_argument(
        "--weak-pident-threshold",
        type=float,
        default=45.0,
        help="Hits below this pident are treated as weak for reciprocal checking",
    )
    parser.add_argument(
        "--exclude-loc",
        action="store_true",
        help="Exclude LOC* genes",
    )
    parser.add_argument(
        "--exclude-gene",
        action="append",
        default=[],
        help="Gene symbol to exclude. Can be used multiple times.",
    )
    parser.add_argument(
        "--include-failed-species",
        default="",
        help="Comma-separated species to include even without mapping.json",
    )
    parser.add_argument(
        "--sort-species",
        choices=["score", "alphabetical", "groups"],
        default="score",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom title for the PNG/table",
    )
    parser.add_argument(
        "--reference-species",
        default="homo_sapiens",
        help="Reference species row to add manually",
    )
    parser.add_argument(
        "--reference-position",
        choices=["top", "middle", "bottom"],
        default="top",
        help="Where to place the reference species",
    )
    parser.add_argument(
        "--species-labels-json",
        default=None,
        help="Optional JSON mapping internal species names to display labels",
    )
    parser.add_argument(
        "--species-groups-json",
        default=str(DEFAULT_SPECIES_GROUPS_JSON),
        help="JSON file with species groups",
    )
    parser.add_argument(
        "--species-icons-dir",
        default=str(DEFAULT_SPECIES_ICONS_DIR),
        help="Directory containing species PNGs",
    )
    parser.add_argument(
        "--collapse-sftp",
        action="store_true",
        help="Collapse SFTP family genes into one SFTP_cluster column",
    )
    parser.add_argument(
        "--exclude-species",
        action="append",
        default=[],
        help="Species to remove from final table. Can be used multiple times.",
    )
    return parser.parse_args()


def load_species_labels(path: str | None) -> dict[str, str]:
    if not path:
        return {}

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Species labels JSON not found: {path}")

    return json.loads(path.read_text(encoding="utf-8"))


def load_species_groups(path: str | None) -> dict[str, list[str]]:
    if not path:
        return {}

    path = Path(path)
    if not path.exists():
        return {}

    return json.loads(path.read_text(encoding="utf-8"))


def normalize_species_name(name: str) -> str:
    return str(name).strip()


def display_species_name(name: str, labels: dict[str, str] | None) -> str:
    if labels is None:
        return name.replace("_", " ")
    return labels.get(name, name.replace("_", " "))


def read_single_fasta_length(path: Path) -> int | None:
    if not path.exists():
        return None

    length = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            length += len(line)

    return length if length > 0 else None


def load_reciprocal_lookup(mapping_json: Path) -> dict[str, dict]:
    reciprocal_json = mapping_json.parent / "reciprocal" / "reciprocal.json"
    if not reciprocal_json.exists():
        return {}

    try:
        rows = json.loads(reciprocal_json.read_text(encoding="utf-8"))
    except Exception:
        return {}

    lookup = {}
    for row in rows:
        neighbor = row.get("neighbor")
        if neighbor:
            lookup[neighbor] = row
    return lookup


def should_skip_gene(
    gene: str | None, exclude_loc: bool, exclude_genes: set[str]
) -> bool:
    if gene is None:
        return True

    gene = gene.upper()

    if exclude_loc and gene.startswith("LOC"):
        return True
    if gene in exclude_genes:
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


def guess_query_length(entry: dict) -> int | None:
    if entry.get("method") == "blastn":
        return None
    return entry.get("query_protein_length")


def classify_hit(
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

    query_len = guess_query_length(entry)
    coverage = None

    if query_len is not None:
        try:
            query_len = int(query_len)
            if query_len > 0:
                coverage = length / query_len
        except Exception:
            coverage = None

    if coverage is not None and coverage < min_coverage:
        return "~", "low_coverage", coverage

    if require_reciprocal_for_weak and pident < weak_pident_threshold:
        reciprocal = reciprocal_lookup.get(entry.get("neighbor_name", ""))
        if reciprocal is None:
            return "~", "no_reciprocal_info", coverage
        if not bool(reciprocal.get("reciprocal_pass", False)):
            return "~", "reciprocal_failed", coverage

    return "✓", "significant", coverage


def build_filtered_rows(
    mapping_files: list[Path],
    include_failed_species: list[str],
    reference_species: str,
    exclude_loc: bool,
    exclude_genes: set[str],
    min_pident: float,
    max_evalue: float,
    min_len: int,
    min_coverage: float,
    require_reciprocal_for_weak: bool,
    weak_pident_threshold: float,
):
    species_to_hits: dict[str, dict[str, str]] = {}
    species_to_reasons: dict[str, dict[str, str]] = {}
    gene_order: list[str] = []
    detail_rows: list[dict] = []

    for mapping_file in mapping_files:
        with mapping_file.open("r", encoding="utf-8") as handle:
            data = json.load(handle)

        species = normalize_species_name(data["target_assembly"])
        reciprocal_lookup = load_reciprocal_lookup(mapping_file)

        hits = {}
        reasons = {}
        current_order = []

        for entry in data.get("hits", []):
            gene = entry.get("neighbor_name")
            if should_skip_gene(
                gene, exclude_loc=exclude_loc, exclude_genes=exclude_genes
            ):
                continue

            current_order.append(gene)

            symbol, reason, coverage = classify_hit(
                entry=entry,
                reciprocal_lookup=reciprocal_lookup,
                min_pident=min_pident,
                max_evalue=max_evalue,
                min_len=min_len,
                min_coverage=min_coverage,
                require_reciprocal_for_weak=require_reciprocal_for_weak,
                weak_pident_threshold=weak_pident_threshold,
            )

            hits[gene] = symbol
            reasons[gene] = reason

            detail_rows.append(
                {
                    "species": species,
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

        if not gene_order:
            gene_order = current_order
        else:
            for gene in current_order:
                if gene not in gene_order:
                    gene_order.append(gene)

        species_to_hits[species] = hits
        species_to_reasons[species] = reasons

    for species in include_failed_species:
        species = normalize_species_name(species)
        species_to_hits.setdefault(species, {})
        species_to_reasons.setdefault(species, {})

    reference_species = normalize_species_name(reference_species)
    if gene_order:
        species_to_hits[reference_species] = {gene: "✓" for gene in gene_order}
        species_to_reasons[reference_species] = {
            gene: "reference" for gene in gene_order
        }

    if not gene_order:
        raise ValueError("No genes found after filtering.")

    table_rows = []
    for species in species_to_hits:
        row = {"species": species}
        for gene in gene_order:
            row[gene] = species_to_hits[species].get(gene, "-")
        table_rows.append(row)

    df = pd.DataFrame(table_rows).set_index("species")

    reason_rows = []
    for species in species_to_reasons:
        row = {"species": species}
        for gene in gene_order:
            row[gene] = species_to_reasons[species].get(gene, "")
        reason_rows.append(row)

    reasons_df = pd.DataFrame(reason_rows).set_index("species")
    details_df = pd.DataFrame(detail_rows)

    return df, reasons_df, details_df


def score_row(row: pd.Series) -> int:
    return sum(1 for value in row if value == "✓")


def collapse_sftp_block(df: pd.DataFrame) -> pd.DataFrame:
    sftp_cols = [col for col in df.columns if col in SFTP_GENES]
    if not sftp_cols:
        return df

    other_cols = [col for col in df.columns if col not in SFTP_GENES]

    def summarize(row):
        values = [row[col] for col in sftp_cols]
        if any(value == "✓" for value in values):
            return "✓"
        if any(value in {"~", "x"} for value in values):
            return "~"
        return "-"

    out = df.copy()
    out["SFTP_cluster"] = out.apply(summarize, axis=1)

    first_idx = min(df.columns.get_loc(col) for col in sftp_cols)
    ordered_cols = other_cols.copy()
    ordered_cols.insert(first_idx, "SFTP_cluster")

    return out[ordered_cols]


def collapse_sftp_reasons(reasons_df: pd.DataFrame) -> pd.DataFrame:
    sftp_cols = [col for col in reasons_df.columns if col in SFTP_GENES]
    if not sftp_cols:
        return reasons_df

    other_cols = [col for col in reasons_df.columns if col not in SFTP_GENES]

    def summarize(row):
        values = [row[col] for col in sftp_cols if row[col]]
        if not values:
            return ""
        if "significant" in values or "reference" in values:
            return "significant"
        if "low_coverage" in values:
            return "low_coverage"
        if "reciprocal_failed" in values:
            return "reciprocal_failed"
        if "weak_thresholds" in values:
            return "weak_thresholds"
        return values[0]

    out = reasons_df.copy()
    out["SFTP_cluster"] = out.apply(summarize, axis=1)

    first_idx = min(reasons_df.columns.get_loc(col) for col in sftp_cols)
    ordered_cols = other_cols.copy()
    ordered_cols.insert(first_idx, "SFTP_cluster")

    return out[ordered_cols]


def sort_species(df: pd.DataFrame, mode: str) -> list[str]:
    if mode == "alphabetical":
        return sorted(df.index.tolist())
    return sorted(df.index.tolist(), key=lambda sp: score_row(df.loc[sp]), reverse=True)


def order_species_by_groups(
    df: pd.DataFrame, species_groups: dict[str, list[str]]
) -> list[str]:
    present = set(df.index.tolist())
    ordered = []
    used = set()

    for _, group_species in species_groups.items():
        for species in group_species:
            if species in present and species not in used:
                ordered.append(species)
                used.add(species)

    remaining = [species for species in df.index.tolist() if species not in used]
    remaining = sorted(remaining, key=lambda sp: score_row(df.loc[sp]), reverse=True)
    ordered.extend(remaining)
    return ordered


def count_symbols(df: pd.DataFrame) -> dict[str, int]:
    counts = {"✓": 0, "~": 0, "x": 0, "-": 0}
    for value in df.to_numpy().flatten():
        if value in counts:
            counts[value] += 1
    return counts


def get_species_icon_path(species: str, icons_dir: Path | None) -> Path | None:
    if icons_dir is None:
        return None

    candidates = [
        icons_dir / f"{species}.png",
        icons_dir / f"{species}.PNG",
        icons_dir / f"{species}.jpg",
        icons_dir / f"{species}.jpeg",
        icons_dir / f"{species}.webp",
        icons_dir / f"{species.capitalize()}.png",
        icons_dir / f"{species.capitalize()}.PNG",
    ]

    for path in candidates:
        if path.exists():
            return path

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


def get_group_spans(
    species_order: list[str], species_groups: dict[str, list[str]]
) -> list[tuple[str, int, int]]:
    spans = []
    species_pos = {species: i for i, species in enumerate(species_order)}

    for group_name, group_species in species_groups.items():
        positions = [
            species_pos[species] for species in group_species if species in species_pos
        ]
        if not positions:
            continue
        spans.append((group_name, min(positions), max(positions) + 1))

    spans.sort(key=lambda item: item[1])
    return spans


def group_color(group_name: str) -> str:
    name = group_name.lower()

    if "water" in name or "aquatic" in name:
        return "#9ecae1"
    if "amphib" in name:
        return "#6baed6"
    if "fish" in name:
        return "#3182bd"
    if "rept" in name:
        return "#a1d99b"
    if "bird" in name:
        return "#fdd0a2"
    if "marsupial" in name or "monotreme" in name:
        return "#d9b3ff"
    if "mammal" in name:
        return "#f4cccc"

    return "#d9d9d9"


def draw_png(
    df: pd.DataFrame,
    out_png: Path,
    title: str,
    species_labels: dict[str, str],
    species_icons_dir: Path | None = None,
    species_groups: dict[str, list[str]] | None = None,
):
    species_order = list(df.index)
    gene_order = list(df.columns)

    n_rows = len(species_order)
    n_cols = len(gene_order)

    show_icons = species_icons_dir is not None
    show_groups = species_groups is not None and len(species_groups) > 0

    x0 = 3.6
    if show_icons:
        x0 = 2.8
    if show_groups:
        x0 += 1.6

    y0 = 0.8
    cell_w = 1.0
    cell_h = 0.65

    fig_w = max(10, x0 + n_cols * 1.0 + 3.0)
    fig_h = max(5, 2.0 + n_rows * 0.60)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    ax.axis("off")
    ax.set_xlim(0, x0 + n_cols * cell_w + 3.0)
    ax.set_ylim(0, y0 + n_rows * cell_h + 1.8)

    ax.text(
        x0 + (n_cols * cell_w) / 2,
        y0 + n_rows * cell_h + 1.1,
        title,
        ha="center",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )

    ax.text(
        x0 - 1.4,
        y0 + n_rows * cell_h + 0.35,
        f"Species ({len(df.index)})",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
    )

    for j, gene in enumerate(gene_order):
        x = x0 + j * cell_w + cell_w / 2
        y = y0 + n_rows * cell_h + 0.35
        ax.text(x, y, gene, ha="center", va="center", fontsize=9.5, fontweight="bold")

    if show_groups:
        spans = get_group_spans(species_order, species_groups)

        bar_x0 = x0 - 2.45
        bar_x1 = x0 - 2.10
        label_x = x0 - 2.55

        for group_name, start_idx, end_idx in spans:
            color = group_color(group_name)
            top_y = y0 + (n_rows - start_idx) * cell_h
            bottom_y = y0 + (n_rows - end_idx) * cell_h
            mid_y = (top_y + bottom_y) / 2

            rect = Rectangle(
                (bar_x0, bottom_y),
                bar_x1 - bar_x0,
                top_y - bottom_y,
                facecolor=color,
                edgecolor="none",
                zorder=1,
            )
            ax.add_patch(rect)

            ax.text(
                label_x,
                mid_y,
                group_name.replace("_", " "),
                ha="right",
                va="center",
                fontsize=9,
                fontweight="bold",
            )

            if end_idx < n_rows:
                sep_y = y0 + (n_rows - end_idx) * cell_h
                ax.plot(
                    [x0 - 0.05, x0 + n_cols * cell_w],
                    [sep_y, sep_y],
                    color=color,
                    linewidth=2.0,
                    solid_capstyle="butt",
                    zorder=6,
                )

    for i, species in enumerate(species_order):
        y = y0 + (n_rows - 1 - i) * cell_h

        icon_path = get_species_icon_path(species, species_icons_dir)
        if icon_path is not None:
            try:
                img = mpimg.imread(icon_path)
                draw_icon_preserve_aspect(
                    ax=ax,
                    img=img,
                    box_x0=x0 - 1.60,
                    box_x1=x0 - 0.22,
                    box_y0=y + 0.04,
                    box_y1=y + cell_h - 0.04,
                    zorder=5,
                )
            except Exception:
                ax.text(
                    x0 - 0.12,
                    y + cell_h / 2,
                    display_species_name(species, species_labels),
                    ha="right",
                    va="center",
                    fontsize=9.5,
                    fontweight="bold",
                )
        else:
            ax.text(
                x0 - 0.12,
                y + cell_h / 2,
                display_species_name(species, species_labels),
                ha="right",
                va="center",
                fontsize=9.5,
                fontweight="bold",
            )

        for j, gene in enumerate(gene_order):
            x = x0 + j * cell_w
            symbol = df.loc[species, gene]

            if symbol == "✓":
                color = "#a9d18e"
            elif symbol == "~":
                color = "#ffe699"
            elif symbol == "x":
                color = "#f4cccc"
            else:
                color = "#f2f2f2"

            ax.add_patch(
                Rectangle(
                    (x, y),
                    cell_w,
                    cell_h,
                    facecolor=color,
                    edgecolor="#8a8a8a",
                    linewidth=0.8,
                )
            )

    ax.add_patch(
        Rectangle(
            (x0, y0),
            n_cols * cell_w,
            n_rows * cell_h,
            fill=False,
            edgecolor="black",
            linewidth=1.2,
        )
    )

    counts = count_symbols(df)

    lx = x0 + n_cols * cell_w + 0.6
    ly = y0 + n_rows * cell_h - 0.2
    ax.text(lx, ly + 0.45, "Legend", fontsize=10, fontweight="bold", ha="left")

    legend_items = [
        ("#a9d18e", f"significant hit ({counts['✓']})"),
        ("#ffe699", f"weak hit ({counts['~']})"),
        ("#f2f2f2", f"no hit ({counts['-']})"),
    ]

    for idx, (color, label) in enumerate(legend_items):
        yy = ly - idx * 0.35
        ax.add_patch(
            Rectangle(
                (lx, yy), 0.28, 0.20, facecolor=color, edgecolor="gray", linewidth=0.6
            )
        )
        ax.text(lx + 0.38, yy + 0.10, label, ha="left", va="center", fontsize=8.5)

    plt.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def ensure_included_species_present(
    ordered: list[str],
    include_failed_species: list[str],
    reference_species: str,
) -> list[str]:
    seen = set(ordered)
    for species in include_failed_species:
        if species != reference_species and species not in seen:
            ordered.append(species)
            seen.add(species)
    return ordered


def filter_out_species(df: pd.DataFrame, species_to_remove: list[str]) -> pd.DataFrame:
    if not species_to_remove:
        return df
    keep = [species for species in df.index if species not in set(species_to_remove)]
    return df.loc[keep]


def main():
    args = parse_args()

    run_folders = [Path(path) for path in args.run_folder]
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    exclude_genes = set(DEFAULT_EXCLUDE_GENES)
    exclude_genes.update(gene.upper() for gene in args.exclude_gene)

    include_failed_species = [
        normalize_species_name(species)
        for species in args.include_failed_species.split(",")
        if species.strip()
    ]
    exclude_species = [
        normalize_species_name(species) for species in args.exclude_species
    ]

    mapping_files = collect_mapping_files(run_folders, args.gene)
    species_labels = load_species_labels(args.species_labels_json)
    species_groups = load_species_groups(args.species_groups_json)
    species_icons_dir = Path(args.species_icons_dir) if args.species_icons_dir else None

    df, reasons_df, details_df = build_filtered_rows(
        mapping_files=mapping_files,
        include_failed_species=include_failed_species,
        reference_species=args.reference_species,
        exclude_loc=args.exclude_loc,
        exclude_genes=exclude_genes,
        min_pident=args.min_pident,
        max_evalue=args.max_evalue,
        min_len=args.min_len,
        min_coverage=args.min_coverage,
        require_reciprocal_for_weak=args.require_reciprocal_for_weak,
        weak_pident_threshold=args.weak_pident_threshold,
    )

    if args.sort_species == "alphabetical":
        ordered = sorted(df.index)
    elif args.sort_species == "groups" and species_groups:
        ordered = order_species_by_groups(df, species_groups)
    else:
        ordered = sort_species(df, "score")

    ref = args.reference_species
    if ref in ordered:
        ordered.remove(ref)
        if args.reference_position == "top":
            ordered = [ref] + ordered
        elif args.reference_position == "bottom":
            ordered = ordered + [ref]
        else:
            mid = len(ordered) // 2
            ordered = ordered[:mid] + [ref] + ordered[mid:]

    ordered = ensure_included_species_present(
        ordered=ordered,
        include_failed_species=include_failed_species,
        reference_species=ref,
    )

    for species in include_failed_species:
        if species not in df.index:
            df.loc[species] = ["-"] * len(df.columns)
        if species not in reasons_df.index:
            reasons_df.loc[species] = [""] * len(reasons_df.columns)

    ordered_existing = [species for species in ordered if species in df.index]
    df = df.loc[ordered_existing]
    reasons_df = reasons_df.loc[ordered_existing]

    if exclude_species:
        df = filter_out_species(df, exclude_species)
        reasons_df = filter_out_species(reasons_df, exclude_species)

    if args.collapse_sftp:
        df = collapse_sftp_block(df)
        reasons_df = collapse_sftp_reasons(reasons_df)

    out_tsv = output_prefix.with_suffix(".tsv")
    out_png = output_prefix.with_suffix(".png")
    out_reasons = output_prefix.with_name(output_prefix.name + "_reasons").with_suffix(
        ".tsv"
    )
    out_details = output_prefix.with_name(output_prefix.name + "_details").with_suffix(
        ".tsv"
    )

    title = args.title if args.title else f"{args.gene} strict synteny table"

    df.to_csv(out_tsv, sep="\t")
    reasons_df.to_csv(out_reasons, sep="\t")
    details_df.to_csv(out_details, sep="\t", index=False)

    print(f"\nStrict TSV written to: {out_tsv}")
    print(f"Strict PNG written to: {out_png}")
    print(f"Reason TSV written to: {out_reasons}")
    print(f"Details TSV written to: {out_details}\n")
    print(df)

    if species_icons_dir is not None:
        missing_icons = [
            species
            for species in df.index
            if get_species_icon_path(species, species_icons_dir) is None
        ]
        if missing_icons:
            print("\nSpecies without matching icon files:")
            for species in missing_icons:
                print(f"  {species}")

    draw_png(
        df=df,
        out_png=out_png,
        title=title,
        species_labels=species_labels,
        species_icons_dir=species_icons_dir,
        species_groups=species_groups if args.sort_species == "groups" else None,
    )


if __name__ == "__main__":
    main()
