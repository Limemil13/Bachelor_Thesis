#!/usr/bin/env python3
"""Build thesis-wide evidence tables and static overview figures.

The figures keep measurements on their own scales.  The integrated dashboard
normalizes columns only for colour; every cell is annotated with the original
value and no combined score is calculated.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[3]
GENES = ["BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN"]
COLORS = {
    "BGN": "#4C78A8",
    "DCN": "#8C564B",
    "FMOD": "#F58518",
    "PRELP": "#ECA82C",
    "EPYC": "#B279A2",
    "LUM": "#54A24B",
    "OGN": "#E45756",
    "OMD": "#9D9DA1",
}


def read_tsv(relative: str) -> pd.DataFrame:
    # Resolve all analysis inputs from the repository root.
    return pd.read_csv(ROOT / relative, sep="\t")


def save_figure(fig: plt.Figure, base: Path) -> None:
    # Save identical content as PNG for preview and SVG for scalable use.
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def percent(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    # Centralize percentage calculation so panels use the same convention.
    return 100.0 * numerator.astype(float) / denominator.astype(float)


def orthology_supported_structure_rows(
    representatives: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Separate comparative ortholog rows from annotation-conflict provenance.

    The chicken BGN-labelled transcript is structurally valid as an annotated
    locus, but its translated product is ASPN-like.  It can therefore be kept
    for provenance without being counted as a chicken BGN ortholog in the
    cross-species summary.
    """
    notes = representatives["note"].fillna("").astype(str)
    excluded_mask = notes.str.contains("locus_structure_only", regex=False)
    retained = representatives.loc[~excluded_mask].copy()
    excluded = representatives.loc[excluded_mask].copy()
    return retained, excluded


def write_orthology_supported_structure_audit(
    representatives: pd.DataFrame, excluded: pd.DataFrame
) -> None:
    """Write filtered structure tables without deleting the raw locus record."""
    tables = ROOT / "analyses/gene_structure/tables"
    extended = ROOT / "analyses/gene_structure/extended/tables"
    representatives.to_csv(
        tables / "gene_structure_representative_orthology_supported.tsv",
        sep="\t",
        index=False,
    )

    excluded_keys = set(
        zip(excluded["query_gene"], excluded["species_or_file"], strict=False)
    )

    qc = pd.read_csv(extended / "full_gene_structure_qc.tsv", sep="\t")
    qc_keys = list(zip(qc["gene"], qc["species"], strict=False))
    qc = qc.loc[[key not in excluded_keys for key in qc_keys]].copy()
    qc.to_csv(
        extended / "full_gene_structure_qc_orthology_supported.tsv",
        sep="\t",
        index=False,
    )

    blocks = pd.read_csv(
        extended / "coding_exon_splice_phase_5species.tsv", sep="\t"
    )
    block_keys = list(zip(blocks["gene"], blocks["species"], strict=False))
    blocks = blocks.loc[[key not in excluded_keys for key in block_keys]].copy()
    junctions = blocks[blocks["intron_phase_after"].notna()].copy()
    splice_rows: list[dict[str, object]] = []
    for (gene, cds_index), group in junctions.groupby(
        ["gene", "cds_index"], sort=False
    ):
        phases = sorted({int(value) for value in group["intron_phase_after"]})
        positions = group["coding_boundary_aa"].astype(float)
        normalized = group["normalized_cds_position_percent"].astype(float)
        splice_rows.append(
            {
                "gene": gene,
                "junction_after_cds": int(cds_index),
                "species_count": group["species"].nunique(),
                "intron_phases": ",".join(map(str, phases)),
                "phase_conserved_across_species": "yes" if len(phases) == 1 else "no",
                "boundary_position_min_aa": round(positions.min(), 3),
                "boundary_position_max_aa": round(positions.max(), 3),
                "boundary_position_range_aa": round(
                    positions.max() - positions.min(), 3
                ),
                "normalized_position_min_percent": round(normalized.min(), 3),
                "normalized_position_max_percent": round(normalized.max(), 3),
                "normalized_position_range_percentage_points": round(
                    normalized.max() - normalized.min(), 3
                ),
                "all_gff_phase_checks_pass": (
                    "yes" if group["phase_consistent"].eq("yes").all() else "no"
                ),
                "boundary_conserved_within_3aa": (
                    "yes" if positions.max() - positions.min() <= 3 else "no"
                ),
            }
        )
    pd.DataFrame(splice_rows).to_csv(
        extended / "splice_junction_conservation_orthology_supported.tsv",
        sep="\t",
        index=False,
    )

    isoforms = pd.read_csv(extended / "isoform_sensitivity_5species.tsv", sep="\t")
    isoform_keys = list(zip(isoforms["gene"], isoforms["species"], strict=False))
    isoforms = isoforms.loc[[key not in excluded_keys for key in isoform_keys]].copy()
    isoforms.to_csv(
        extended / "isoform_sensitivity_orthology_supported.tsv",
        sep="\t",
        index=False,
    )

    domains = pd.read_csv(extended / "representative_pfam_domain_hits.tsv", sep="\t")
    domain_keys = list(zip(domains["gene"], domains["species"], strict=False))
    domains = domains.loc[[key not in excluded_keys for key in domain_keys]].copy()
    domains.to_csv(
        extended / "representative_pfam_domain_hits_orthology_supported.tsv",
        sep="\t",
        index=False,
    )


def build_gene_structure_summary() -> pd.DataFrame:
    # Aggregate transcript-level structure and domain mappings to one row per gene.
    all_representatives = read_tsv(
        "analyses/gene_structure/tables/gene_structure_representative_5species.tsv"
    )
    representatives, excluded = orthology_supported_structure_rows(
        all_representatives
    )
    excluded_out = (
        ROOT
        / "analyses/gene_structure/tables/gene_structure_orthology_exclusions.tsv"
    )
    excluded.to_csv(excluded_out, sep="\t", index=False)
    write_orthology_supported_structure_audit(representatives, excluded)
    # Use the extended exon table because it is regenerated with the current
    # panel (including DCN).  The older compact CDS-block table predates the
    # DCN promotion and is retained only as historical output.
    blocks = read_tsv(
        "analyses/gene_structure/extended/tables/coding_exon_splice_phase_5species.tsv"
    )
    blocks["query_gene"] = blocks["gene"]
    blocks["cds_block_index"] = blocks["cds_index"]
    blocks["coding_exon_length_bp"] = blocks["length_bp"]
    excluded_keys = set(
        zip(excluded["query_gene"], excluded["species_or_file"], strict=False)
    )
    if excluded_keys:
        block_keys = list(zip(blocks["gene"], blocks["species"], strict=False))
        blocks = blocks.loc[[key not in excluded_keys for key in block_keys]].copy()
    exon_cv = (
        blocks.groupby(["query_gene", "cds_block_index"])["coding_exon_length_bp"]
        .agg(lambda x: float(x.std(ddof=0) / x.mean()) if x.mean() else np.nan)
        .groupby("query_gene")
        .median()
    )

    rows = []
    for gene, group in representatives.groupby("query_gene", sort=False):
        rows.append(
            {
                "gene": gene,
                "species_count": group["species_or_file"].nunique(),
                "cds_exon_count": int(group["cds_exon_count"].mode().iloc[0]),
                "cds_exon_count_range": (
                    f"{int(group['cds_exon_count'].min())}-"
                    f"{int(group['cds_exon_count'].max())}"
                ),
                "protein_length_min_aa": int(group["protein_length_est_aa"].min()),
                "protein_length_max_aa": int(group["protein_length_est_aa"].max()),
                "protein_length_cv_percent": round(
                    100
                    * group["protein_length_est_aa"].std(ddof=0)
                    / group["protein_length_est_aa"].mean(),
                    2,
                ),
                "gene_span_min_bp": int(group["gene_span_bp"].min()),
                "gene_span_max_bp": int(group["gene_span_bp"].max()),
                "gene_span_fold_range": round(
                    group["gene_span_bp"].max() / group["gene_span_bp"].min(), 2
                ),
                "median_cds_intron_min_bp": int(
                    group["median_cds_intron_length_bp"].min()
                ),
                "median_cds_intron_max_bp": int(
                    group["median_cds_intron_length_bp"].max()
                ),
                "median_coding_exon_length_cv": round(float(exon_cv.loc[gene]), 3),
                "structural_outlier_count": int(
                    group["structural_flag_clean"].fillna("").astype(str).ne("").sum()
                ),
            }
        )
    result = pd.DataFrame(rows)
    out = ROOT / "analyses/gene_structure/tables/gene_structure_extended_summary.tsv"
    result.to_csv(out, sep="\t", index=False)
    return result


def plot_gene_structure() -> None:
    # Compare gene spans, coding-exon counts and protein lengths across species.
    all_reps = read_tsv(
        "analyses/gene_structure/tables/gene_structure_representative_5species.tsv"
    )
    reps, excluded = orthology_supported_structure_rows(all_reps)
    gene_order = [g for g in GENES + ["OMD"] if g in set(reps["query_gene"])]
    species_order = ["human", "mouse", "cow", "chicken", "zebrafish"]

    fig, axes = plt.subplots(
        1, 3, figsize=(16, 6), gridspec_kw={"width_ratios": [1.1, 1, 1]}
    )
    exon = reps.pivot(
        index="query_gene", columns="species_or_file", values="cds_exon_count"
    ).reindex(index=gene_order, columns=species_order)
    exon_values = exon.to_numpy(dtype=float)
    exon_cmap = LinearSegmentedColormap.from_list("exons", ["#F3F7FB", "#4C78A8"])
    exon_cmap.set_bad("#D9D9D9")
    axes[0].imshow(np.ma.masked_invalid(exon_values), cmap=exon_cmap, aspect="auto")
    axes[0].set_xticks(
        range(len(species_order)), species_order, rotation=35, ha="right"
    )
    axes[0].set_yticks(range(len(gene_order)), gene_order)
    for i in range(exon_values.shape[0]):
        for j in range(exon_values.shape[1]):
            label = "excluded" if np.isnan(exon_values[i, j]) else f"{exon_values[i, j]:.0f}"
            axes[0].text(j, i, label, ha="center", va="center", fontsize=8)
    axes[0].set_title("A  Coding-exon count")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")

    for gene in gene_order:
        sub = (
            reps[reps["query_gene"].eq(gene)]
            .set_index("species_or_file")
            .reindex(species_order)
        )
        axes[1].plot(
            species_order,
            sub["protein_length_est_aa"],
            marker="o",
            linewidth=2,
            label=gene,
            color=COLORS[gene],
        )
        axes[2].plot(
            species_order,
            sub["gene_span_bp"],
            marker="o",
            linewidth=2,
            label=gene,
            color=COLORS[gene],
        )
    axes[1].set_title("B  Estimated protein length")
    axes[1].set_ylabel("Amino acids")
    axes[2].set_title("C  Genomic gene span")
    axes[2].set_ylabel("Base pairs (log scale)")
    axes[2].set_yscale("log")
    for ax in axes[1:]:
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", alpha=0.25)
    axes[2].legend(ncol=2, fontsize=8, frameon=False)
    fig.suptitle(
        "Orthology-supported gene structures conserve coding-exon number while locus span varies",
        fontsize=15,
        y=1.02,
    )
    fig.text(
        0.5,
        -0.02,
        "Representative protein-coding transcripts; chicken BGN excluded because its translated product is ASPN-like",
        ha="center",
        fontsize=9,
    )
    save_figure(
        fig,
        ROOT / "analyses/gene_structure/figures/gene_structure_conservation_overview",
    )


def plot_protein_overview() -> None:
    # Summarize identity, domains, signal peptides and sequence-QC status without
    # combining them into a single confidence score.
    summary = (
        read_tsv(
            "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_gene_summary.tsv"
        )
        .rename(columns={"Gene": "gene"})
        .set_index("gene")
        .loc[GENES]
        .reset_index()
    )
    msa = (
        read_tsv("analyses/protein_analysis/tables/msa_conservation_gene_summary.tsv")
        .set_index("gene")
        .loc[GENES]
    )
    proteins = read_tsv(
        "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv"
    )
    qc = (
        proteins.groupby(["Gene", "Protein QC status"])
        .size()
        .unstack(fill_value=0)
        .reindex(GENES)
    )
    qc_statuses = ["Supported", "Watch", "Needs inspection"]
    for status in qc_statuses:
        if status not in qc:
            qc[status] = 0

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = [COLORS[g] for g in GENES]
    axes[0, 0].bar(
        GENES, summary["Mean within-gene pairwise identity percent"], color=colors
    )
    axes[0, 0].set_ylim(0, 80)
    axes[0, 0].set_ylabel("Mean pairwise identity (%)")
    axes[0, 0].set_title("A  Cross-species MSA identity")

    x = np.arange(len(GENES))
    width = 0.36
    axes[0, 1].bar(
        x - width / 2,
        msa["mean_alignment_column_occupancy_percent"],
        width,
        label="Column occupancy",
        color="#72B7B2",
    )
    axes[0, 1].bar(
        x + width / 2,
        msa["mean_modal_conservation_high_occupancy_percent"],
        width,
        label="Modal conservation",
        color="#4C78A8",
    )
    axes[0, 1].set_xticks(x, GENES)
    axes[0, 1].set_ylim(0, 100)
    axes[0, 1].set_ylabel("Percent")
    axes[0, 1].set_title("B  Alignment occupancy and conservation")
    axes[0, 1].legend(frameon=False)

    bottom = np.zeros(len(GENES))
    status_colors = {
        "Supported": "#59A14F",
        "Watch": "#F2CF5B",
        "Needs inspection": "#E15759",
    }
    for status in qc_statuses:
        if not qc[status].sum():
            continue
        values = qc[status].to_numpy()
        axes[1, 0].bar(
            GENES, values, bottom=bottom, color=status_colors[status], label=status
        )
        bottom += values
    axes[1, 0].set_ylabel("Proteins")
    axes[1, 0].set_title("C  Integrated protein-QC calls")
    axes[1, 0].set_ylim(0, max(bottom) * 1.25)
    axes[1, 0].legend(frameon=False, fontsize=9)

    domain_pct = percent(summary["SLRP-like domain count"], summary["Protein count"])
    signal_pct = percent(summary["SignalP-positive count"], summary["Protein count"])
    pending_counts = (
        proteins.assign(_pending=proteins["SignalP prediction"].eq("pending"))
        .groupby("Gene")["_pending"]
        .sum()
        .reindex(GENES, fill_value=0)
    )
    if pending_counts.sum():
        raise ValueError("The final protein figure requires complete SignalP coverage")
    small_width = 0.32
    axes[1, 1].bar(
        x - small_width / 2,
        domain_pct,
        small_width,
        label="SLRP-like Pfam",
        color="#B279A2",
    )
    axes[1, 1].bar(
        x + small_width / 2,
        signal_pct,
        small_width,
        label="SignalP positive",
        color="#54A24B",
    )
    axes[1, 1].set_xticks(x, GENES)
    axes[1, 1].set_ylim(0, 120)
    axes[1, 1].set_ylabel("Proteins positive (%)")
    axes[1, 1].set_title("D  Family architecture and secretion")
    axes[1, 1].legend(frameon=False)

    for ax in axes.flat:
        ax.grid(axis="y", alpha=0.25)
    protein_count = int(summary["Protein count"].sum())
    fig.suptitle(
        f"Protein conservation and quality-control overview ({protein_count} canonical proteins)",
        fontsize=16,
    )
    fig.tight_layout()
    save_figure(
        fig, ROOT / "analyses/protein_analysis/figures/protein_conservation_overview"
    )


def plot_expression_details() -> None:
    # Keep local TPM values and independent growth-plate ranks in separate panels.
    age = read_tsv(
        "analyses/expression/gse114919/tables/gse114919_slrp_tibia_age_contrasts.tsv"
    )
    age = age[age["gene"].isin(GENES)].copy()
    age["condition"] = age["species"].str.title() + " " + age["zone"]
    age_order = ["Mouse PZ", "Mouse HZ", "Rat PZ", "Rat HZ"]
    age_matrix = age.pivot(
        index="gene", columns="condition", values="four_minus_one_week_delta"
    ).loc[GENES, age_order]

    bgee = read_tsv(
        "analyses/expression/tables/bgee_slrp_relevant_expression_summary.tsv"
    )
    bgee = bgee[
        bgee["thesis_gene"].isin(GENES) & bgee["query_mode"].eq("anatomy_all_data")
    ].copy()
    species_order = ["mouse", "zebrafish", "chicken"]
    bgee_columns = [
        ("direct_cartilage_or_growth_plate_call_count", "Direct cartilage"),
        ("bone_or_joint_call_count", "Bone/joint"),
        ("skeletal_development_call_count", "Skeletal dev."),
    ]

    fig, axes = plt.subplots(
        1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [1.05, 1.4]}
    )
    limit = float(np.nanmax(np.abs(age_matrix.to_numpy())))
    image = axes[0].imshow(
        age_matrix.to_numpy(), cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto"
    )
    axes[0].set_xticks(range(len(age_order)), age_order, rotation=30, ha="right")
    axes[0].set_yticks(range(len(GENES)), GENES)
    for i, gene in enumerate(GENES):
        for j, condition in enumerate(age_order):
            value = age_matrix.loc[gene, condition]
            text_color = "white" if abs(value) >= 0.58 * limit else "#303030"
            axes[0].text(
                j,
                i,
                f"{value:+.2f}",
                ha="center",
                va="center",
                fontsize=9,
                color=text_color,
            )
    axes[0].set_title("A  Four-week minus one-week expression")
    colorbar = fig.colorbar(image, ax=axes[0], fraction=0.045, pad=0.03)
    colorbar.set_label("Change on authors' processed scale")

    sub_axes = axes[1].inset_axes([0, 0, 1, 1])
    axes[1].axis("off")
    combined = []
    labels = []
    boundaries = []
    cursor = 0
    for column, title in bgee_columns:
        pivot = bgee.pivot(index="thesis_gene", columns="species", values=column).loc[
            GENES, species_order
        ]
        combined.append(pivot.to_numpy(dtype=float))
        labels.extend([f"{title}\n{s.title()}" for s in species_order])
        cursor += len(species_order)
        boundaries.append(cursor)
    bgee_matrix = np.concatenate(combined, axis=1)
    bgee_max = float(np.nanmax(bgee_matrix))
    sub_axes.imshow(bgee_matrix, cmap="YlGn", aspect="auto")
    sub_axes.set_xticks(range(len(labels)), labels, rotation=40, ha="right", fontsize=8)
    sub_axes.set_yticks(range(len(GENES)), GENES)
    for i in range(bgee_matrix.shape[0]):
        for j in range(bgee_matrix.shape[1]):
            text_color = "white" if bgee_matrix[i, j] >= 0.58 * bgee_max else "#303030"
            sub_axes.text(
                j,
                i,
                f"{bgee_matrix[i, j]:.0f}",
                ha="center",
                va="center",
                fontsize=8,
                color=text_color,
            )
    for boundary in boundaries[:-1]:
        sub_axes.axvline(boundary - 0.5, color="white", linewidth=3)
    sub_axes.set_title(
        "B  Curated Bgee call counts (database coverage, not absence tests)"
    )
    fig.suptitle("Additional expression-result interpretation", fontsize=16)
    fig.text(
        0.5,
        -0.03,
        "Age contrasts are descriptive, not formal differential expression. Zero Bgee calls mean no matching database evidence, not zero biological expression.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    save_figure(
        fig, ROOT / "analyses/expression/figures/expression_age_and_bgee_overview"
    )


def plot_synteny() -> None:
    # Display SynVoy call classes and candidate counts from finalized review tables.
    evidence = read_tsv("analyses/synteny/tables/synvoy_gene_species_evidence.tsv")
    updated_fmod_calls = {
        "catshark": "NONE",
        "zebrafish": "MEDIUM",
        "whale_shark": "MEDIUM",
    }
    for species_name, confidence in updated_fmod_calls.items():
        mask = evidence["gene"].eq("FMOD") & evidence["species"].eq(species_name)
        evidence.loc[mask, "best_confidence"] = confidence
    genes = [g for g in GENES if g in set(evidence["gene"])]
    species = [
        "mouse",
        "cow",
        "dog",
        "opossum",
        "chicken",
        "anole",
        "frog",
        "coelacanth",
        "spotted_gar",
        "zebrafish",
        "catshark",
        "elephant_shark",
        "whale_shark",
        "amphioxus",
    ]
    code = {"NONE": 0, "MEDIUM": 1, "HIGH": 2}
    pivot = evidence.pivot(
        index="gene", columns="species", values="best_confidence"
    ).loc[genes, species]
    values = pivot.apply(lambda column: column.map(code)).astype(float)
    labels = pivot.replace({"NONE": "–", "MEDIUM": "M", "HIGH": "H"})

    status_code = {"rejected": 0, "ambiguous": 1, "tentative": 2, "accepted": 3}
    status_labels = {
        "rejected": "R",
        "ambiguous": "?",
        "tentative": "T",
        "accepted": "A",
    }
    status_pivot = evidence.pivot(
        index="gene", columns="species", values="review_status"
    ).loc[genes, species]
    status_values = status_pivot.apply(lambda column: column.map(status_code)).astype(
        float
    )
    status_display = status_pivot.replace(status_labels)

    fig, axes = plt.subplots(
        2, 1, figsize=(16, 9), gridspec_kw={"height_ratios": [1, 1]}
    )
    cmap = ListedColormap(["#E5E5E5", "#F2CF5B", "#59A14F"])
    heat = values.to_numpy(dtype=float)
    axes[0].imshow(heat, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    axes[0].set_xticks(range(len(species)), species, rotation=35, ha="right")
    axes[0].set_yticks(range(len(genes)), genes)
    label_values = labels.to_numpy()
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            axes[0].text(
                j, i, label_values[i, j], ha="center", va="center", fontweight="bold"
            )
    axes[0].set_title(
        "A  Best retained SynVoy gene-of-interest confidence per target genome"
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=35)

    review_cmap = ListedColormap(["#E15759", "#B7B7B7", "#F2CF5B", "#59A14F"])
    review_heat = status_values.to_numpy(dtype=float)
    axes[1].imshow(review_heat, cmap=review_cmap, vmin=0, vmax=3, aspect="auto")
    axes[1].set_xticks(range(len(species)), species, rotation=35, ha="right")
    axes[1].set_yticks(range(len(genes)), genes)
    review_labels = status_display.to_numpy()
    for i in range(review_heat.shape[0]):
        for j in range(review_heat.shape[1]):
            axes[1].text(
                j,
                i,
                review_labels[i, j],
                ha="center",
                va="center",
                fontweight="bold",
            )
    axes[1].set_title("B  Completed integrated locus-review decision")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    fig.suptitle(
        "SynVoy overview: seven completed genes across 14 non-human genomes",
        fontsize=16,
        y=1.02,
    )
    fig.text(
        0.5,
        -0.01,
        "H/M/– = automated SynVoy result; A = accepted, T = tentative, ? = ambiguous, R = rejected. Opossum uses the reviewed matched-assembly reruns; BGN remains ambiguous because its expected interval is mostly assembly gap.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    save_figure(
        fig, ROOT / "analyses/synteny/figures/synvoy_species_confidence_heatmap"
    )


def count_newick_tips(path: Path) -> int:
    # Count terminal labels in a Newick tree without reconstructing the tree.
    if not path.exists():
        return 0
    text = path.read_text(encoding="utf-8")
    return len(re.findall(r"(?<=[(,])([^(),:;]+):", text))


def parse_model(path: Path) -> str:
    # Extract the selected substitution model from an IQ-TREE report.
    if not path.exists():
        return "not run"
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"Best-fit model according to BIC:\s*(\S+)", text)
    return match.group(1) if match else "not parsed"


def build_phylogeny_inventory() -> pd.DataFrame:
    # Record the retained compact and combined trees and their basic metadata.
    review = read_tsv(
        "analyses/phylogenetics/combined_trees/seven_gene_tree/compact_tree_gene_clade_review.tsv"
    ).set_index("gene")
    rows = []
    for gene in GENES:
        small_tree = (
            ROOT / f"analyses/phylogenetics/compact_panel/trees/{gene}/{gene}.treefile"
        )
        small_report = small_tree.with_suffix(".iqtree")
        rows.append(
            {
                "gene": gene,
                "compact_tip_count": count_newick_tips(small_tree),
                "compact_model": parse_model(small_report),
                "combined_tree_gene_tip_count": int(review.loc[gene, "tip_count"]),
                "combined_tree_largest_pure_split": int(
                    review.loc[gene, "largest_pure_gene_split_tip_count"]
                ),
                "combined_tree_unrooted_monophyletic": review.loc[
                    gene, "unrooted_monophyletic"
                ],
                "combined_tree_outside_tip": review.loc[
                    gene, "tips_outside_largest_pure_gene_split"
                ],
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(
        ROOT / "analyses/overview/tables/phylogeny_result_inventory.tsv",
        sep="\t",
        index=False,
    )
    return result


def plot_phylogeny_overview(inventory: pd.DataFrame) -> None:
    # Visualize tree coverage and selected models from the inventory table.
    review = (
        read_tsv(
            "analyses/phylogenetics/combined_trees/seven_gene_tree/compact_tree_gene_clade_review.tsv"
        )
        .set_index("gene")
        .loc[GENES]
    )
    fraction = (
        100
        * review["largest_pure_gene_split_tip_count"].astype(float)
        / review["tip_count"].astype(float)
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].bar(GENES, fraction, color=[COLORS[g] for g in GENES])
    axes[0].set_ylim(0, 105)
    axes[0].set_ylabel("Tips in largest pure gene split (%)")
    combined_tip_count = int(review["tip_count"].sum())
    axes[0].set_title(
        f"A  Gene-clade separation in the valid {combined_tip_count}-tip combined tree"
    )
    for i, gene in enumerate(GENES):
        outside = review.loc[gene, "tips_outside_largest_pure_gene_split"]
        label = (
            "complete"
            if pd.isna(outside) or str(outside).strip() == ""
            else str(outside).split("|")[1]
        )
        axes[0].text(
            i,
            fraction.iloc[i] - 7,
            label,
            ha="center",
            va="top",
            fontsize=8,
            rotation=90,
        )

    axes[1].bar(
        GENES,
        inventory["compact_tip_count"],
        color=[COLORS[g] for g in GENES],
    )
    axes[1].set_ylabel("Proteins in compact tree")
    axes[1].set_title("B  Curated proteins retained in each per-gene tree")
    for i, row in inventory.iterrows():
        axes[1].text(
            i,
            row["compact_tip_count"] + 0.2,
            row["compact_model"],
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=35,
        )
    for ax in axes:
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Phylogenetic result status", fontsize=16)
    fig.text(
        0.5,
        -0.02,
        "All displayed trees are unrooted. The compact trees use the curated candidate panel; model labels are taken from the IQ-TREE reports.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    save_figure(
        fig, ROOT / "analyses/phylogenetics/figures/phylogenetic_result_overview"
    )


def plot_integrated_dashboard() -> pd.DataFrame:
    # Normalize columns only to control colour intensity; printed cell labels stay
    # on their original scales and no overall ranking is calculated.
    integrated = (
        read_tsv("analyses/overview/seven_gene_integrated_summary.tsv")
        .set_index("gene")
        .loc[GENES]
    )
    synvoy = read_tsv("analyses/synteny/tables/synvoy_gene_species_evidence.tsv")
    updated_fmod_calls = {
        "catshark": "NONE",
        "zebrafish": "MEDIUM",
        "whale_shark": "MEDIUM",
    }
    for species_name, confidence in updated_fmod_calls.items():
        mask = synvoy["gene"].eq("FMOD") & synvoy["species"].eq(species_name)
        synvoy.loc[mask, "best_confidence"] = confidence
    syn_high = (
        synvoy.assign(high=synvoy["best_confidence"].eq("HIGH"))
        .groupby("gene")["high"]
        .mean()
        * 100
    )

    matrix = pd.DataFrame(index=GENES)
    matrix["Local P0 rank"] = integrated["local_p0_rank_of_9"]
    matrix["Mouse GP rank"] = integrated["mouse_growth_plate_mean_rank"]
    matrix["Rat GP rank"] = integrated["rat_growth_plate_mean_rank"]
    matrix["MSA identity %"] = integrated["mean_pairwise_msa_identity_percent"]
    matrix["Pfam positive %"] = percent(
        integrated["slrp_like_domain_count"], integrated["canonical_protein_count"]
    )
    # Retain the tested-protein denominator so this panel remains valid if a
    # future sequence is added before its SignalP result is available.
    signalp_tested = integrated["signalp_tested_count"].astype(float)
    matrix["SignalP positive % (tested)"] = np.where(
        signalp_tested.gt(0),
        100.0 * integrated["signalp_positive_count"].astype(float) / signalp_tested,
        np.nan,
    )
    matrix["ProtSpace centroid %"] = percent(
        integrated["protspace_centroid_correct"], integrated["canonical_protein_count"]
    )
    matrix["Pure tree split %"] = percent(
        integrated["compact_tree_largest_pure_gene_split"],
        integrated["canonical_protein_count"],
    )
    matrix["SynVoy high genomes %"] = syn_high.reindex(GENES)

    display = matrix.copy()
    normalized = pd.DataFrame(index=matrix.index)
    for column in matrix:
        values = matrix[column].astype(float)
        if column.endswith("rank") or "rank" in column.lower():
            score = (values.max() - values) / (values.max() - values.min())
        else:
            score = (values - values.min()) / (values.max() - values.min())
        normalized[column] = score

    annotations = display.copy().astype(object)
    for col in annotations:
        annotations[col] = annotations[col].map(
            lambda x: "NA" if pd.isna(x) else f"{x:.1f}"
        )

    fig, ax = plt.subplots(figsize=(15, 6))
    heat = ax.imshow(
        normalized.to_numpy(dtype=float), cmap="YlGnBu", vmin=0, vmax=1, aspect="auto"
    )
    ax.set_xticks(
        range(len(normalized.columns)), normalized.columns, rotation=35, ha="right"
    )
    ax.set_yticks(range(len(normalized.index)), normalized.index)
    annotation_values = annotations.to_numpy()
    normalized_values = normalized.to_numpy(dtype=float)
    for i in range(normalized.shape[0]):
        for j in range(normalized.shape[1]):
            text_color = "white" if normalized_values[i, j] >= 0.58 else "#303030"
            ax.text(
                j,
                i,
                annotation_values[i, j],
                ha="center",
                va="center",
                fontsize=9,
                color=text_color,
            )
    colorbar = fig.colorbar(heat, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Within-column relative support (colour only)")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35)
    ax.set_title(
        "Seven-gene evidence dashboard — original values shown; no combined score"
    )
    fig.text(
        0.5,
        -0.04,
        "Lower ranks are better. SynVoy values are the fraction of genomes whose best retained call is HIGH; all calls still require prioritized review.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    save_figure(fig, ROOT / "analyses/overview/figures/seven_gene_evidence_dashboard")
    matrix.to_csv(
        ROOT / "analyses/overview/tables/seven_gene_evidence_dashboard_values.tsv",
        sep="\t",
        index_label="gene",
    )
    return matrix


def plot_workflow() -> None:
    # Draw a static overview of how independent evidence streams feed synthesis.
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    boxes = [
        (
            0.5,
            6.6,
            3.0,
            1.4,
            "1  Candidate prioritization",
            "Local P0 RNA + direct mouse/rat\ngrowth-plate data + literature",
        ),
        (
            4.3,
            6.6,
            3.0,
            1.4,
            "2  Canonical proteins",
            "Reciprocal BLAST/annotation checks\n103 proteins across seven genes",
        ),
        (
            8.1,
            6.6,
            3.0,
            1.4,
            "3  Molecular conservation",
            "Pfam/HMMER + SignalP\nMAFFT, MSA QC, identity, logos",
        ),
        (
            11.9,
            6.6,
            3.0,
            1.4,
            "4  Protein-space QC",
            "ProtT5/ProtSpace projections\nleave-one-out centroid review",
        ),
        (
            2.4,
            3.7,
            3.2,
            1.4,
            "5  Phylogenetics",
            "Compact per-gene trees +\ncombined 103-tip tree",
        ),
        (
            6.4,
            3.7,
            3.2,
            1.4,
            "6  Gene structure",
            "Representative NCBI GFF3 transcript\nCDS exons, introns, spans; five species",
        ),
        (
            10.4,
            3.7,
            3.2,
            1.4,
            "7  Synteny",
            "SynVoy genomic neighborhoods\n7 genes complete",
        ),
        (
            6.4,
            0.8,
            3.2,
            1.4,
            "8  Integrated interpretation",
            "Conserved secreted SLRP core +\nlineage-specific uncertainty/outliers",
        ),
    ]
    for x, y, w, h, title, body in boxes:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            linewidth=1.5,
            edgecolor="#2F4B7C",
            facecolor="#EEF4FB",
        )
        ax.add_patch(patch)
        ax.text(x + 0.15, y + h - 0.3, title, fontsize=11, fontweight="bold", va="top")
        ax.text(x + 0.15, y + h - 0.68, body, fontsize=9, va="top")

    arrows = [
        ((3.5, 7.3), (4.3, 7.3)),
        ((7.3, 7.3), (8.1, 7.3)),
        ((11.1, 7.3), (11.9, 7.3)),
        ((5.8, 6.6), (4.3, 5.1)),
        ((7.0, 6.6), (8.0, 5.1)),
        ((9.6, 6.6), (12.0, 5.1)),
        ((4.0, 3.7), (7.0, 2.2)),
        ((8.0, 3.7), (8.0, 2.2)),
        ((12.0, 3.7), (9.6, 2.2)),
    ]
    for start, end in arrows:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#666666"},
        )
    ax.text(
        8,
        8.55,
        "Bachelor-thesis analysis workflow and evidence integration",
        ha="center",
        fontsize=17,
        fontweight="bold",
    )
    ax.text(
        8,
        0.25,
        "Expression prioritizes biological relevance; sequence, tree, gene-structure and synteny evidence test evolutionary identity independently.",
        ha="center",
        fontsize=9,
    )
    save_figure(fig, ROOT / "analyses/overview/figures/thesis_workflow_overview")


def build_scope_inventory() -> None:
    # Record the unit and species scope of every analysis to prevent accidental
    # claims that all methods used an identical panel.
    rows = [
        {
            "analysis": "Local expression",
            "scope": "9 SLRPs; 3 P0 Prx1-lineage chondroprogenitor samples",
            "status": "complete; descriptive",
            "source_of_truth": "analyses/expression/tables/mouse_expression_descriptive_metadata_corrected.tsv",
            "main_caveat": "not anatomically isolated growth-plate zones",
        },
        {
            "analysis": "Direct growth-plate expression",
            "scope": "9 SLRPs; mouse and rat; 4 tibial age-by-zone conditions; 5 replicates each",
            "status": "complete",
            "source_of_truth": "analyses/expression/gse114919/tables/gse114919_slrp_tibia_cross_condition_summary.tsv",
            "main_caveat": "processed values interpreted within dataset/species, not pooled with TPM",
        },
        {
            "analysis": "Canonical protein/domain/MSA/SignalP",
            "scope": "103 proteins; 7 genes; 13-16 species per gene",
            "status": "complete; SignalP calls available for all 103 proteins",
            "source_of_truth": "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv",
            "main_caveat": "deep-lineage predicted proteins require manual review",
        },
        {
            "analysis": "ProtSpace",
            "scope": "103 canonical SLRPs plus 4 controls",
            "status": "complete",
            "source_of_truth": "analyses/protein_analysis/protspace/tables/protspace_canonical_embedding_qc.tsv",
            "main_caveat": "2-D projections are supportive QC, not proof of orthology",
        },
        {
            "analysis": "Compact phylogeny",
            "scope": "seven per-gene trees and one combined 103-tip tree",
            "status": "complete; OGN amphioxus removed from confident ortholog panel on 2026-08-25",
            "source_of_truth": "analyses/phylogenetics/combined_trees/seven_gene_tree/compact_tree_gene_clade_review.tsv",
            "main_caveat": "unrooted; no direction of evolution without justified outgroup",
        },
        {
            "analysis": "Gene structure",
            "scope": "39 orthology-supported representatives for 8 genes (7 panel genes plus OMD); human, mouse, cow, chicken, zebrafish",
            "status": "complete including full transcript, splice phase, isoform sensitivity, and domain-exon mapping",
            "source_of_truth": "analyses/gene_structure/extended/tables/extended_gene_structure_gene_summary.tsv",
            "main_caveat": "five reference species; chicken BGN is retained only as excluded annotation-conflict provenance; UTR and isoform results are annotation/transcript-choice dependent",
        },
        {
            "analysis": "Pairwise coding constraint",
            "scope": "8 genes (7 panel genes plus OMD); human versus mouse, cow, chicken, and zebrafish",
            "status": "complete as supplementary NG86 screen",
            "source_of_truth": "analyses/evolutionary_rates/tables/pairwise_dn_ds_human_reference.tsv",
            "main_caveat": "whole-sequence pairwise method; deep dS saturation prevents branch/site inference",
        },
        {
            "analysis": "SynVoy synteny",
            "scope": "7 completed genes x 14 target genomes = 98 gene-genome rows",
            "status": "seven automated runs and all 98 P1/P2/P3 locus reviews complete",
            "source_of_truth": "analyses/synteny/tables/synvoy_gene_species_evidence.tsv",
            "main_caveat": "matched-assembly opossum reruns were integrated after row-level review; opossum BGN remains unresolved because its expected interval is mostly ambiguous sequence",
        },
    ]
    pd.DataFrame(rows).to_csv(
        ROOT / "analyses/overview/tables/analysis_scope_inventory.tsv",
        sep="\t",
        index=False,
    )


def main() -> None:
    # Regenerate all overview tables and figures from their primary result files.
    plt.style.use("seaborn-v0_8-whitegrid")
    (ROOT / "analyses/overview/tables").mkdir(parents=True, exist_ok=True)
    build_gene_structure_summary()
    plot_gene_structure()
    plot_expression_details()
    plot_protein_overview()
    plot_synteny()
    inventory = build_phylogeny_inventory()
    plot_phylogeny_overview(inventory)
    plot_integrated_dashboard()
    plot_workflow()
    build_scope_inventory()
    print("Built thesis-wide synthesis tables and figures.")


if __name__ == "__main__":
    main()
