#!/usr/bin/env python3
"""Summarize protein and synteny evidence by sampled evolutionary lineage."""

from __future__ import annotations

import csv
import shutil
from collections import Counter
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
GENES = ("BGN", "DCN", "EPYC", "FMOD", "LUM", "OGN", "PRELP")

SPECIES = (
    {
        "protein": None,
        "synteny": "amphioxus",
        "display": "Amphioxus‡",
        "scientific": "Branchiostoma spp.",
        "lineage": "Cephalochordate",
        "role": "Deep chordate uncertainty boundary",
        "finding": "No confident canonical protein; four synteny rows ambiguous and three rejected. Protein and synteny inputs represent different Branchiostoma species.",
    },
    {
        "protein": "lamprey",
        "synteny": None,
        "display": "Sea lamprey",
        "scientific": "Petromyzon marinus",
        "lineage": "Jawless vertebrate",
        "role": "Early vertebrate diversification test",
        "finding": "Only BGN-, LUM- and PRELP-like proteins were retained; LUM remains tentative and no SynVoy comparison was run.",
    },
    {
        "protein": "callorhinchus_milii",
        "synteny": "elephant_shark",
        "display": "Elephant shark",
        "scientific": "Callorhinchus milii",
        "lineage": "Cartilaginous fish",
        "role": "Early jawed-vertebrate comparison",
        "finding": "Five canonical proteins were retained. The FMOD-like locus is tentative, but its 684-aa protein is a likely compound SLRP model and is excluded from the confident panel.",
    },
    {
        "protein": "scyliorhinus_canicula",
        "synteny": "catshark",
        "display": "Small-spotted catshark",
        "scientific": "Scyliorhinus canicula",
        "lineage": "Cartilaginous fish",
        "role": "Independent shark annotation comparison",
        "finding": "All seven genes have canonical candidates, but BGN is partial and several loci remain tentative despite recognizable SLRP architecture.",
    },
    {
        "protein": "whale_shark",
        "synteny": "whale_shark",
        "display": "Whale shark",
        "scientific": "Rhincodon typus",
        "lineage": "Cartilaginous fish",
        "role": "Independent shark protein comparison",
        "finding": "All seven genes have candidates, but partial BGN and divergent EPYC models reduce confidence; three of seven synteny rows were accepted.",
    },
    {
        "protein": "lepisosteus_oculatus",
        "synteny": "spotted_gar",
        "display": "Spotted gar",
        "scientific": "Lepisosteus oculatus",
        "lineage": "Non-teleost ray-finned fish",
        "role": "Bridge before teleost-specific duplication",
        "finding": "All seven genes were recovered and five synteny rows accepted. FMOD required an alternative retained locus and OGN remains tentative because SignalP is negative.",
    },
    {
        "protein": "zebrafish",
        "synteny": "zebrafish",
        "display": "Zebrafish",
        "scientific": "Danio rerio",
        "lineage": "Teleost",
        "role": "Teleost duplication and divergence test",
        "finding": "All seven canonical genes were retained; fmoda and fmodb provide complementary co-ortholog evidence rather than a simple one-to-one FMOD assignment.",
    },
    {
        "protein": "latimeria_chalumnae",
        "synteny": "coelacanth",
        "display": "Coelacanth",
        "scientific": "Latimeria chalumnae",
        "lineage": "Lobe-finned fish",
        "role": "Bridge toward tetrapod lineages",
        "finding": "All seven genes were recovered and five synteny rows accepted, including canonical EPYC and LUM loci in the conserved class-III genomic block.",
    },
    {
        "protein": "frog",
        "synteny": "frog",
        "display": "Western clawed frog",
        "scientific": "Xenopus tropicalis",
        "lineage": "Amphibian",
        "role": "Early tetrapod comparison",
        "finding": "All seven proteins were recovered; four synteny rows were accepted and three remained tentative because protein QC was not uniformly complete.",
    },
    {
        "protein": "anolis_carolinensis",
        "synteny": "anole",
        "display": "Green anole",
        "scientific": "Anolis carolinensis",
        "lineage": "Reptile",
        "role": "Amniote reptile comparison",
        "finding": "All seven proteins were recovered with higher median identity than the fish panel; four synteny rows were accepted and three tentative.",
    },
    {
        "protein": "chicken",
        "synteny": "chicken",
        "display": "Chicken",
        "scientific": "Gallus gallus",
        "lineage": "Bird",
        "role": "Avian amniote comparison",
        "finding": "Six canonical proteins were retained. The historical BGN-locus product is ASPN-like and was rejected as BGN, demonstrating paralog/annotation conflict rather than skeletal gene loss.",
    },
    {
        "protein": "monodelphis_domestica",
        "synteny": "opossum",
        "display": "Opossum†",
        "scientific": "Monodelphis domestica",
        "lineage": "Marsupial",
        "role": "Marsupial mammal comparison",
        "finding": "Six canonical proteins were retained, but all seven synteny rows are technically uninterpretable because the paired FASTA and GFF sequence versions differ.",
    },
    {
        "protein": "mouse",
        "synteny": "mouse",
        "display": "Mouse",
        "scientific": "Mus musculus",
        "lineage": "Placental mammal",
        "role": "Experimental and expression reference",
        "finding": "All seven proteins were recovered with high identity; six of seven synteny rows were accepted. Mouse also anchors the direct growth-plate expression and phenotype evidence.",
    },
    {
        "protein": "rattus_norvegicus",
        "synteny": None,
        "display": "Rat",
        "scientific": "Rattus norvegicus",
        "lineage": "Placental mammal",
        "role": "Independent growth-plate expression reference",
        "finding": "Six canonical proteins were retained and rat contributes direct growth-plate-zone expression, but it was not included in the SynVoy target set.",
    },
    {
        "protein": "bos_taurus",
        "synteny": "cow",
        "display": "Cattle",
        "scientific": "Bos taurus",
        "lineage": "Placental mammal",
        "role": "Well-annotated mammalian comparison",
        "finding": "All seven proteins were recovered with high identity and six of seven synteny rows were accepted.",
    },
    {
        "protein": "canis_lupus_familiaris",
        "synteny": "dog",
        "display": "Dog",
        "scientific": "Canis lupus familiaris",
        "lineage": "Placental mammal",
        "role": "Well-annotated mammalian comparison",
        "finding": "All seven proteins were recovered with the highest non-human median identity; four synteny rows were accepted and three tentative after protein QC integration.",
    },
    {
        "protein": "human",
        "synteny": None,
        "display": "Human",
        "scientific": "Homo sapiens",
        "lineage": "Placental mammal",
        "role": "Protein query and home-locus reference",
        "finding": "Human provides the reference proteins and SynVoy home loci; it is an analytical anchor, not an evolutionary ancestor of the other sampled species.",
    },
)

LINEAGE_COLORS = {
    "Cephalochordate": "#8c8c8c",
    "Jawless vertebrate": "#9467bd",
    "Cartilaginous fish": "#1f77b4",
    "Non-teleost ray-finned fish": "#2ca02c",
    "Teleost": "#66a61e",
    "Lobe-finned fish": "#1b9e77",
    "Amphibian": "#e6ab02",
    "Reptile": "#d95f02",
    "Bird": "#e7298a",
    "Marsupial": "#a6761d",
    "Placental mammal": "#d62728",
}
STATUS_COLORS = {
    "accepted": "#2ca25f",
    "tentative": "#f0a202",
    "ambiguous": "#7b8da1",
    "rejected": "#c73e4d",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: float | None, digits: int = 1) -> str:
    return "not available" if value is None else f"{value:.{digits}f}"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_rows() -> list[dict[str, str]]:
    manifest = read_tsv(
        ROOT / "analyses/protein_analysis/candidates/canonical_candidate_manifest.tsv"
    )
    protein_qc = read_tsv(
        ROOT / "analyses/protein_analysis/tables/"
        "protein_conservation_domain_msa_signalp_summary.tsv"
    )
    synteny = read_tsv(
        ROOT / "analyses/synteny/tables/synvoy_gene_species_evidence.tsv"
    )

    output: list[dict[str, str]] = []
    for species in SPECIES:
        protein_name = species["protein"]
        synteny_name = species["synteny"]
        manifest_rows = [row for row in manifest if row["species"] == protein_name]
        qc_rows = [row for row in protein_qc if row["Species"] == protein_name]
        synteny_rows = [row for row in synteny if row["species"] == synteny_name]

        identities = [float(row["forward_pident"]) for row in manifest_rows]
        coverages = [float(row["forward_qcovs"]) for row in manifest_rows]
        gaps = [float(row["Gap percentage"]) for row in qc_rows]
        statuses = Counter(row["review_status"] for row in synteny_rows)
        qc_statuses = Counter(row["Protein QC status"] for row in qc_rows)
        signalp_positive = sum(row["SignalP prediction"] == "SP" for row in qc_rows)
        signalp_negative = sum(
            row["SignalP prediction"] not in {"SP", "pending"} for row in qc_rows
        )
        signalp_pending = sum(row["SignalP prediction"] == "pending" for row in qc_rows)
        genes = sorted(row["gene"] for row in manifest_rows)
        missing = sorted(set(GENES) - set(genes))
        input_qc = sorted({row["input_pair_qc"] for row in synteny_rows})

        output.append(
            {
                "display_species": species["display"],
                "scientific_name": species["scientific"],
                "lineage": species["lineage"],
                "panel_role": species["role"],
                "canonical_proteins": str(len(manifest_rows)),
                "canonical_genes": ",".join(genes) or "none",
                "missing_from_canonical_panel": ",".join(missing) or "none",
                "median_identity_to_human_percent": fmt(
                    median(identities) if identities else None
                ),
                "identity_range_percent": (
                    f"{min(identities):.1f}-{max(identities):.1f}"
                    if identities
                    else "not available"
                ),
                "median_query_coverage_percent": fmt(
                    median(coverages) if coverages else None
                ),
                "mean_msa_gap_percent": fmt(sum(gaps) / len(gaps) if gaps else None),
                "protein_qc_supported": str(qc_statuses["Supported"]),
                "protein_qc_watch": str(qc_statuses["Watch"]),
                "protein_qc_needs_signalp": str(qc_statuses["Needs SignalP rerun"]),
                "signalp_positive": str(signalp_positive),
                "signalp_negative": str(signalp_negative),
                "signalp_pending": str(signalp_pending),
                "synvoy_tested_genes": str(len(synteny_rows)),
                "synteny_accepted": str(statuses["accepted"]),
                "synteny_tentative": str(statuses["tentative"]),
                "synteny_ambiguous": str(statuses["ambiguous"]),
                "synteny_rejected": str(statuses["rejected"]),
                "synvoy_input_qc": ",".join(input_qc) if input_qc else "not run",
                "species_level_finding": species["finding"],
            }
        )
    return output


def plot(rows: list[dict[str, str]], output_stem: Path) -> None:
    labels = [row["display_species"] for row in rows]
    y = np.arange(len(rows))
    identities = [
        np.nan
        if row["median_identity_to_human_percent"] == "not available"
        else float(row["median_identity_to_human_percent"])
        for row in rows
    ]
    protein_counts = [int(row["canonical_proteins"]) for row in rows]
    colors = [LINEAGE_COLORS[row["lineage"]] for row in rows]

    fig, (ax_identity, ax_synteny) = plt.subplots(
        1,
        2,
        figsize=(14, 9.5),
        sharey=True,
        gridspec_kw={"width_ratios": (1.15, 1)},
    )
    ax_identity.barh(y, np.nan_to_num(identities), color=colors, edgecolor="white")
    for index, (identity, count) in enumerate(
        zip(identities, protein_counts, strict=True)
    ):
        if np.isnan(identity):
            ax_identity.text(
                2, index, "no confident canonical protein", va="center", fontsize=8
            )
        else:
            ax_identity.text(
                min(identity + 1.2, 96),
                index,
                f"{identity:.1f}%  (n={count})",
                va="center",
                fontsize=8,
            )
    ax_identity.set_xlim(0, 108)
    ax_identity.set_xlabel("Median amino-acid identity to human reference (%)")
    ax_identity.set_yticks(y, labels)
    ax_identity.invert_yaxis()
    ax_identity.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax_identity.set_axisbelow(True)
    ax_identity.set_title(
        "A  Curated canonical protein panel", loc="left", fontweight="bold"
    )

    left = np.zeros(len(rows))
    for status in ("accepted", "tentative", "ambiguous", "rejected"):
        values = np.array([int(row[f"synteny_{status}"]) for row in rows])
        ax_synteny.barh(
            y,
            values,
            left=left,
            label=status.capitalize(),
            color=STATUS_COLORS[status],
            edgecolor="white",
        )
        left += values
    for index, row in enumerate(rows):
        if int(row["synvoy_tested_genes"]) == 0:
            ax_synteny.text(0.15, index, "not in SynVoy panel", va="center", fontsize=8)
    ax_synteny.set_xlim(0, 7.1)
    ax_synteny.set_xlabel("Reviewed gene loci (maximum 7)")
    ax_synteny.set_title(
        "B  Completed synteny decisions", loc="left", fontweight="bold"
    )
    ax_synteny.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax_synteny.set_axisbelow(True)
    ax_synteny.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.11),
        ncol=4,
        frameon=False,
        fontsize=8,
    )

    fig.suptitle(
        "Species-lineage summary of curated SLRP evidence",
        x=0.06,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.012,
        "† Opossum synteny is technically uninterpretable because FASTA/GFF sequence versions differ.  "
        "‡ Amphioxus protein and synteny searches used different Branchiostoma species; no absence claim is made.",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.22, right=0.98, top=0.92, bottom=0.12, wspace=0.12)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(output_stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def main() -> None:
    rows = build_rows()
    table = ROOT / "analyses/overview/tables/species_lineage_evidence_summary.tsv"
    figure = ROOT / "analyses/overview/figures/species_lineage_evidence_summary"
    write_tsv(table, rows)
    plot(rows, figure)
    thesis_directory = ROOT / "thesis"
    if thesis_directory.is_dir():
        thesis_figure = (
            thesis_directory / "figures/overview/species_lineage_evidence_summary.png"
        )
        thesis_figure.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(figure.with_suffix(".png"), thesis_figure)
    print(f"Wrote {len(rows)} species rows to {table.relative_to(ROOT)}")
    print(f"Wrote {figure.relative_to(ROOT)}.png/.svg")


if __name__ == "__main__":
    main()
