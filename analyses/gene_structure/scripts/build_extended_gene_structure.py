#!/usr/bin/env python3
"""Build representative transcript models and splice-phase summaries.

The analysis extracts exon, CDS, and derived UTR structure from the selected
NCBI GFF3/genome pairs. It also reconstructs CDS/protein sequences, checks
splice-frame continuity, summarizes coding isoforms, and renders the associated
figures.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parents[3]
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN", "OMD")
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
GENE_COLORS = {
    "BGN": "#4C78A8",
    "DCN": "#9C6ADE",
    "FMOD": "#F58518",
    "PRELP": "#ECA82C",
    "EPYC": "#B279A2",
    "LUM": "#54A24B",
    "OGN": "#E45756",
    "OMD": "#8C8C8C",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


class IndexedFasta:
    """Small read-only FASTA accessor using a samtools .fai file."""

    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.handle = fasta.open("rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with fasta.with_suffix(fasta.suffix + ".fai").open(encoding="utf-8") as handle:
            for line in handle:
                name, length, offset, line_bases, line_width, *_ = line.rstrip().split(
                    "\t"
                )
                self.index[name] = (
                    int(length),
                    int(offset),
                    int(line_bases),
                    int(line_width),
                )

    def fetch(self, name: str, start: int, end: int) -> str:
        """Return a 1-based inclusive interval."""
        length, offset, line_bases, line_width = self.index[name]
        if start < 1 or end > length or start > end:
            raise ValueError(
                f"Invalid FASTA interval {name}:{start}-{end} (length {length})"
            )
        pos = start - 1
        remaining = end - start + 1
        output = bytearray()
        while remaining:
            line_pos = pos % line_bases
            take = min(remaining, line_bases - line_pos)
            byte_pos = offset + (pos // line_bases) * line_width + line_pos
            self.handle.seek(byte_pos)
            output.extend(self.handle.read(take))
            pos += take
            remaining -= take
        return output.decode("ascii").upper()

    def close(self) -> None:
        self.handle.close()


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def oriented_interval(
    start: int, end: int, transcript_start: int, transcript_end: int, strand: str
) -> tuple[int, int]:
    """Return 0-based half-open transcript-oriented genomic coordinates."""
    if strand == "+":
        return start - transcript_start, end - transcript_start + 1
    return transcript_end - end, transcript_end - start + 1


def selected_gff_features(
    gff: Path, transcript_ids: set[str]
) -> dict[str, list[dict[str, object]]]:
    pattern = re.compile(
        "|".join(
            re.escape(value) for value in sorted(transcript_ids, key=len, reverse=True)
        )
    )
    found: dict[str, list[dict[str, object]]] = defaultdict(list)
    with gff.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") or not pattern.search(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2].lower() not in {
                "mrna",
                "transcript",
                "exon",
                "cds",
            }:
                continue
            attrs = parse_attrs(parts[8])
            feature_id = attrs.get("ID", "")
            parents = attrs.get("Parent", "").split(",")
            matches = [tx for tx in transcript_ids if tx == feature_id or tx in parents]
            for tx in matches:
                found[tx].append(
                    {
                        "seqid": parts[0],
                        "source": parts[1],
                        "type": parts[2].lower(),
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "phase": parts[7],
                        "attrs": attrs,
                    }
                )
    return found


def subtract_intervals(
    start: int, end: int, occupied: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    cursor = start
    for left, right in sorted(occupied):
        left, right = max(left, start), min(right, end)
        if left > right:
            continue
        if cursor < left:
            result.append((cursor, left - 1))
        cursor = max(cursor, right + 1)
    if cursor <= end:
        result.append((cursor, end))
    return result


def safe_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def modal(values: list[int]) -> int | None:
    return Counter(values).most_common(1)[0][0] if values else None


def build_isoform_sensitivity(
    all_transcript_rows: list[dict[str, str]], representatives: list[dict[str, str]]
) -> list[dict[str, object]]:
    selected = {
        (row["query_gene"], row["species_or_file"]): row for row in representatives
    }
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in all_transcript_rows:
        grouped[(row["query_gene"], row["species_or_file"])].append(row)

    output: list[dict[str, object]] = []
    for key in sorted(grouped, key=lambda x: (GENES.index(x[0]), SPECIES.index(x[1]))):
        rows = grouped[key]
        valid_counts = [
            value
            for row in rows
            if (value := safe_int(row["cds_exon_count"])) and value > 0
        ]
        valid_lengths = [
            value
            for row in rows
            if (value := safe_int(row["protein_length_est_aa"])) and value > 0
        ]
        selected_row = selected[key]
        selected_count = int(selected_row["cds_exon_count"])
        selected_length = int(selected_row["protein_length_est_aa"])
        count_set = sorted(set(valid_counts))
        length_set = sorted(set(valid_lengths))
        output.append(
            {
                "gene": key[0],
                "species": key[1],
                "annotated_transcript_count": len(rows),
                "protein_coding_transcript_count": len(valid_counts),
                "cds_exon_count_set": ",".join(map(str, count_set)),
                "protein_length_range_aa": (
                    f"{min(length_set)}-{max(length_set)}" if length_set else ""
                ),
                "selected_transcript": selected_row["transcript_id"],
                "selected_cds_exon_count": selected_count,
                "selected_protein_length_aa": selected_length,
                "selected_matches_modal_cds_count": (
                    "yes" if selected_count == modal(valid_counts) else "no"
                ),
                "alternative_coding_structure": (
                    "yes" if len(count_set) > 1 or len(length_set) > 1 else "no"
                ),
                "interpretation": (
                    "coding isoform choice can change exon count and/or protein length"
                    if len(count_set) > 1 or len(length_set) > 1
                    else "annotated coding isoforms preserve exon count and protein length"
                ),
            }
        )
    return output


def plot_full_structure(features: pd.DataFrame, qc: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(17, 21), constrained_layout=True)
    axes = axes.ravel()
    for ax, gene in zip(axes, GENES, strict=False):
        gene_features = features[features["gene"].eq(gene)]
        gene_qc = qc[qc["gene"].eq(gene)].set_index("species")
        for row_index, species in enumerate(SPECIES):
            group = gene_features[gene_features["species"].eq(species)]
            span = float(gene_qc.loc[species, "transcript_span_bp"])
            y = len(SPECIES) - 1 - row_index
            ax.hlines(y, 0, 100, color="#9A9A9A", linewidth=0.8, zorder=1)
            for _, feature in group[group["feature_type"].eq("exon")].iterrows():
                x = 100 * feature["oriented_start_bp"] / span
                width = 100 * feature["length_bp"] / span
                ax.add_patch(
                    Rectangle((x, y - 0.12), width, 0.24, color="#B9B9B9", zorder=2)
                )
            for feature_type, color, height in (
                ("five_prime_UTR", "#76A5CF", 0.18),
                ("three_prime_UTR", "#E3A45B", 0.18),
                ("CDS", GENE_COLORS[gene], 0.42),
            ):
                for _, feature in group[
                    group["feature_type"].eq(feature_type)
                ].iterrows():
                    x = 100 * feature["oriented_start_bp"] / span
                    width = 100 * feature["length_bp"] / span
                    ax.add_patch(
                        Rectangle(
                            (x, y - height / 2), width, height, color=color, zorder=3
                        )
                    )
            ax.text(102, y, f"{span / 1000:.1f} kb", va="center", fontsize=8)
        cds_count = int(gene_qc["cds_block_count"].mode().iloc[0])
        exon_count = int(gene_qc["exon_count"].mode().iloc[0])
        ax.set_title(
            f"{gene}  ({exon_count} exons; {cds_count} CDS exons)",
            loc="left",
            fontweight="bold",
        )
        ax.set_yticks(range(len(SPECIES)), list(reversed(SPECIES)))
        ax.set_xlim(-1, 116)
        ax.set_ylim(-0.7, len(SPECIES) - 0.3)
        ax.set_xlabel("Position through representative transcript (%)")
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=0.15)
    axes[-1].axis("off")
    legend = [
        Patch(color="#B9B9B9", label="Exon extent"),
        Patch(color="#76A5CF", label="5′ UTR"),
        Patch(color="#E3A45B", label="3′ UTR"),
        Patch(color="#4C78A8", label="CDS (gene-specific colour)"),
    ]
    fig.legend(
        handles=legend,
        loc="lower right",
        bbox_to_anchor=(0.97, 0.04),
        frameon=False,
        ncol=2,
    )
    fig.suptitle(
        "Full representative transcript structure across five vertebrates",
        fontsize=17,
        fontweight="bold",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_splice_phase(blocks: pd.DataFrame, qc: pd.DataFrame, output: Path) -> None:
    phase_colors = {0: "#4C78A8", 1: "#F58518", 2: "#54A24B"}
    phase_markers = {0: "o", 1: "s", 2: "^"}
    fig, axes = plt.subplots(4, 2, figsize=(17, 19), constrained_layout=True)
    axes = axes.ravel()
    for ax, gene in zip(axes, GENES, strict=False):
        gene_blocks = blocks[blocks["gene"].eq(gene)]
        gene_qc = qc[qc["gene"].eq(gene)].set_index("species")
        for row_index, species in enumerate(SPECIES):
            y = len(SPECIES) - 1 - row_index
            group = gene_blocks[gene_blocks["species"].eq(species)].sort_values(
                "cds_index"
            )
            protein_length = float(gene_qc.loc[species, "translated_length_aa"])
            ax.hlines(y, 0, protein_length, color="#BEBEBE", linewidth=8, zorder=1)
            for _, block in group.iterrows():
                x0 = (block["cumulative_cds_start_bp"] - 1) / 3
                x1 = block["cumulative_cds_end_bp"] / 3
                alpha = 0.50 if int(block["cds_index"]) % 2 else 0.85
                ax.hlines(
                    y,
                    x0,
                    x1,
                    color=GENE_COLORS[gene],
                    linewidth=8,
                    alpha=alpha,
                    zorder=2,
                )
                if str(block["intron_phase_after"]) not in {"", "nan"}:
                    phase = int(float(block["intron_phase_after"]))
                    ax.scatter(
                        x1,
                        y,
                        s=38,
                        marker=phase_markers[phase],
                        color=phase_colors[phase],
                        edgecolor="white",
                        linewidth=0.5,
                        zorder=4,
                    )
        ax.set_title(f"{gene}: coding-exon boundaries", loc="left", fontweight="bold")
        ax.set_yticks(range(len(SPECIES)), list(reversed(SPECIES)))
        ax.set_xlim(-5, max(gene_qc["translated_length_aa"]) * 1.04)
        ax.set_ylim(-0.7, len(SPECIES) - 0.3)
        ax.set_xlabel("Position in translated protein (amino acids)")
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=0.15)
    axes[-1].axis("off")
    legend = [
        plt.Line2D(
            [0],
            [0],
            marker=phase_markers[p],
            linestyle="",
            color=phase_colors[p],
            label=f"Intron phase {p}",
        )
        for p in (0, 1, 2)
    ]
    fig.legend(
        handles=legend, loc="lower right", bbox_to_anchor=(0.94, 0.08), frameon=False
    )
    fig.suptitle(
        "Coding-exon boundaries and intron phases are conserved beyond exon count",
        fontsize=17,
        fontweight="bold",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--synvoy-root",
        type=Path,
        required=True,
        help="SynVoy root containing pro_panel/genomes/gff and fna",
    )
    parser.add_argument(
        "--representatives",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/tables/gene_structure_representative_5species.tsv",
    )
    parser.add_argument(
        "--all-transcripts",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/raw/main_slrp_gene_structure_5species.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "analyses/gene_structure/extended",
    )
    args = parser.parse_args()

    representatives = read_tsv(args.representatives)
    expected_keys = {(gene, species) for gene in GENES for species in SPECIES}
    observed_keys = {
        (row["query_gene"], row["species_or_file"]) for row in representatives
    }
    if observed_keys != expected_keys:
        raise ValueError(
            f"Representative set does not contain the expected 35 rows: {expected_keys - observed_keys}"
        )

    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in representatives:
        by_species[row["species_or_file"]].append(row)

    feature_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    cds_fastas: list[tuple[str, str]] = []
    protein_fastas: list[tuple[str, str]] = []

    for species in SPECIES:
        gff = args.synvoy_root / f"pro_panel/genomes/gff/{species}.gff"
        fasta_path = args.synvoy_root / f"pro_panel/genomes/fna/{species}.fna"
        transcript_ids = {row["transcript_id"] for row in by_species[species]}
        parsed = selected_gff_features(gff, transcript_ids)
        fasta = IndexedFasta(fasta_path)
        try:
            for representative in by_species[species]:
                gene = representative["query_gene"]
                tx = representative["transcript_id"]
                items = parsed.get(tx, [])
                transcript = next(
                    (item for item in items if item["type"] in {"mrna", "transcript"}),
                    None,
                )
                exons = [item for item in items if item["type"] == "exon"]
                cds = [item for item in items if item["type"] == "cds"]
                if transcript is None or not exons or not cds:
                    raise ValueError(
                        f"Incomplete selected model in {gff}: {gene} {species} {tx}"
                    )
                strand = str(transcript["strand"])
                ordered_exons = sorted(
                    exons, key=lambda x: int(x["start"]), reverse=strand == "-"
                )
                ordered_cds = sorted(
                    cds, key=lambda x: int(x["start"]), reverse=strand == "-"
                )
                tx_start, tx_end = int(transcript["start"]), int(transcript["end"])
                protein_accessions = sorted(
                    {
                        str(item["attrs"].get("protein_id", ""))
                        for item in ordered_cds
                        if item["attrs"].get("protein_id")
                    }
                )
                protein_accession = ";".join(protein_accessions)

                exon_index_by_feature: dict[int, int] = {}
                for exon_index, exon in enumerate(ordered_exons, 1):
                    exon_index_by_feature[id(exon)] = exon_index
                    oriented_start, oriented_end = oriented_interval(
                        int(exon["start"]), int(exon["end"]), tx_start, tx_end, strand
                    )
                    base = {
                        "gene": gene,
                        "species": species,
                        "transcript_id": tx,
                        "protein_accession": protein_accession,
                        "seqid": transcript["seqid"],
                        "strand": strand,
                        "feature_type": "exon",
                        "exon_index": exon_index,
                        "cds_index": "",
                        "genomic_start": exon["start"],
                        "genomic_end": exon["end"],
                        "oriented_start_bp": oriented_start,
                        "oriented_end_bp": oriented_end,
                        "length_bp": int(exon["end"]) - int(exon["start"]) + 1,
                        "gff_phase": "",
                    }
                    feature_rows.append(base)
                    overlapping = [
                        (int(block["start"]), int(block["end"]))
                        for block in ordered_cds
                        if int(block["start"]) <= int(exon["end"])
                        and int(block["end"]) >= int(exon["start"])
                    ]
                    for utr_start, utr_end in subtract_intervals(
                        int(exon["start"]), int(exon["end"]), overlapping
                    ):
                        u_start, u_end = oriented_interval(
                            utr_start, utr_end, tx_start, tx_end, strand
                        )
                        cds_oriented_starts = [
                            oriented_interval(
                                int(block["start"]),
                                int(block["end"]),
                                tx_start,
                                tx_end,
                                strand,
                            )[0]
                            for block in ordered_cds
                        ]
                        feature_type = (
                            "five_prime_UTR"
                            if u_end <= min(cds_oriented_starts)
                            else "three_prime_UTR"
                        )
                        feature_rows.append(
                            {
                                **base,
                                "feature_type": feature_type,
                                "genomic_start": utr_start,
                                "genomic_end": utr_end,
                                "oriented_start_bp": u_start,
                                "oriented_end_bp": u_end,
                                "length_bp": utr_end - utr_start + 1,
                            }
                        )

                sequences: list[str] = []
                cumulative = 0
                for cds_index, block in enumerate(ordered_cds, 1):
                    genomic_start, genomic_end = int(block["start"]), int(block["end"])
                    sequence = fasta.fetch(
                        str(block["seqid"]), genomic_start, genomic_end
                    )
                    if strand == "-":
                        sequence = revcomp(sequence)
                    sequences.append(sequence)
                    length = len(sequence)
                    overlapping_exon = next(
                        exon
                        for exon in ordered_exons
                        if int(exon["start"]) <= genomic_start
                        and int(exon["end"]) >= genomic_end
                    )
                    exon_index = exon_index_by_feature[id(overlapping_exon)]
                    oriented_start, oriented_end = oriented_interval(
                        genomic_start, genomic_end, tx_start, tx_end, strand
                    )
                    cumulative_start = cumulative + 1
                    cumulative += length
                    intron_phase = (
                        cumulative % 3 if cds_index < len(ordered_cds) else ""
                    )
                    next_gff_phase = (
                        int(ordered_cds[cds_index]["phase"])
                        if cds_index < len(ordered_cds)
                        and str(ordered_cds[cds_index]["phase"]) in {"0", "1", "2"}
                        else ""
                    )
                    expected_next_phase = (
                        (3 - int(intron_phase)) % 3 if intron_phase != "" else ""
                    )
                    next_block = (
                        ordered_cds[cds_index] if cds_index < len(ordered_cds) else None
                    )
                    intron_length = (
                        max(
                            0,
                            abs(int(next_block["start"]) - genomic_end) - 1
                            if strand == "+"
                            else abs(genomic_start - int(next_block["end"])) - 1,
                        )
                        if next_block is not None
                        else ""
                    )
                    feature_rows.append(
                        {
                            "gene": gene,
                            "species": species,
                            "transcript_id": tx,
                            "protein_accession": protein_accession,
                            "seqid": transcript["seqid"],
                            "strand": strand,
                            "feature_type": "CDS",
                            "exon_index": exon_index,
                            "cds_index": cds_index,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "oriented_start_bp": oriented_start,
                            "oriented_end_bp": oriented_end,
                            "length_bp": length,
                            "gff_phase": block["phase"],
                        }
                    )
                    block_rows.append(
                        {
                            "gene": gene,
                            "species": species,
                            "transcript_id": tx,
                            "protein_accession": protein_accession,
                            "strand": strand,
                            "cds_index": cds_index,
                            "exon_index": exon_index,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "length_bp": length,
                            "gff_phase": block["phase"],
                            "cumulative_cds_start_bp": cumulative_start,
                            "cumulative_cds_end_bp": cumulative,
                            "coding_boundary_aa": round(cumulative / 3, 3)
                            if intron_phase != ""
                            else "",
                            "normalized_cds_position_percent": (
                                round(
                                    100
                                    * cumulative
                                    / int(representative["cds_length_bp"]),
                                    3,
                                )
                                if intron_phase != ""
                                else ""
                            ),
                            "intron_phase_after": intron_phase,
                            "expected_next_gff_phase": expected_next_phase,
                            "next_gff_phase": next_gff_phase,
                            "phase_consistent": (
                                "yes"
                                if intron_phase != ""
                                and expected_next_phase == next_gff_phase
                                else "not applicable"
                                if intron_phase == ""
                                else "no"
                            ),
                            "intron_length_bp": intron_length,
                        }
                    )

                cds_sequence = "".join(sequences)
                translation = str(Seq(cds_sequence).translate())
                internal_stops = translation[:-1].count("*")
                translated = translation.removesuffix("*")
                five_utr = sum(
                    int(row["length_bp"])
                    for row in feature_rows
                    if row["gene"] == gene
                    and row["species"] == species
                    and row["feature_type"] == "five_prime_UTR"
                )
                three_utr = sum(
                    int(row["length_bp"])
                    for row in feature_rows
                    if row["gene"] == gene
                    and row["species"] == species
                    and row["feature_type"] == "three_prime_UTR"
                )
                phase_rows = [
                    row
                    for row in block_rows
                    if row["gene"] == gene
                    and row["species"] == species
                    and row["phase_consistent"] != "not applicable"
                ]
                qc_rows.append(
                    {
                        "gene": gene,
                        "species": species,
                        "transcript_id": tx,
                        "protein_accession": protein_accession,
                        "strand": strand,
                        "transcript_span_bp": tx_end - tx_start + 1,
                        "exon_count": len(ordered_exons),
                        "cds_block_count": len(ordered_cds),
                        "cds_length_bp": len(cds_sequence),
                        "five_prime_utr_length_bp": five_utr,
                        "three_prime_utr_length_bp": three_utr,
                        "cds_length_divisible_by_three": "yes"
                        if len(cds_sequence) % 3 == 0
                        else "no",
                        "translation_starts_with_methionine": "yes"
                        if translated.startswith("M")
                        else "no",
                        "terminal_stop_present": "yes"
                        if translation.endswith("*")
                        else "no",
                        "internal_stop_count": internal_stops,
                        "translated_length_aa": len(translated),
                        "all_splice_phases_consistent_with_gff": (
                            "yes"
                            if all(
                                row["phase_consistent"] == "yes" for row in phase_rows
                            )
                            else "no"
                        ),
                        "matches_previous_exon_count": (
                            "yes"
                            if len(ordered_exons) == int(representative["exon_count"])
                            else "no"
                        ),
                        "matches_previous_cds_count": (
                            "yes"
                            if len(ordered_cds) == int(representative["cds_exon_count"])
                            else "no"
                        ),
                        "matches_previous_cds_length": (
                            "yes"
                            if len(cds_sequence) == int(representative["cds_length_bp"])
                            else "no"
                        ),
                    }
                )
                identifier = f"{gene}|{species}|{protein_accession or tx}"
                cds_fastas.append((identifier, cds_sequence))
                protein_fastas.append((identifier, translated))
        finally:
            fasta.close()

    feature_fields = [
        "gene",
        "species",
        "transcript_id",
        "protein_accession",
        "seqid",
        "strand",
        "feature_type",
        "exon_index",
        "cds_index",
        "genomic_start",
        "genomic_end",
        "oriented_start_bp",
        "oriented_end_bp",
        "length_bp",
        "gff_phase",
    ]
    block_fields = [
        "gene",
        "species",
        "transcript_id",
        "protein_accession",
        "strand",
        "cds_index",
        "exon_index",
        "genomic_start",
        "genomic_end",
        "length_bp",
        "gff_phase",
        "cumulative_cds_start_bp",
        "cumulative_cds_end_bp",
        "coding_boundary_aa",
        "normalized_cds_position_percent",
        "intron_phase_after",
        "expected_next_gff_phase",
        "next_gff_phase",
        "phase_consistent",
        "intron_length_bp",
    ]
    qc_fields = [
        "gene",
        "species",
        "transcript_id",
        "protein_accession",
        "strand",
        "transcript_span_bp",
        "exon_count",
        "cds_block_count",
        "cds_length_bp",
        "five_prime_utr_length_bp",
        "three_prime_utr_length_bp",
        "cds_length_divisible_by_three",
        "translation_starts_with_methionine",
        "terminal_stop_present",
        "internal_stop_count",
        "translated_length_aa",
        "all_splice_phases_consistent_with_gff",
        "matches_previous_exon_count",
        "matches_previous_cds_count",
        "matches_previous_cds_length",
    ]
    tables = args.out_dir / "tables"
    figures = args.out_dir / "figures"
    sequences = args.out_dir / "sequences"
    write_tsv(
        tables / "full_transcript_features_5species.tsv", feature_rows, feature_fields
    )
    write_tsv(
        tables / "coding_exon_splice_phase_5species.tsv", block_rows, block_fields
    )
    write_tsv(tables / "full_gene_structure_qc.tsv", qc_rows, qc_fields)

    splice_summary_rows: list[dict[str, object]] = []
    block_df = pd.DataFrame(block_rows)
    junction_df = block_df[block_df["intron_phase_after"].ne("")].copy()
    for (gene, cds_index), group in junction_df.groupby(
        ["gene", "cds_index"], sort=False
    ):
        phases = sorted({int(value) for value in group["intron_phase_after"]})
        positions = group["coding_boundary_aa"].astype(float)
        normalized = group["normalized_cds_position_percent"].astype(float)
        splice_summary_rows.append(
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
                "all_gff_phase_checks_pass": "yes"
                if group["phase_consistent"].eq("yes").all()
                else "no",
                "boundary_conserved_within_3aa": "yes"
                if positions.max() - positions.min() <= 3
                else "no",
            }
        )
    splice_summary_fields = list(splice_summary_rows[0])
    write_tsv(
        tables / "splice_junction_conservation_summary.tsv",
        splice_summary_rows,
        splice_summary_fields,
    )

    isoform_rows = build_isoform_sensitivity(
        read_tsv(args.all_transcripts), representatives
    )
    write_tsv(
        tables / "isoform_sensitivity_5species.tsv", isoform_rows, list(isoform_rows[0])
    )

    gene_summary_rows: list[dict[str, object]] = []
    splice_df = pd.DataFrame(splice_summary_rows)
    isoform_df = pd.DataFrame(isoform_rows)
    qc_df = pd.DataFrame(qc_rows)
    for gene in GENES:
        splice_gene = splice_df[splice_df["gene"].eq(gene)]
        isoform_gene = isoform_df[isoform_df["gene"].eq(gene)]
        qc_gene = qc_df[qc_df["gene"].eq(gene)]
        gene_summary_rows.append(
            {
                "gene": gene,
                "species_count": qc_gene["species"].nunique(),
                "representative_exon_count": int(qc_gene["exon_count"].mode().iloc[0]),
                "representative_cds_exon_count": int(
                    qc_gene["cds_block_count"].mode().iloc[0]
                ),
                "splice_junction_count": len(splice_gene),
                "phase_conserved_junction_count": int(
                    splice_gene["phase_conserved_across_species"].eq("yes").sum()
                ),
                "boundary_within_3aa_junction_count": int(
                    splice_gene["boundary_conserved_within_3aa"].eq("yes").sum()
                ),
                "median_boundary_range_aa": round(
                    float(splice_gene["boundary_position_range_aa"].median()), 3
                ),
                "species_with_alternative_coding_structure": int(
                    isoform_gene["alternative_coding_structure"].eq("yes").sum()
                ),
                "five_prime_utr_length_range_bp": f"{qc_gene['five_prime_utr_length_bp'].min()}-{qc_gene['five_prime_utr_length_bp'].max()}",
                "three_prime_utr_length_range_bp": f"{qc_gene['three_prime_utr_length_bp'].min()}-{qc_gene['three_prime_utr_length_bp'].max()}",
                "complete_start_count": int(
                    qc_gene["translation_starts_with_methionine"].eq("yes").sum()
                ),
                "terminal_stop_count": int(
                    qc_gene["terminal_stop_present"].eq("yes").sum()
                ),
                "internal_stop_total": int(qc_gene["internal_stop_count"].sum()),
                "interpretation": (
                    "splice phase and approximate coding boundaries conserved across all tested representatives"
                    if splice_gene["phase_conserved_across_species"].eq("yes").all()
                    else "at least one splice phase differs across representatives"
                ),
            }
        )
    write_tsv(
        tables / "extended_gene_structure_gene_summary.tsv",
        gene_summary_rows,
        list(gene_summary_rows[0]),
    )

    sequences.mkdir(parents=True, exist_ok=True)
    with (sequences / "representative_cds_5species.fna").open(
        "w", encoding="utf-8"
    ) as handle:
        for identifier, sequence in cds_fastas:
            handle.write(f">{identifier}\n{sequence}\n")
    with (sequences / "representative_proteins_5species.faa").open(
        "w", encoding="utf-8"
    ) as handle:
        for identifier, sequence in protein_fastas:
            handle.write(f">{identifier}\n{sequence}\n")

    plot_full_structure(
        pd.DataFrame(feature_rows), qc_df, figures / "full_gene_structure_5species"
    )
    plot_splice_phase(block_df, qc_df, figures / "coding_exon_splice_phase_map")

    print(f"Wrote {len(feature_rows)} exon/CDS/UTR feature rows")
    print(
        f"Wrote {len(block_rows)} coding-exon rows and {len(splice_summary_rows)} junction summaries"
    )
    print(f"Reconstructed {len(cds_fastas)} CDS/protein pairs")


if __name__ == "__main__":
    main()
