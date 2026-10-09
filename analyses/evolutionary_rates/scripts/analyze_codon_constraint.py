#!/usr/bin/env python3
"""Estimate pairwise coding-sequence constraint

Protein sequences are aligned with MAFFT, CDSs are back-translated to codon
alignments, and human-versus-target dN/dS is estimated with the Nei--Gojobori
(NG86) method implemented in Biopython.  This is a supplementary descriptive
constraint analysis, not a branch/site selection test.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import BiopythonExperimentalWarning
from Bio.Seq import Seq
from matplotlib.lines import Line2D

warnings.simplefilter("ignore", BiopythonExperimentalWarning)
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN", "OMD")
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
SPECIES_COLORS = {
    "mouse": "#4C78A8",
    "cow": "#F58518",
    "chicken": "#54A24B",
    "zebrafish": "#B279A2",
}
EXCLUDED_COMPARISONS = {
    ("BGN", "chicken"): (
        "excluded because the selected chicken BGN-labelled model has the known "
        "ASPN-like orthology/annotation conflict"
    ),
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(chunks)
                identifier = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if identifier is not None:
        records[identifier] = "".join(chunks)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n{sequence}\n")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def backtranslate(aligned_protein: str, cds: str, identifier: str) -> str:
    codons = [cds[index : index + 3] for index in range(0, len(cds), 3)]
    translated = str(Seq(cds).translate())
    if translated.endswith("*"):
        translated = translated[:-1]
        codons = codons[:-1]
    protein_ungapped = aligned_protein.replace("-", "")
    if translated != protein_ungapped:
        raise ValueError(f"Protein/CDS translation mismatch for {identifier}")
    output: list[str] = []
    codon_index = 0
    for residue in aligned_protein:
        if residue == "-":
            output.append("---")
        else:
            output.append(codons[codon_index])
            codon_index += 1
    return "".join(output)


def pairwise_ungapped_codons(sequence_a: str, sequence_b: str) -> tuple[str, str]:
    clean_a: list[str] = []
    clean_b: list[str] = []
    for index in range(0, len(sequence_a), 3):
        codon_a = sequence_a[index : index + 3]
        codon_b = sequence_b[index : index + 3]
        if len(codon_a) != 3 or len(codon_b) != 3:
            continue
        if set(codon_a + codon_b) - set("ACGT"):
            continue
        if str(Seq(codon_a).translate()) == "*" or str(Seq(codon_b).translate()) == "*":
            continue
        clean_a.append(codon_a)
        clean_b.append(codon_b)
    return "".join(clean_a), "".join(clean_b)


def finite(value: float) -> bool:
    return not (math.isnan(value) or math.isinf(value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cds",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/extended/sequences/representative_cds_5species.fna",
    )
    parser.add_argument(
        "--proteins",
        type=Path,
        default=ROOT
        / "analyses/gene_structure/extended/sequences/representative_proteins_5species.faa",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "analyses/evolutionary_rates",
    )
    parser.add_argument("--mafft", default="mafft")
    args = parser.parse_args()

    cds_records = read_fasta(args.cds)
    protein_records = read_fasta(args.proteins)
    expected_count = len(GENES) * len(SPECIES)
    if set(cds_records) != set(protein_records) or len(cds_records) != expected_count:
        raise ValueError(
            f"Expected matching {expected_count}-sequence CDS and protein FASTAs"
        )

    by_gene: dict[str, dict[str, str]] = defaultdict(dict)
    for identifier in cds_records:
        gene, species, _ = identifier.split("|", 2)
        by_gene[gene][species] = identifier

    align_dir = args.out_dir / "alignments"
    pair_rows: list[dict[str, object]] = []
    for gene in GENES:
        identifiers = {species: by_gene[gene][species] for species in SPECIES}
        protein_input = align_dir / f"{gene}_representative_proteins.faa"
        protein_alignment = align_dir / f"{gene}_representative_proteins_aligned.faa"
        codon_alignment_path = align_dir / f"{gene}_representative_codon_alignment.fna"
        write_fasta(
            protein_input,
            {
                identifier: protein_records[identifier]
                for identifier in identifiers.values()
            },
        )
        result = subprocess.run(
            [args.mafft, "--quiet", "--auto", str(protein_input)],
            check=True,
            capture_output=True,
            text=True,
        )
        protein_alignment.write_text(result.stdout, encoding="utf-8")
        aligned_proteins = read_fasta(protein_alignment)
        codon_alignment = {
            identifier: backtranslate(
                aligned_proteins[identifier], cds_records[identifier], identifier
            )
            for identifier in identifiers.values()
        }
        write_fasta(codon_alignment_path, codon_alignment)

        human_id = identifiers["human"]
        for species in SPECIES[1:]:
            target_id = identifiers[species]
            human_clean, target_clean = pairwise_ungapped_codons(
                codon_alignment[human_id], codon_alignment[target_id]
            )
            codon_count = len(human_clean) // 3
            aa_human = str(Seq(human_clean).translate())
            aa_target = str(Seq(target_clean).translate())
            aa_identity = (
                100
                * sum(a == b for a, b in zip(aa_human, aa_target, strict=False))
                / codon_count
            )
            nt_identity = (
                100
                * sum(a == b for a, b in zip(human_clean, target_clean, strict=False))
                / len(human_clean)
            )
            status = "interpretable"
            note = "descriptive pairwise NG86 estimate"
            exclusion_note = EXCLUDED_COMPARISONS.get((gene, species))
            if exclusion_note:
                d_n, d_s, omega = math.nan, math.nan, math.nan
                status = "excluded orthology conflict"
                note = exclusion_note
            else:
                try:
                    d_n, d_s = cal_dn_ds(
                        CodonSeq(human_clean), CodonSeq(target_clean), method="NG86"
                    )
                    d_n, d_s = float(d_n), float(d_s)
                    if not finite(d_n) or not finite(d_s) or d_s <= 0:
                        omega = math.nan
                        status = "not interpretable"
                        note = "dS was non-positive, infinite, or undefined"
                    else:
                        omega = d_n / d_s
                        if d_s >= 2:
                            status = "caution"
                            note = "high dS may indicate synonymous-site saturation"
                except (ValueError, ZeroDivisionError, OverflowError) as error:
                    d_n, d_s, omega = math.nan, math.nan, math.nan
                    status = "not interpretable"
                    note = f"NG86 failed: {type(error).__name__}"
            pair_rows.append(
                {
                    "gene": gene,
                    "reference_species": "human",
                    "target_species": species,
                    "human_accession": human_id.split("|", 2)[2],
                    "target_accession": target_id.split("|", 2)[2],
                    "aligned_ungapped_codon_count": codon_count,
                    "amino_acid_identity_percent": round(aa_identity, 3),
                    "nucleotide_identity_percent": round(nt_identity, 3),
                    "dN_NG86": "" if not finite(d_n) else round(d_n, 6),
                    "dS_NG86": "" if not finite(d_s) else round(d_s, 6),
                    "dN_dS_omega": "" if not finite(omega) else round(omega, 6),
                    "estimate_status": status,
                    "interpretation": (
                        "omega below 1; compatible with purifying constraint"
                        if finite(omega) and omega < 1
                        else "omega at or above 1; requires cautious follow-up"
                        if finite(omega)
                        else "omega unavailable"
                    ),
                    "note": note,
                }
            )

    tables = args.out_dir / "tables"
    fields = [
        "gene",
        "reference_species",
        "target_species",
        "human_accession",
        "target_accession",
        "aligned_ungapped_codon_count",
        "amino_acid_identity_percent",
        "nucleotide_identity_percent",
        "dN_NG86",
        "dS_NG86",
        "dN_dS_omega",
        "estimate_status",
        "interpretation",
        "note",
    ]
    write_tsv(tables / "pairwise_dn_ds_human_reference.tsv", pair_rows, fields)

    pair_df = pd.DataFrame(pair_rows)
    pair_df["omega_numeric"] = pd.to_numeric(pair_df["dN_dS_omega"], errors="coerce")
    summary_rows: list[dict[str, object]] = []
    for gene in GENES:
        group = pair_df[pair_df["gene"].eq(gene)]
        interpretable = group[group["omega_numeric"].notna()]
        summary_rows.append(
            {
                "gene": gene,
                "target_comparison_count": len(group),
                "interpretable_omega_count": len(interpretable),
                "caution_or_unavailable_count": int(
                    group["estimate_status"].ne("interpretable").sum()
                ),
                "median_amino_acid_identity_percent": round(
                    float(group["amino_acid_identity_percent"].median()), 3
                ),
                "median_dN_dS_omega": (
                    round(float(interpretable["omega_numeric"].median()), 6)
                    if len(interpretable)
                    else ""
                ),
                "maximum_dN_dS_omega": (
                    round(float(interpretable["omega_numeric"].max()), 6)
                    if len(interpretable)
                    else ""
                ),
                "all_interpretable_omega_below_one": (
                    "yes"
                    if len(interpretable) and interpretable["omega_numeric"].lt(1).all()
                    else "no"
                ),
                "interpretation": (
                    "all interpretable human-target comparisons are compatible with purifying constraint"
                    if len(interpretable) and interpretable["omega_numeric"].lt(1).all()
                    else "pairwise estimates are mixed or insufficient"
                ),
            }
        )
    write_tsv(
        tables / "codon_constraint_gene_summary.tsv",
        summary_rows,
        list(summary_rows[0]),
    )

    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(GENES))
    offsets = {"mouse": -0.24, "cow": -0.08, "chicken": 0.08, "zebrafish": 0.24}
    for species in SPECIES[1:]:
        species_rows = (
            pair_df[pair_df["target_species"].eq(species)]
            .set_index("gene")
            .reindex(GENES)
        )
        caution = species_rows["estimate_status"].eq("caution")
        ax.scatter(
            (x + offsets[species])[~caution.to_numpy()],
            species_rows.loc[~caution, "omega_numeric"],
            s=60,
            color=SPECIES_COLORS[species],
            label=species,
            zorder=3,
        )
        ax.scatter(
            (x + offsets[species])[caution.to_numpy()],
            species_rows.loc[caution, "omega_numeric"],
            s=72,
            facecolors="none",
            edgecolors=SPECIES_COLORS[species],
            linewidths=1.8,
            zorder=4,
        )
    ax.axhline(1, color="#666666", linestyle="--", linewidth=1.1)
    observed = pair_df["omega_numeric"].dropna()
    upper = max(1.1, float(observed.max()) * 1.12) if len(observed) else 1.1
    ax.set_ylim(0, upper)
    ax.set_xticks(x, GENES)
    ax.set_ylabel("Pairwise dN/dS (NG86)")
    ax.set_xlabel("Gene (human reference versus target species)")
    ax.set_title(
        "Interpretable pairwise dN/dS estimates are below one", fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    handles.extend(
        [
            Line2D([0], [0], color="#666666", linestyle="--", linewidth=1.1),
            Line2D(
                [0],
                [0],
                marker="o",
                color="#666666",
                markerfacecolor="none",
                linewidth=0,
                markersize=7,
            ),
        ]
    )
    labels.extend(["dN/dS = 1", "high dS caution"])
    ax.legend(handles, labels, frameon=False, ncol=3, loc="upper left")
    figures = args.out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    base = figures / "pairwise_dn_ds_human_reference"
    fig.savefig(
        base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(
        f"Wrote {len(pair_rows)} human-target NG86 comparisons across {len(GENES)} genes"
    )


if __name__ == "__main__":
    main()
