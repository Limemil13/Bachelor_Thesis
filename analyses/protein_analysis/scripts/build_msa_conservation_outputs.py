"""Generate canonical per-gene MSA conservation and species-comparison outputs."""

from __future__ import annotations

import csv
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from statistics import mean, median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from canonical_dataset import parse_header, read_fasta

BASE = Path(__file__).resolve().parents[1]
ALIGNMENT_DIR = BASE / "alignments" / "canonical"
OUT_DIR = ALIGNMENT_DIR / "conservation"
SUMMARY_OUT = BASE / "tables" / "msa_conservation_gene_summary.tsv"
GENES = ("BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM")


def pairwise_identity(seq_a: str, seq_b: str) -> tuple[float, int]:
    comparable = [
        (a, b) for a, b in zip(seq_a, seq_b, strict=False) if a != "-" and b != "-"
    ]
    if not comparable:
        return 0.0, 0
    matches = sum(a == b for a, b in comparable)
    return 100.0 * matches / len(comparable), len(comparable)


def rolling_mean(values: list[float], window: int = 15) -> list[float]:
    radius = window // 2
    out: list[float] = []
    for index in range(len(values)):
        left = max(0, index - radius)
        right = min(len(values), index + radius + 1)
        out.append(sum(values[left:right]) / (right - left))
    return out


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []

    for gene in GENES:
        alignment = ALIGNMENT_DIR / f"{gene}_canonical_aligned.faa"
        records = read_fasta(alignment)
        parsed: list[tuple[str, str, str]] = []
        for header, sequence in records:
            parsed_gene, species, accession, _ = parse_header(header)
            if parsed_gene != gene:
                raise ValueError(f"{alignment} contains {parsed_gene}, expected {gene}")
            parsed.append((species, accession, sequence))

        lengths = {len(sequence) for _, _, sequence in parsed}
        if len(lengths) != 1:
            raise ValueError(f"Unequal alignment lengths in {alignment}")
        alignment_length = lengths.pop()
        n_sequences = len(parsed)

        column_rows: list[dict[str, object]] = []
        consensus_residues: list[str] = []
        occupancies: list[float] = []
        modal_conservation: list[float] = []

        for column_index in range(alignment_length):
            residues = [sequence[column_index] for _, _, sequence in parsed]
            nongap = [residue for residue in residues if residue != "-"]
            counts = Counter(nongap)
            occupancy = len(nongap) / n_sequences
            if counts:
                consensus, consensus_count = counts.most_common(1)[0]
                modal_fraction = consensus_count / len(nongap)
                entropy = -sum(
                    (count / len(nongap)) * math.log2(count / len(nongap))
                    for count in counts.values()
                )
            else:
                consensus = "-"
                modal_fraction = 0.0
                entropy = 0.0

            occupancies.append(occupancy)
            modal_conservation.append(modal_fraction)
            if occupancy >= 0.5:
                consensus_residues.append(consensus)

            column_rows.append(
                {
                    "gene": gene,
                    "alignment_position": column_index + 1,
                    "consensus_residue": consensus,
                    "occupancy_fraction": round(occupancy, 4),
                    "modal_residue_fraction_nongap": round(modal_fraction, 4),
                    "shannon_entropy_bits": round(entropy, 4),
                    "nongap_sequence_count": len(nongap),
                }
            )

        write_tsv(
            OUT_DIR / f"{gene}_column_conservation.tsv",
            column_rows,
            list(column_rows[0]),
        )
        with (OUT_DIR / f"{gene}_majority_consensus.faa").open(
            "w", encoding="utf-8"
        ) as handle:
            handle.write(
                f">{gene}_majority_consensus columns_with_occupancy_at_least_0.5\n"
            )
            consensus_sequence = "".join(consensus_residues)
            for start in range(0, len(consensus_sequence), 60):
                handle.write(consensus_sequence[start : start + 60] + "\n")

        species = [record[0] for record in parsed]
        if len(species) != len(set(species)):
            raise ValueError(f"Species are not unique within {gene}")
        matrix: dict[tuple[int, int], float] = {}
        pair_values: list[float] = []
        for left, right in combinations(range(n_sequences), 2):
            identity, _ = pairwise_identity(parsed[left][2], parsed[right][2])
            matrix[(left, right)] = identity
            matrix[(right, left)] = identity
            pair_values.append(identity)

        matrix_rows: list[dict[str, object]] = []
        for row_index, row_species in enumerate(species):
            row: dict[str, object] = {"species": row_species}
            for column_index, column_species in enumerate(species):
                row[column_species] = (
                    100.0
                    if row_index == column_index
                    else round(matrix[(row_index, column_index)], 2)
                )
            matrix_rows.append(row)
        write_tsv(
            OUT_DIR / f"{gene}_species_pairwise_identity_percent.tsv",
            matrix_rows,
            ["species", *species],
        )

        positions = list(range(1, alignment_length + 1))
        fig, ax = plt.subplots(figsize=(12, 4.5))
        ax.plot(
            positions,
            occupancies,
            color="#6b7280",
            alpha=0.35,
            linewidth=0.8,
            label="Occupancy",
        )
        ax.plot(
            positions,
            rolling_mean(occupancies),
            color="#2563eb",
            linewidth=1.5,
            label="Occupancy (15-column mean)",
        )
        ax.plot(
            positions,
            rolling_mean(modal_conservation),
            color="#dc2626",
            linewidth=1.5,
            label="Modal amino-acid conservation (15-column mean)",
        )
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Fraction")
        ax.set_title(f"{gene}: canonical protein-alignment conservation")
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"{gene}_conservation_profile.png", dpi=300)
        plt.close(fig)

        high_occupancy = [
            modal_conservation[index]
            for index, occupancy in enumerate(occupancies)
            if occupancy >= 0.8
        ]
        summary_rows.append(
            {
                "gene": gene,
                "sequence_count": n_sequences,
                "alignment_length": alignment_length,
                "consensus_length_occupancy_ge_0.5": len(consensus_residues),
                "mean_pairwise_identity_percent": round(mean(pair_values), 2),
                "median_pairwise_identity_percent": round(median(pair_values), 2),
                "minimum_pairwise_identity_percent": round(min(pair_values), 2),
                "mean_alignment_column_occupancy_percent": round(
                    100 * mean(occupancies), 2
                ),
                "mean_modal_conservation_high_occupancy_percent": round(
                    100 * mean(high_occupancy), 2
                ),
                "fully_conserved_ungapped_columns": sum(
                    occupancy == 1.0 and conservation == 1.0
                    for occupancy, conservation in zip(
                        occupancies, modal_conservation, strict=False
                    )
                ),
            }
        )

    write_tsv(SUMMARY_OUT, summary_rows, list(summary_rows[0]))
    print(f"Wrote {SUMMARY_OUT}")
    print(
        f"Wrote consensus, identity, column, and profile outputs for {len(GENES)} genes"
    )


if __name__ == "__main__":
    main()
