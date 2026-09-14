from __future__ import annotations

from dataclasses import dataclass

from synteny_pipeline.genome.gff_region import GeneRecord


@dataclass(frozen=True)
class NeighborBlock:
    upstream: list[GeneRecord]
    focal: GeneRecord
    downstream: list[GeneRecord]


def pick_middle_gene(genes: list[GeneRecord]) -> GeneRecord:
    """Pick the middle gene in genomic order."""
    if not genes:
        raise ValueError("No genes provided.")
    return genes[len(genes) // 2]


def neighbors_around_focal(
    genes: list[GeneRecord],
    focal: GeneRecord,
    n_each_side: int,
) -> NeighborBlock:
    """Return upstream and downstream neighbors around one focal gene."""
    if focal not in genes:
        raise ValueError("Focal gene not in gene list.")

    idx = genes.index(focal)

    return NeighborBlock(
        upstream=genes[max(0, idx - n_each_side) : idx],
        focal=focal,
        downstream=genes[idx + 1 : idx + 1 + n_each_side],
    )
