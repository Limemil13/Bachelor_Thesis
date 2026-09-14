from __future__ import annotations

from dataclasses import dataclass

from synteny_pipeline.genome.gff_region import GeneRecord
from synteny_pipeline.synteny.neighbors import NeighborBlock


@dataclass(frozen=True)
class SyntenyProfile:
    assembly: str
    focal_name: str
    upstream: list[str]
    downstream: list[str]


def gene_label(gene: GeneRecord) -> str:
    if gene.gene_name:
        return gene.gene_name.strip().lower()
    if gene.gene_id:
        return gene.gene_id.strip().lower()
    return f"{gene.contig}:{gene.start_1based}-{gene.end_1based}".lower()


def profile_from_block(assembly: str, block: NeighborBlock) -> SyntenyProfile:
    return SyntenyProfile(
        assembly=assembly,
        focal_name=gene_label(block.focal),
        upstream=[gene_label(gene) for gene in block.upstream],
        downstream=[gene_label(gene) for gene in block.downstream],
    )


@dataclass(frozen=True)
class SyntenyComparison:
    neighbor_sets: dict[str, set[str]]
    pairwise_overlap: dict[tuple[str, str], int]


def compare_profiles(profiles: list[SyntenyProfile]) -> SyntenyComparison:
    neighbor_sets = {
        profile.assembly: set(profile.upstream + profile.downstream)
        for profile in profiles
    }

    pairwise_overlap: dict[tuple[str, str], int] = {}
    assemblies = [profile.assembly for profile in profiles]

    for i in range(len(assemblies)):
        for j in range(i + 1, len(assemblies)):
            a = assemblies[i]
            b = assemblies[j]
            pairwise_overlap[(a, b)] = len(neighbor_sets[a] & neighbor_sets[b])

    return SyntenyComparison(
        neighbor_sets=neighbor_sets,
        pairwise_overlap=pairwise_overlap,
    )
