from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pyfaidx import Fasta


class ContigNotFoundError(ValueError):
    pass


class InvalidCoordinatesError(ValueError):
    pass


@dataclass(frozen=True)
class RegionSpec:
    contig: str
    start_1based: int
    end_1based: int

    def length(self) -> int:
        return self.end_1based - self.start_1based + 1


def extract_region_to_fasta(
    genome_fasta: Path,
    region: RegionSpec,
    out_fasta: Path,
    *,
    line_width: int = 60,
) -> Path:
    # Extract a 1-based inclusive region from a genome FASTA
    genome_fasta = genome_fasta.resolve()
    out_fasta = out_fasta.resolve()
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    if region.start_1based <= 0 or region.end_1based <= 0:
        raise InvalidCoordinatesError("Coordinates must be positive.")

    if region.end_1based < region.start_1based:
        raise InvalidCoordinatesError("end_1based must be >= start_1based.")

    fa = Fasta(str(genome_fasta), rebuild=False)

    if region.contig not in fa.keys():
        some_contigs = list(fa.keys())[:10]
        raise ContigNotFoundError(
            f"Contig '{region.contig}' not found in FASTA. Example contigs: {some_contigs}"
        )

    contig_len = len(fa[region.contig])
    if region.end_1based > contig_len:
        raise InvalidCoordinatesError(
            f"Region end ({region.end_1based}) exceeds contig length ({contig_len})."
        )

    start0 = region.start_1based - 1
    end0 = region.end_1based
    seq = fa[region.contig][start0:end0].seq

    header = f"{region.contig}:{region.start_1based}-{region.end_1based}"
    with out_fasta.open("w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            handle.write(seq[i : i + line_width] + "\n")

    return out_fasta
