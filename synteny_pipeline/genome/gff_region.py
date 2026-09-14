from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeneRecord:
    contig: str
    start_1based: int
    end_1based: int
    strand: str
    gene_id: str | None
    gene_name: str | None


def _parse_gff_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_str.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


def genes_in_region_fast(
    gff_path: Path,
    contig: str,
    start_1based: int,
    end_1based: int,
) -> list[GeneRecord]:
    # Return gene features overlapping the requested region
    if start_1based <= 0 or end_1based <= 0 or end_1based < start_1based:
        raise ValueError(
            "Invalid coordinates: end must be >= start and both must be positive."
        )

    records: list[GeneRecord] = []

    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, _, feature, start, end, _, strand, _, attributes = fields
            if feature != "gene" or seqid != contig:
                continue

            gene_start = int(start)
            gene_end = int(end)

            if gene_end < start_1based or gene_start > end_1based:
                continue

            attrs = _parse_gff_attributes(attributes)
            records.append(
                GeneRecord(
                    contig=seqid,
                    start_1based=gene_start,
                    end_1based=gene_end,
                    strand=strand,
                    gene_id=attrs.get("ID"),
                    gene_name=(
                        attrs.get("Name")
                        or attrs.get("gene")
                        or attrs.get("gene_name")
                        or attrs.get("locus_tag")
                    ),
                )
            )

    records.sort(key=lambda gene: gene.start_1based)
    return records
