from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FoundGene:
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


def find_gene_by_name_fast(gff_path: Path, gene_name: str) -> FoundGene | None:
    # Return the first matching gene entry from a GFF file
    target = gene_name.strip().lower()

    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, _, feature, start, end, _, strand, _, attributes = fields
            if feature != "gene":
                continue

            attrs = _parse_gff_attributes(attributes)
            candidate = (
                attrs.get("Name")
                or attrs.get("gene")
                or attrs.get("gene_name")
                or attrs.get("locus_tag")
            )

            if candidate and candidate.lower() == target:
                return FoundGene(
                    contig=seqid,
                    start_1based=int(start),
                    end_1based=int(end),
                    strand=strand,
                    gene_id=attrs.get("ID"),
                    gene_name=candidate,
                )

    return None
