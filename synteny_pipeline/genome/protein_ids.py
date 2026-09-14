from __future__ import annotations

from pathlib import Path


def _parse_attrs(attr_field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_field.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


def protein_id_for_gene_symbol_first(gff_path: Path, gene_symbol: str) -> str | None:
    # Return the first protein_id linked to the requested gene symbol
    target = gene_symbol.lower()
    transcript_ids = set()

    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature = cols[2]
            if feature not in {"mRNA", "transcript"}:
                continue

            attrs = _parse_attrs(cols[8])
            name = (
                attrs.get("gene") or attrs.get("gene_name") or attrs.get("Name") or ""
            ).lower()

            if name == target:
                transcript_id = attrs.get("ID")
                if transcript_id:
                    transcript_ids.add(transcript_id)

    if not transcript_ids:
        return None

    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue

            attrs = _parse_attrs(cols[8])
            parent = attrs.get("Parent", "")
            protein_id = attrs.get("protein_id")

            if protein_id and parent in transcript_ids:
                return protein_id

    return None
