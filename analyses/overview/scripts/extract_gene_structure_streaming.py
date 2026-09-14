#!/usr/bin/env python3
"""Memory-bounded GFF gene-structure extraction for one or more focal genes.

This produces the same thesis-facing columns as SynVoy's historical extractor,
but scans a GFF in three passes and retains only matching genes, transcripts,
exons, and CDS features.  It is suitable for multi-gigabyte NCBI annotations.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

GENE_ALIASES = {
    "DCN": ("DCN", "decorin"),
    "LUM": ("LUM", "lumican"),
}
TRANSCRIPT_TYPES = {"mrna", "transcript", "rna", "primary_transcript"}


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def feature_id(attrs: dict[str, str]) -> str | None:
    return (
        attrs.get("ID")
        or attrs.get("transcript_id")
        or attrs.get("gene_id")
        or attrs.get("Name")
    )


def feature_parent(attrs: dict[str, str]) -> str | None:
    return attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("gene_id")


def matches_gene(attrs: dict[str, str], gene: str) -> bool:
    text = " ".join(attrs.values()).lower()
    return any(
        re.search(r"(^|[^a-z0-9])" + re.escape(alias.lower()) + r"([^a-z0-9]|$)", text)
        for alias in GENE_ALIASES[gene]
    )


def features(path: Path):
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, kind, start, end, _score, strand, phase, attr_text = parts
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attrs = parse_attrs(attr_text)
            yield {
                "seqid": seqid,
                "type": kind,
                "start": start_i,
                "end": end_i,
                "strand": strand,
                "phase": phase,
                "attrs": attrs,
                "id": feature_id(attrs),
                "parent": feature_parent(attrs),
            }


def infer_species(path: Path) -> str:
    return re.sub(r"\.gff3?$", "", path.name)


def summarize(path: Path, gene: str) -> list[dict[str, object]]:
    gene_hits = [
        item
        for item in features(path)
        if item["type"].lower() == "gene" and matches_gene(item["attrs"], gene)
    ]
    if not gene_hits:
        raise ValueError(f"No explicit {gene} gene feature found in {path}")

    gene_ids = {item["id"] for item in gene_hits if item["id"]}
    transcripts_by_gene: dict[str, list[dict]] = defaultdict(list)
    for item in features(path):
        parents = set((item["parent"] or "").split(","))
        if item["type"].lower() in TRANSCRIPT_TYPES:
            for gene_id in gene_ids & parents:
                transcripts_by_gene[gene_id].append(item)

    targets: list[tuple[dict, dict]] = []
    for gene_hit in gene_hits:
        gene_id = (
            gene_hit["id"]
            or gene_hit["attrs"].get("gene_id")
            or gene_hit["attrs"].get("Name")
            or "unknown_gene"
        )
        transcripts = transcripts_by_gene.get(gene_id) or [gene_hit]
        targets.extend((gene_hit, transcript) for transcript in transcripts)

    tx_ids = {transcript["id"] or gene_hit["id"] for gene_hit, transcript in targets}
    children: dict[str, list[dict]] = defaultdict(list)
    for item in features(path):
        if item["type"].lower() not in {"exon", "cds"}:
            continue
        for parent in (item["parent"] or "").split(","):
            if parent in tx_ids:
                children[parent].append(item)

    species = infer_species(path)
    rows: list[dict[str, object]] = []
    for gene_hit, transcript in targets:
        gene_id = (
            gene_hit["id"]
            or gene_hit["attrs"].get("gene_id")
            or gene_hit["attrs"].get("Name")
            or "unknown_gene"
        )
        transcript_id = transcript["id"] or gene_id
        items = children.get(transcript_id, [])
        exons = sorted(
            (item for item in items if item["type"].lower() == "exon"),
            key=lambda item: item["start"],
        )
        cds = sorted(
            (item for item in items if item["type"].lower() == "cds"),
            key=lambda item: item["start"],
        )
        cds_length = sum(abs(item["end"] - item["start"]) + 1 for item in cds)
        protein_length = cds_length // 3 if cds_length else ""
        if cds:
            cds_start = min(item["start"] for item in cds)
            cds_end = max(item["end"] for item in cds)
        else:
            cds_start, cds_end = gene_hit["start"], gene_hit["end"]
        introns = [
            max(0, right["start"] - left["end"] - 1)
            for left, right in zip(cds, cds[1:], strict=False)
        ]
        notes: list[str] = []
        if not cds:
            notes.append("no_CDS_found")
        if len(cds) <= 1:
            notes.append("single_or_no_CDS_exon")
        if protein_length and not 150 <= int(protein_length) <= 600:
            notes.append("unusual_protein_length_estimate")
        if len(transcripts_by_gene.get(gene_id, [])) > 1:
            notes.append("multiple_transcripts_for_gene")
        rows.append(
            {
                "query_gene": gene,
                "species_or_file": species,
                "matched_gene_id": gene_id,
                "matched_gene_name": gene_hit["attrs"].get("Name")
                or gene_hit["attrs"].get("gene")
                or gene_hit["attrs"].get("gene_name")
                or gene_id,
                "transcript_id": transcript_id,
                "seqid": gene_hit["seqid"],
                "strand": gene_hit["strand"],
                "gene_start": gene_hit["start"],
                "gene_end": gene_hit["end"],
                "gene_span_bp": abs(gene_hit["end"] - gene_hit["start"]) + 1,
                "cds_start": cds_start,
                "cds_end": cds_end,
                "cds_span_bp": abs(cds_end - cds_start) + 1,
                "exon_count": len(exons),
                "cds_exon_count": len(cds),
                "cds_length_bp": cds_length,
                "protein_length_est_aa": protein_length,
                "intron_count_between_cds": max(0, len(cds) - 1),
                "median_cds_intron_length_bp": sorted(introns)[len(introns) // 2]
                if introns
                else "",
                "note": ";".join(notes),
                "gff_file": str(path),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-list", required=True, type=Path)
    parser.add_argument("--gene", required=True, choices=tuple(GENE_ALIASES))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    paths = [
        Path(line.strip())
        for line in args.gff_list.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = [row for path in paths for row in summarize(path, args.gene)]
    fields = [
        "query_gene",
        "species_or_file",
        "matched_gene_id",
        "matched_gene_name",
        "transcript_id",
        "seqid",
        "strand",
        "gene_start",
        "gene_end",
        "gene_span_bp",
        "cds_start",
        "cds_end",
        "cds_span_bp",
        "exon_count",
        "cds_exon_count",
        "cds_length_bp",
        "protein_length_est_aa",
        "intron_count_between_cds",
        "median_cds_intron_length_bp",
        "note",
        "gff_file",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} {args.gene} gene-structure rows to {args.out}")


if __name__ == "__main__":
    main()
