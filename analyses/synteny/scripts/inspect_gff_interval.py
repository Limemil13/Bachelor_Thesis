#!/usr/bin/env python3
"""Inspect annotation features and nearby genes around one genomic interval."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.parse import unquote


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def symbol(attributes: dict[str, str]) -> str:
    value = (
        attributes.get("gene")
        or attributes.get("Name")
        or attributes.get("locus_tag")
        or attributes.get("ID", "unnamed")
    )
    return value.removeprefix("gene-")


def overlaps(start: int, end: int, query_start: int, query_end: int) -> bool:
    return start <= query_end and end >= query_start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff", type=Path)
    parser.add_argument("seqid")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("--flank-count", type=int, default=10)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    genes: list[dict[str, object]] = []
    interval_features: list[dict[str, object]] = []
    with args.gff.open(encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != args.seqid:
                continue
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            attributes = parse_attributes(fields[8])
            record = {
                "seqid": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": feature_start,
                "end": feature_end,
                "strand": fields[6],
                "symbol": symbol(attributes),
                "id": attributes.get("ID", ""),
                "parent": attributes.get("Parent", ""),
                "description": attributes.get("description", ""),
                "product": attributes.get("product", ""),
                "protein_id": attributes.get("protein_id", ""),
            }
            if fields[2] == "gene":
                genes.append(record.copy())
            if overlaps(feature_start, feature_end, args.start, args.end):
                interval_features.append(record.copy())

    genes.sort(key=lambda item: (int(item["start"]), int(item["end"])))
    before = [gene for gene in genes if int(gene["end"]) < args.start]
    after = [gene for gene in genes if int(gene["start"]) > args.end]
    overlapping_genes = [
        gene
        for gene in genes
        if overlaps(int(gene["start"]), int(gene["end"]), args.start, args.end)
    ]
    neighborhood = (
        before[-args.flank_count :] + overlapping_genes + after[: args.flank_count]
    )
    for gene in neighborhood:
        if overlaps(int(gene["start"]), int(gene["end"]), args.start, args.end):
            gene["role"] = "overlap"
        elif int(gene["end"]) < args.start:
            gene["role"] = "upstream"
        else:
            gene["role"] = "downstream"

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {
        args.output_prefix.with_name(args.output_prefix.name + "_features.tsv"): interval_features,
        args.output_prefix.with_name(args.output_prefix.name + "_neighbors.tsv"): neighborhood,
    }
    for path, rows in outputs.items():
        fields = [
            "seqid",
            "source",
            "feature",
            "start",
            "end",
            "strand",
            "symbol",
            "id",
            "parent",
            "description",
            "product",
            "protein_id",
        ]
        if rows and "role" in rows[0]:
            fields.append("role")
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {path}")

    overlap_names = ", ".join(str(item["symbol"]) for item in overlapping_genes)
    print(f"Overlapping genes: {overlap_names or 'none'}")


if __name__ == "__main__":
    main()
