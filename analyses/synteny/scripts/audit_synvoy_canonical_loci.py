#!/usr/bin/env python3
"""Compare SynVoy candidate coordinates with canonical protein loci in target GFFs.

This is a rapid confirmation tool for P2/P3 rows. It does not replace the
multi-evidence P1 review: it asks the narrower, reproducible question of whether
the protein accession used in the thesis maps to the same annotated locus as
the best post-ownership SynVoy candidate.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.rstrip().split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def split_accessions(value: str) -> set[str]:
    return {part.strip() for part in value.split(",") if part.strip()}


def collect_protein_intervals(
    gff: Path, wanted: set[str]
) -> dict[str, dict[str, object]]:
    intervals: dict[str, dict[str, object]] = {}
    with gff.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"CDS", "polypeptide"}:
                continue
            attributes = parse_attributes(fields[8])
            accessions: set[str] = set()
            for key in ("protein_id", "Name", "ID"):
                accessions.update(split_accessions(attributes.get(key, "")))
            matched = wanted.intersection(accessions)
            if not matched:
                continue
            seqid, start, end = fields[0], int(fields[3]), int(fields[4])
            for accession in matched:
                record = intervals.setdefault(
                    accession,
                    {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "product": attributes.get("product", ""),
                    },
                )
                if record["seqid"] != seqid:
                    record["multiple_seqids"] = True
                record["start"] = min(int(record["start"]), start)
                record["end"] = max(int(record["end"]), end)
                if not record["product"] and attributes.get("product"):
                    record["product"] = attributes["product"]
    return intervals


def attach_containing_genes(
    gff: Path, intervals: dict[str, dict[str, object]]
) -> None:
    by_seqid: dict[str, list[tuple[str, dict[str, object]]]] = defaultdict(list)
    for accession, record in intervals.items():
        by_seqid[str(record["seqid"])].append((accession, record))
    if not by_seqid:
        return

    with gff.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in by_seqid:
                continue
            start, end = int(fields[3]), int(fields[4])
            attributes = parse_attributes(fields[8])
            for _, record in by_seqid[fields[0]]:
                if start <= int(record["start"]) and end >= int(record["end"]):
                    record["gene_start"] = start
                    record["gene_end"] = end
                    record["gene_symbol"] = (
                        attributes.get("gene")
                        or attributes.get("Name")
                        or attributes.get("locus_tag")
                        or attributes.get("ID", "")
                    )
                    record["gene_id"] = attributes.get("GeneID", "")


def relation(
    candidate_seqid: str,
    candidate_start: str,
    candidate_end: str,
    canonical: dict[str, object] | None,
) -> tuple[str, str]:
    if canonical is None:
        return "canonical accession not found in target GFF", ""
    if not candidate_seqid or not candidate_start or not candidate_end:
        return "no retained SynVoy coordinate", ""
    canonical_seqid = str(canonical["seqid"])
    if canonical_seqid != candidate_seqid:
        return "different sequence", ""
    left, right = int(candidate_start), int(candidate_end)
    canonical_left = int(canonical.get("gene_start", canonical["start"]))
    canonical_right = int(canonical.get("gene_end", canonical["end"]))
    if left <= canonical_right and right >= canonical_left:
        return "same locus", "0"
    distance = max(canonical_left - right, left - canonical_right, 0)
    return "different locus on same sequence", str(distance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-table", required=True, type=Path)
    parser.add_argument("--gff-dir", required=True, type=Path)
    parser.add_argument("--input-audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    evidence = read_tsv(args.evidence_table)
    input_status = {
        row["species"]: row["status"] for row in read_tsv(args.input_audit)
    }
    wanted_by_species: dict[str, set[str]] = defaultdict(set)
    for row in evidence:
        accession = row["canonical_protein_accession"]
        if accession:
            wanted_by_species[row["species"]].add(accession)

    loci_by_species: dict[str, dict[str, dict[str, object]]] = {}
    for species, wanted in sorted(wanted_by_species.items()):
        gff = args.gff_dir / f"{species}.gff"
        if not gff.is_file():
            loci_by_species[species] = {}
            continue
        intervals = collect_protein_intervals(gff, wanted)
        attach_containing_genes(gff, intervals)
        loci_by_species[species] = intervals

    output: list[dict[str, str]] = []
    for row in evidence:
        accession = row["canonical_protein_accession"]
        canonical = loci_by_species.get(row["species"], {}).get(accession)
        coordinate_relation, distance = relation(
            row["best_chromosome"], row["best_start"], row["best_end"], canonical
        )
        target_input_status = input_status.get(row["species"], "UNKNOWN")
        if target_input_status != "PASS":
            conclusion = "invalid input pair; rerun required"
        elif not accession:
            conclusion = "no canonical protein accession available"
        elif canonical is None:
            conclusion = "canonical accession absent from exact target GFF"
        elif coordinate_relation == "same locus":
            conclusion = "confirmed same canonical locus"
        elif row["review_status"]:
            conclusion = f"resolved by detailed review: {row['review_status']}"
        else:
            conclusion = "coordinate conflict; escalate to detailed review"

        output.append(
            {
                "gene": row["gene"],
                "species": row["species"],
                "review_priority": row["review_priority"],
                "input_status": target_input_status,
                "canonical_protein_accession": accession,
                "canonical_seqid": str(canonical.get("seqid", "")) if canonical else "",
                "canonical_start": str(canonical.get("gene_start", canonical.get("start", ""))) if canonical else "",
                "canonical_end": str(canonical.get("gene_end", canonical.get("end", ""))) if canonical else "",
                "canonical_gene_symbol": str(canonical.get("gene_symbol", "")) if canonical else "",
                "canonical_product": str(canonical.get("product", "")) if canonical else "",
                "synvoy_confidence": row["best_confidence"],
                "synvoy_seqid": row["best_chromosome"],
                "synvoy_start": row["best_start"],
                "synvoy_end": row["best_end"],
                "coordinate_relation": coordinate_relation,
                "distance_bp": distance,
                "detailed_review_status": row["review_status"],
                "rapid_conclusion": conclusion,
            }
        )

    fields = list(output[0])
    write_tsv(args.output, output, fields)
    counts: dict[str, int] = defaultdict(int)
    for row in output:
        counts[row["rapid_conclusion"]] += 1
    print(f"Audited {len(output)} gene-by-species rows")
    for conclusion, count in sorted(counts.items()):
        print(f"{count}\t{conclusion}")


if __name__ == "__main__":
    main()
