#!/usr/bin/env python3
"""List pre- and post-ownership SynVoy candidates for priority review rows."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def confidence_rank(value: str) -> int:
    return {"HIGH": 3, "MEDIUM": 2, "LOW": 1}.get(value, 0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", required=True, type=Path)
    parser.add_argument("--synvoy-results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    queue = read_tsv(args.queue)
    by_run: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in queue:
        by_run[row["run"]].append(row)

    output: list[dict[str, object]] = []
    for run, review_rows in sorted(by_run.items()):
        report_path = args.synvoy_results / run / "synvoy_report.json"
        with report_path.open(encoding="utf-8") as handle:
            report = json.load(handle)
        records_by_species: dict[str, list[dict]] = defaultdict(list)
        for record in report["goi_dedup"]["records"]:
            species = Path(record["genome"]).stem
            records_by_species[species].append(record)

        for review in review_rows:
            records = records_by_species.get(review["species"], [])
            records.sort(
                key=lambda item: (
                    confidence_rank(item.get("confidence", "")),
                    float(item.get("identity") or 0),
                ),
                reverse=True,
            )
            for rank, record in enumerate(records, start=1):
                if record.get("confidence") not in {"HIGH", "MEDIUM"}:
                    continue
                output.append(
                    {
                        "gene": review["gene"],
                        "species": review["species"],
                        "rank_within_report": rank,
                        "chromosome": record.get("chrom", ""),
                        "start": record.get("start", ""),
                        "end": record.get("end", ""),
                        "identity_percent": record.get("identity", ""),
                        "confidence": record.get("confidence", ""),
                        "original_goi_class": record.get("original_goi_class", ""),
                        "final_goi_class": record.get("goi_class", ""),
                        "owning_gene": record.get("owning_gene", ""),
                        "inferred_paralog": record.get("inferred_paralog", ""),
                        "model_status": record.get("model_status", ""),
                    }
                )

    if not output:
        raise RuntimeError("No HIGH/MEDIUM priority candidates were found")
    write_tsv(args.output, output)
    print(f"Wrote {len(output)} candidates to {args.output}")


if __name__ == "__main__":
    main()
