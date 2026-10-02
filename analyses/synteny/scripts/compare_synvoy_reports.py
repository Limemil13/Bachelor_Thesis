#!/usr/bin/env python3
"""Compare final GOI candidates in two SynVoy reports, by target genome."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from build_synvoy_review import candidate_is_goi, candidate_sort_key


def candidate_rows(path: Path) -> dict[str, list[dict]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    grouped = defaultdict(list)
    for record in report["goi_dedup"]["records"]:
        if candidate_is_goi(record):
            grouped[Path(record["genome"]).stem].append(record)
    for records in grouped.values():
        records.sort(key=candidate_sort_key, reverse=True)
    return grouped


def describe(records: list[dict]) -> dict[str, str]:
    best = records[0] if records else {}
    return {
        "candidate_count": str(len(records)),
        "high_count": str(sum(r["confidence"] == "HIGH" for r in records)),
        "medium_count": str(sum(r["confidence"] == "MEDIUM" for r in records)),
        "confidence": best.get("confidence", "NONE"),
        "locus": (f"{best['chrom']}:{best['start']}-{best['end']}" if best else ""),
        "identity_percent": best.get("identity", ""),
        "owning_gene": best.get("owning_gene", ""),
        "model_status": best.get("model_status", ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old_report", type=Path)
    parser.add_argument("new_report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    old = candidate_rows(args.old_report)
    new = candidate_rows(args.new_report)
    rows = []
    for species in sorted(set(old) | set(new)):
        before = describe(old.get(species, []))
        after = describe(new.get(species, []))
        rows.append(
            {
                "species": species,
                **{f"old_{key}": value for key, value in before.items()},
                **{f"new_{key}": value for key, value in after.items()},
                "best_candidate_changed": (
                    before["locus"] != after["locus"]
                    or before["confidence"] != after["confidence"]
                    or before["owning_gene"] != after["owning_gene"]
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Compared {len(rows)} genomes; {sum(r['best_candidate_changed'] for r in rows)} best candidates changed."
    )


if __name__ == "__main__":
    main()
