"""Create a conservative musculoskeletal summary from fetched Bgee calls."""

from __future__ import annotations

import csv
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
TABLE_DIR = BASE / "tables"
CALLS = TABLE_DIR / "bgee_slrp_expression_calls_all.tsv"
MAPPING = TABLE_DIR / "bgee_slrp_gene_mapping.tsv"

DIRECT = re.compile(r"growth plate|cartilage|chondr", re.IGNORECASE)
BONE_JOINT = re.compile(
    r"\bbone\b|bone element|epiphys|\bjoint\b|femur|tibia|humerus|condyle",
    re.IGNORECASE,
)
DEVELOPMENT = re.compile(
    r"skeletal system|notochord|sclerotome|somite|intervertebral|vertebr|tendon|ligament",
    re.IGNORECASE,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def category(name: str) -> str:
    # Assign anatomy labels to broad evidence classes. Muscle and bone marrow
    # are excluded because their names can otherwise match broader terms.
    lowered = name.casefold()
    if "skeletal muscle" in lowered or "bone marrow" in lowered:
        return "none"
    if DIRECT.search(name):
        return "direct_cartilage_or_growth_plate"
    if BONE_JOINT.search(name):
        return "bone_or_joint"
    if DEVELOPMENT.search(name):
        return "skeletal_development_or_connective_tissue"
    return "none"


def main() -> None:
    # First retain musculoskeletal calls, then summarize their presence for each
    # requested species, gene and query mode.
    calls = read_tsv(CALLS)
    mappings = read_tsv(MAPPING)
    # Preserve individual Bgee calls so the summary remains auditable.
    relevant: list[dict[str, str]] = []
    for row in calls:
        label = category(row["anatomical_entity_name"])
        if label != "none":
            relevant.append({**row, "relevance_category": label})
    relevant.sort(
        key=lambda row: (
            row["species"],
            row["thesis_gene"],
            row["query_mode"],
            {
                "direct_cartilage_or_growth_plate": 0,
                "bone_or_joint": 1,
                "skeletal_development_or_connective_tissue": 2,
            }[row["relevance_category"]],
            -float(row["expression_score"] or 0),
        )
    )
    write_tsv(TABLE_DIR / "bgee_slrp_relevant_expression_calls.tsv", relevant)

    # An empty subset means no matching downloaded call, not biological absence.
    summary: list[dict[str, str]] = []
    for mapping in mappings:
        for mode in ("anatomy_all_data", "stage_rnaseq"):
            subset = [
                row
                for row in relevant
                if row["species"] == mapping["species"]
                and row["thesis_gene"] == mapping["thesis_gene"]
                and row["query_mode"] == mode
            ]
            by_category = {
                label: [row for row in subset if row["relevance_category"] == label]
                for label in (
                    "direct_cartilage_or_growth_plate",
                    "bone_or_joint",
                    "skeletal_development_or_connective_tissue",
                )
            }
            # Prefer direct cartilage calls, then broader bone and developmental
            # terms. The highest expression score is selected within that class.
            prioritized = (
                by_category["direct_cartilage_or_growth_plate"]
                or by_category["bone_or_joint"]
                or by_category["skeletal_development_or_connective_tissue"]
            )
            best = max(
                prioritized,
                key=lambda row: float(row["expression_score"] or 0),
                default=None,
            )
            if not mapping["bgee_gene_id"]:
                interpretation = "unresolved Bgee gene mapping; missing evidence, not biological absence"
            elif best:
                interpretation = "curated relevant expression support found"
            else:
                interpretation = (
                    "no relevant Bgee call; missing evidence, not biological absence"
                )
            summary.append(
                {
                    "species": mapping["species"],
                    "thesis_gene": mapping["thesis_gene"],
                    "bgee_gene_id": mapping["bgee_gene_id"],
                    "query_mode": mode,
                    "direct_cartilage_or_growth_plate_call_count": str(
                        len(by_category["direct_cartilage_or_growth_plate"])
                    ),
                    "bone_or_joint_call_count": str(len(by_category["bone_or_joint"])),
                    "skeletal_development_call_count": str(
                        len(by_category["skeletal_development_or_connective_tissue"])
                    ),
                    "best_priority_category": best["relevance_category"]
                    if best
                    else "",
                    "best_priority_term": best["anatomical_entity_name"]
                    if best
                    else "",
                    "best_priority_stage": best["developmental_stage_name"]
                    if best
                    else "",
                    "best_priority_expression_score": best["expression_score"]
                    if best
                    else "",
                    "best_priority_quality": best["expression_quality"] if best else "",
                    "best_priority_data_types": best["supporting_data_types"]
                    if best
                    else "",
                    "interpretation": interpretation,
                }
            )
    write_tsv(TABLE_DIR / "bgee_slrp_relevant_expression_summary.tsv", summary)
    print(
        f"Relevant calls after excluding skeletal muscle/bone marrow: {len(relevant)}"
    )
    print(f"Summary rows: {len(summary)}")


if __name__ == "__main__":
    main()
