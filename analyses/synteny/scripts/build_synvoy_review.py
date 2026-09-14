#!/usr/bin/env python3
"""Build a compact, evidence-based review queue from the final SynVoy runs.

The SynVoy report contains many low-confidence family matches and also retains
HIGH/MEDIUM records that its ownership or coverage checks demote. This script
therefore uses the report's final definition of a GOI candidate: HIGH/MEDIUM
records excluding paralogs and coverage-demoted calls. It never turns a SynVoy
candidate into a claimed ortholog without manual/NCBI confirmation.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

RUNS = {
    "BGN": "bgn_human_15species_dev_20260905",
    "DCN": "dcn_human_15species_dev_20260903",
    "EPYC": "epyc_human_15species",
    "FMOD": "fmod_human_15species",
    "OGN": "ogn_human_15species_standard_1408",
    "PRELP": "prelp_human_15species_fresh",
    "LUM": "lum_human_15species",
}

SPECIES_TO_CANONICAL = {
    "mouse": "mouse",
    "cow": "bos_taurus",
    "dog": "canis_lupus_familiaris",
    "opossum": "monodelphis_domestica",
    "chicken": "chicken",
    "anole": "anolis_carolinensis",
    "frog": "frog",
    "zebrafish": "zebrafish",
    "spotted_gar": "lepisosteus_oculatus",
    "coelacanth": "latimeria_chalumnae",
    "elephant_shark": "callorhinchus_milii",
    "catshark": "scyliorhinus_canicula",
    "amphioxus": "amphioxus",
    "whale_shark": "whale_shark",
}

DEEP_LINEAGES = {
    "amphioxus",
    "catshark",
    "elephant_shark",
    "whale_shark",
    "coelacanth",
    "spotted_gar",
}

MANUAL_FIELDS = (
    "coordinate_to_accession_status",
    "review_status",
    "confirmed_gene_symbol",
    "confirmed_protein_accession",
    "neighbor_order_consistent",
    "phylogeny_consistent",
    "reviewer",
    "review_date",
    "manual_notes",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def number(value, default=0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def candidate_is_goi(record: dict) -> bool:
    return (
        record.get("confidence") in {"HIGH", "MEDIUM"}
        and record.get("goi_class") != "paralog_not_goi"
        and not record.get("coverage_demoted", False)
    )


def candidate_sort_key(record: dict) -> tuple:
    confidence_rank = 2 if record.get("confidence") == "HIGH" else 1
    complete_rank = 1 if record.get("model_status") == "complete" else 0
    return (confidence_rank, number(record.get("identity")), complete_rank)


def load_protein_support(
    conservation_path: Path, protspace_path: Path
) -> dict[tuple[str, str], dict]:
    support: dict[tuple[str, str], dict] = {}
    for row in read_tsv(conservation_path):
        key = (row["Gene"], row["Species"])
        support[key] = {
            "canonical_accession": row["Protein accession"],
            "protein_qc_status": row["Protein QC status"],
            "protein_qc_concerns": row["QC concerns"],
            "signalp_quality": row["SignalP quality"],
            "signalp_prediction": row["SignalP prediction"],
        }
    for row in read_tsv(protspace_path):
        species, gene, accession = row["species_accession"].split("|", 2)
        key = (gene, species)
        item = support.setdefault(key, {"canonical_accession": accession})
        item["protspace_status"] = row["embedding_qc_status"]
        item["protspace_concerns"] = row["embedding_qc_concerns"]
    return support


def classify_review(
    valid: list[dict],
    best: dict | None,
    flag_count: int,
    protein: dict,
    species: str,
) -> tuple[str, str, str, str]:
    reasons: list[str] = []
    if not valid:
        reasons.append("no post-ownership HIGH/MEDIUM GOI candidate")
        return "P1", "D", "detailed-review", "; ".join(reasons)

    high = [x for x in valid if x.get("confidence") == "HIGH"]
    if not high:
        reasons.append("best retained candidate is MEDIUM")
    if len(valid) > 1:
        reasons.append(f"{len(valid)} retained GOI candidates")
    if len(high) > 1:
        reasons.append(f"{len(high)} HIGH candidates")
    if best and best.get("model_status") != "complete":
        reasons.append(f"best model status is {best.get('model_status') or 'unknown'}")
    if best and best.get("goi_class") == "synteny_hull_rescue":
        reasons.append("best call is a synteny-hull rescue")
    if flag_count:
        reasons.append(f"{flag_count} self-consistency flag(s) for species")
    if protein.get("protein_qc_status") not in {None, "", "Supported"}:
        reasons.append(f"protein QC {protein.get('protein_qc_status')}")
    if protein.get("protspace_status") not in {None, "", "Supported"}:
        reasons.append(f"ProtSpace {protein.get('protspace_status')}")
    if species in DEEP_LINEAGES:
        reasons.append("deep-lineage annotation")

    strong = (
        len(valid) == 1
        and len(high) == 1
        and best is not None
        and best.get("model_status") == "complete"
        and not flag_count
        and protein.get("protein_qc_status") in {None, "", "Supported"}
        and protein.get("protspace_status") in {None, "", "Supported"}
        and species not in DEEP_LINEAGES
    )
    if strong:
        return (
            "P3",
            "A",
            "spot-check",
            "single complete HIGH call; no automated exception",
        )
    if not high:
        return "P1", "C", "detailed-review", "; ".join(reasons)
    return "P2", "B", "rapid-confirmation", "; ".join(reasons)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synvoy-results", required=True, type=Path)
    parser.add_argument("--conservation-table", required=True, type=Path)
    parser.add_argument("--protspace-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--input-audit",
        type=Path,
        help="Optional FASTA/GFF sequence-ID audit keyed by species.",
    )
    parser.add_argument(
        "--manual-decisions",
        type=Path,
        help="Optional reviewed manual fields keyed by gene and species.",
    )
    args = parser.parse_args()

    protein_support = load_protein_support(
        args.conservation_table, args.protspace_table
    )
    input_pair_qc: dict[str, str] = {}
    if args.input_audit:
        input_pair_qc = {
            row["species"]: row["status"] for row in read_tsv(args.input_audit)
        }
    # Preserve completed review fields when the automated table is rebuilt
    # after a new or updated SynVoy report.
    previous_manual: dict[tuple[str, str], dict[str, str]] = {}
    previous_table = args.output_dir / "synvoy_gene_species_evidence.tsv"
    if previous_table.is_file():
        for old in read_tsv(previous_table):
            previous_manual[(old["gene"], old["species"])] = {
                field: old.get(field, "") for field in MANUAL_FIELDS
            }
    reviewed_manual: dict[tuple[str, str], dict[str, str]] = {}
    if args.manual_decisions:
        for reviewed in read_tsv(args.manual_decisions):
            key = (reviewed["gene"], reviewed["species"])
            if key in reviewed_manual:
                raise ValueError(f"Duplicate manual decision for {key}")
            reviewed_manual[key] = {
                field: reviewed.get(field, "") for field in MANUAL_FIELDS
            }
    rows: list[dict] = []
    run_rows: list[dict] = []

    for gene, run_name in RUNS.items():
        report_path = args.synvoy_results / run_name / "synvoy_report.json"
        with report_path.open(encoding="utf-8") as handle:
            report = json.load(handle)

        summary = report["summary"]
        qc = report["qc_summary"]
        consistency = report["self_consistency"]["summary"]
        run_rows.append(
            {
                "gene": gene,
                "run": run_name,
                "report": f"results/{run_name}/synvoy_report.json",
                "genomes": qc["total_genomes"],
                "genomes_qc_pass": qc["pass"],
                "final_high_goi": summary["high_confidence_goi"],
                "final_medium_goi": summary["medium_confidence_goi"],
                "low_ambiguous_goi": summary["low_confidence_goi"],
                "self_consistency_flags": consistency["total_flags"],
                "paralog_misassignments_removed": consistency.get(
                    "n_paralog_misassignment", 0
                ),
                "headline": summary["headline"],
            }
        )

        annotations = {x["genome"]: x for x in report["annotations"]["per_genome"]}
        regions = {x["genome"]: x for x in report["regions"]["per_genome"]}
        genome_qc = {Path(x["genome"]).stem: x for x in report["genome_qc"]}
        flags_by_species: dict[str, list[dict]] = defaultdict(list)
        for flag in report["self_consistency"].get("flags", []):
            flags_by_species[Path(flag["genome"]).stem].append(flag)
        records_by_species: dict[str, list[dict]] = defaultdict(list)
        for record in report["goi_dedup"]["records"]:
            records_by_species[Path(record["genome"]).stem].append(record)

        final_valid = [
            record
            for records in records_by_species.values()
            for record in records
            if candidate_is_goi(record)
        ]
        assert (
            sum(x["confidence"] == "HIGH" for x in final_valid)
            == summary["high_confidence_goi"]
        )
        assert (
            sum(x["confidence"] == "MEDIUM" for x in final_valid)
            == summary["medium_confidence_goi"]
        )

        for species in sorted(annotations):
            annotation = annotations[species]
            valid = [x for x in records_by_species[species] if candidate_is_goi(x)]
            valid.sort(key=candidate_sort_key, reverse=True)
            best = valid[0] if valid else None
            high_count = sum(x.get("confidence") == "HIGH" for x in valid)
            medium_count = sum(x.get("confidence") == "MEDIUM" for x in valid)
            raw_counts = annotation.get("goi_confidence_counts", {})
            canonical_species = SPECIES_TO_CANONICAL[species]
            protein = protein_support.get((gene, canonical_species), {})
            species_flags = flags_by_species[species]
            priority, grade, action, reasons = classify_review(
                valid, best, len(species_flags), protein, species
            )
            region = regions.get(species, {})
            qc_item = genome_qc.get(species, {})

            row = {
                "gene": gene,
                "species": species,
                "run": run_name,
                "genome_qc": qc_item.get("status", "UNKNOWN"),
                "input_pair_qc": input_pair_qc.get(species, "not audited"),
                "raw_high_annotations": raw_counts.get("HIGH", 0),
                "raw_medium_annotations": raw_counts.get("MEDIUM", 0),
                "raw_low_ambiguous_annotations": raw_counts.get("LOW", 0),
                "retained_high_goi_candidates": high_count,
                "retained_medium_goi_candidates": medium_count,
                "retained_goi_candidates": len(valid),
                "best_confidence": best.get("confidence", "NONE") if best else "NONE",
                "best_goi_class": best.get("goi_class", "") if best else "",
                "best_model_status": best.get("model_status", "") if best else "",
                "best_identity_percent": best.get("identity", "") if best else "",
                "best_chromosome": best.get("chrom", "") if best else "",
                "best_start": best.get("start", "") if best else "",
                "best_end": best.get("end", "") if best else "",
                "best_owning_gene": best.get("owning_gene", "") if best else "",
                "best_ownership_gap": best.get("ownership_gap", "") if best else "",
                "best_region_score": region.get("best_score", ""),
                "self_consistency_flags": len(species_flags),
                "self_consistency_flag_types": ",".join(
                    sorted({x["type"] for x in species_flags})
                ),
                "canonical_protein_accession": protein.get("canonical_accession", ""),
                "protein_qc_status": protein.get("protein_qc_status", "not available"),
                "protein_qc_concerns": protein.get("protein_qc_concerns", ""),
                "signalp_prediction": protein.get("signalp_prediction", ""),
                "signalp_quality": protein.get("signalp_quality", ""),
                "protspace_status": protein.get("protspace_status", "not available"),
                "protspace_concerns": protein.get("protspace_concerns", ""),
                "evidence_grade": grade,
                "review_priority": priority,
                "manual_action": action,
                "review_reasons": reasons,
                "coordinate_to_accession_status": "manual NCBI/GFF confirmation required",
                "review_status": "",
                "confirmed_gene_symbol": "",
                "confirmed_protein_accession": "",
                "neighbor_order_consistent": "",
                "phylogeny_consistent": "",
                "reviewer": "",
                "review_date": "",
                "manual_notes": "",
            }
            old_manual = previous_manual.get((gene, species), {})
            for field in MANUAL_FIELDS:
                if old_manual.get(field):
                    row[field] = old_manual[field]
            reviewed = reviewed_manual.get((gene, species), {})
            for field in MANUAL_FIELDS:
                if reviewed.get(field):
                    row[field] = reviewed[field]
            rows.append(row)

    priority_rank = {"P1": 1, "P2": 2, "P3": 3}
    rows.sort(
        key=lambda x: (priority_rank[x["review_priority"]], x["gene"], x["species"])
    )
    fields = list(rows[0])
    write_tsv(args.output_dir / "synvoy_gene_species_evidence.tsv", rows, fields)
    detailed_rows = [x for x in rows if x["manual_action"] == "detailed-review"]
    pending_detailed = [x for x in detailed_rows if not x["review_status"]]
    completed_detailed = [x for x in detailed_rows if x["review_status"]]
    write_tsv(
        args.output_dir / "synvoy_manual_review_queue.tsv",
        pending_detailed,
        fields,
    )
    write_tsv(
        args.output_dir / "synvoy_manual_review_completed.tsv",
        completed_detailed,
        fields,
    )
    rapid_rows = [x for x in rows if x["manual_action"] == "rapid-confirmation"]
    pending_rapid = [x for x in rapid_rows if not x["review_status"]]
    completed_rapid = [x for x in rapid_rows if x["review_status"]]
    write_tsv(
        args.output_dir / "synvoy_batch_confirmation_queue.tsv",
        pending_rapid,
        fields,
    )
    write_tsv(
        args.output_dir / "synvoy_batch_confirmation_completed.tsv",
        completed_rapid,
        fields,
    )
    spot_rows = [x for x in rows if x["manual_action"] == "spot-check"]
    pending_spots = [x for x in spot_rows if not x["review_status"]]
    completed_spots = [x for x in spot_rows if x["review_status"]]
    write_tsv(
        args.output_dir / "synvoy_spot_check_queue.tsv",
        pending_spots,
        fields,
    )
    write_tsv(
        args.output_dir / "synvoy_spot_check_completed.tsv",
        completed_spots,
        fields,
    )
    write_tsv(args.output_dir / "synvoy_run_summary.tsv", run_rows, list(run_rows[0]))

    counts = Counter((x["review_priority"], x["evidence_grade"]) for x in rows)
    print(f"Validated {len(rows)} gene x species rows from {len(RUNS)} final reports")
    for (priority, grade), count in sorted(counts.items()):
        print(f"{priority} / grade {grade}: {count}")
    print(f"Detailed exceptions completed: {len(completed_detailed)}")
    print(f"Detailed exceptions pending: {len(pending_detailed)}")
    print(f"Rapid confirmations completed: {len(completed_rapid)}")
    print(f"Rapid confirmations pending: {len(pending_rapid)}")
    print(f"Spot checks completed: {len(completed_spots)}")
    print(f"Spot checks pending: {len(pending_spots)}")


if __name__ == "__main__":
    main()
