#!/usr/bin/env python3
"""Build one thesis-facing evidence row per main-panel SLRP gene."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "OGN", "LUM")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def indexed(path: Path, field: str) -> dict[str, dict[str, str]]:
    return {row[field]: row for row in read_tsv(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-structure-summary", required=True, type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "analyses/overview/six_gene_integrated_summary.tsv",
    )
    args = parser.parse_args()

    candidate = indexed(ROOT / "analyses/overview/candidate_gene_evidence.tsv", "gene")
    protein = indexed(
        ROOT
        / "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_gene_summary.tsv",
        "Gene",
    )
    msa = indexed(
        ROOT / "analyses/protein_analysis/tables/msa_conservation_gene_summary.tsv",
        "gene",
    )
    protspace = indexed(
        ROOT
        / "analyses/protein_analysis/protspace/tables/protspace_canonical_gene_summary.tsv",
        "gene",
    )
    tree = indexed(
        ROOT
        / "analyses/phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv",
        "gene",
    )
    structure = indexed(args.gene_structure_summary, "query_gene")
    synvoy = indexed(ROOT / "analyses/synteny/tables/synvoy_run_summary.tsv", "gene")

    rows: list[dict[str, str]] = []
    for gene in GENES:
        c, p, m, e, t, s = (
            candidate[gene],
            protein[gene],
            msa[gene],
            protspace[gene],
            tree[gene],
            structure[gene],
        )
        sy = synvoy.get(gene)
        rows.append(
            {
                "gene": gene,
                "slrp_class": c["slrp_class"],
                "recommended_role": c["recommended_role"],
                "decision_note": c["decision_note"],
                "local_p0_mean_tpm": c["local_p0_mean_tpm"],
                "local_p0_rank_of_9": c["local_p0_rank_of_9"],
                "mouse_growth_plate_mean_rank": c["mouse_growth_plate_mean_rank"],
                "mouse_conditions_all_replicates": c[
                    "mouse_growth_plate_all_replicates_conditions"
                ],
                "rat_growth_plate_mean_rank": c["rat_growth_plate_mean_rank"],
                "rat_conditions_all_replicates": c[
                    "rat_growth_plate_all_replicates_conditions"
                ],
                "bgee_mouse_direct_cartilage_calls": c[
                    "bgee_mouse_direct_cartilage_call_count"
                ],
                "literature_strength": c["literature_strength"],
                "literature_note": c["literature_note"],
                "literature_url": c["literature_url"],
                "canonical_protein_count": p["Protein count"],
                "mean_forward_identity_percent": p["Mean forward identity percent"],
                "minimum_forward_identity_percent": p[
                    "Minimum forward identity percent"
                ],
                "mean_forward_query_coverage_percent": p[
                    "Mean forward query coverage percent"
                ],
                "mean_pairwise_msa_identity_percent": m[
                    "mean_pairwise_identity_percent"
                ],
                "fully_conserved_ungapped_columns": m[
                    "fully_conserved_ungapped_columns"
                ],
                "slrp_like_domain_count": p["SLRP-like domain count"],
                "signalp_positive_count": p["SignalP-positive count"],
                "signalp_tested_count": p["SignalP-tested count"],
                "signalp_pending_count": p["SignalP-pending count"],
                "protein_watch_count": p["Watch count"],
                "protein_needs_inspection_count": p["Needs-inspection count"],
                "protspace_centroid_correct": e["centroid_assignment_correct_count"],
                "protspace_supported_count": e["supported_count"],
                "protspace_watch_count": e["watch_count"],
                "protspace_needs_inspection_count": e["needs_inspection_count"],
                "gene_structure_species_count": s["n_species"],
                "gene_structure_cds_exon_counts": s["cds_exon_counts"],
                "gene_structure_protein_length_range": f"{s['protein_length_min']}-{s['protein_length_max']}",
                "gene_structure_span_range_bp": f"{s['gene_span_min_bp']}-{s['gene_span_max_bp']}",
                "gene_structure_outliers": s["structural_outlier_species"],
                "compact_tree_unrooted_monophyletic": t["unrooted_monophyletic"],
                "compact_tree_support_label": t["supporting_split_label"],
                "compact_tree_largest_pure_gene_split": t[
                    "largest_pure_gene_split_tip_count"
                ],
                "compact_tree_outside_main_gene_split": t[
                    "tips_outside_largest_pure_gene_split"
                ],
                "synvoy_status": "complete" if sy else "not available",
                "synvoy_final_high_goi": sy["final_high_goi"] if sy else "",
                "synvoy_final_medium_goi": sy["final_medium_goi"] if sy else "",
                "synvoy_self_consistency_flags": sy["self_consistency_flags"]
                if sy
                else "",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output} ({len(rows)} genes)")


if __name__ == "__main__":
    main()
