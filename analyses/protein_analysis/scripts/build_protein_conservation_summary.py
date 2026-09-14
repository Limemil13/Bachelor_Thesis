"""Build canonical per-protein and per-gene thesis conservation summaries."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from canonical_dataset import MANIFEST

BASE = Path(__file__).resolve().parents[1]
DOMAIN_IN = (
    BASE
    / "domains"
    / "domain_scan_canonical"
    / "domain_summary_by_protein_canonical.tsv"
)
MSA_IN = BASE / "alignments" / "canonical" / "msa_quality_summary_canonical.tsv"
SIGNALP_IN = BASE / "signal_peptides" / "tables" / "signalp_summary_clean.tsv"
OUT_DIR = BASE / "tables"
PROTEIN_OUT = OUT_DIR / "protein_conservation_domain_msa_signalp_summary.tsv"
GENE_OUT = OUT_DIR / "protein_conservation_domain_msa_signalp_gene_summary.tsv"
REVIEW_OUT = OUT_DIR / "protein_conservation_review_candidates.tsv"

REFERENCE_ACCESSIONS = {
    "BGN": "XP_016885213.1",
    "DCN": "NP_001911.1",
    "EPYC": "NP_004941.2",
    "FMOD": "XP_047272260.1",
    "OGN": "NP_148935.1",
    "PRELP": "NP_958505.1",
    "LUM": "NP_002336.1",
}

DOG_OGN_KEY = ("OGN", "canis_lupus_familiaris", "XP_038383338.1")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["gene"], row["species"], row["accession"]


def index_unique(
    rows: list[dict[str, str]], label: str
) -> dict[tuple[str, str, str], dict[str, str]]:
    indexed: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        row_key = key(row)
        if row_key in indexed:
            raise ValueError(f"Duplicate {label} key: {row_key}")
        indexed[row_key] = row
    return indexed


def msa_quality(raw_status: str) -> str:
    return {
        "ok": "Good",
        "watch": "Watch",
        "check_high_gap": "Needs inspection",
    }.get(raw_status, "Needs inspection")


def domain_quality(raw_status: str, n_lrr: int) -> str:
    if raw_status == "SLRP_like" and n_lrr >= 2:
        return "Good SLRP-like domain architecture"
    if raw_status == "SLRP_like":
        return "Weak but SLRP-like"
    return "Needs inspection"


def canonical_signalp(
    row_key: tuple[str, str, str], raw: dict[str, str]
) -> dict[str, str]:
    result = dict(raw)
    result["signalp_source"] = (
        raw.get("signalp_source") or "SignalP-6.0 better_run (2026-07-27)"
    )
    result["signalp_rerun_required"] = raw.get("signalp_rerun_required") or "no"

    # Preserve the documented fallback only when the corrected canonical dog
    # sequence has not been included in an explicit canonical SignalP run.
    if row_key == DOG_OGN_KEY and "canonical" not in result["signalp_source"].lower():
        result.update(
            {
                "signal_peptide": "yes",
                "signalp_prediction": "SP",
                "other_score": "",
                "sp_score": "",
                "cleavage_site": "19-20",
                "cleavage_probability": "",
                "signalp_quality": "strong_signal_peptide",
                "strict_check": "yes",
                "watch": "yes",
                "notes": (
                    "Corrected 297-aa dog OGN was reported SignalP-positive in "
                    "OGN_final_diagnosis.md; raw corrected prediction record was not retained"
                ),
                "signalp_source": "documented corrected-dog SignalP-6.0 retest",
                "signalp_rerun_required": "yes",
            }
        )
    return result


def qc_status_and_concerns(
    manifest: dict[str, str],
    domain: dict[str, str],
    msa: dict[str, str],
    signalp: dict[str, str],
) -> tuple[str, str]:
    severe: list[str] = []
    watch: list[str] = []

    if domain["domain_status"] != "SLRP_like":
        severe.append("weak_domain_support")
    if msa["msa_status"] == "check_high_gap":
        severe.append("high_alignment_gap")
    elif msa["msa_status"] == "watch":
        watch.append("moderate_alignment_gap")
    if signalp["signal_peptide"] != "yes":
        severe.append("no_signalp_prediction")
    if signalp["signalp_rerun_required"] == "yes":
        severe.append("canonical_signalp_raw_result_missing")

    qcov = float(manifest["forward_qcovs"]) if manifest["forward_qcovs"] else 0.0
    pident = float(manifest["forward_pident"]) if manifest["forward_pident"] else 0.0
    if qcov < 70:
        severe.append("low_query_coverage")
    elif qcov < 85:
        watch.append("reduced_query_coverage")
    if pident < 40:
        watch.append("low_forward_identity")
    if manifest["reverse_hit"] != REFERENCE_ACCESSIONS[manifest["gene"]]:
        watch.append("reciprocal_hit_differs_from_reference")
    if manifest["sequence_correction"] != "none":
        watch.append("sequence_model_corrected")
    if signalp["watch"] == "yes" and signalp["signal_peptide"] == "yes":
        watch.append("signalp_watch")

    concerns = severe + [item for item in watch if item not in severe]
    if "canonical_signalp_raw_result_missing" in severe:
        status = "Needs SignalP rerun"
    elif severe:
        status = "Needs inspection"
    elif watch:
        status = "Watch"
    else:
        status = "Supported"
    return status, ";".join(concerns)


def main() -> None:
    manifest_rows = read_tsv(MANIFEST)
    domain_rows = read_tsv(DOMAIN_IN)
    msa_rows = read_tsv(MSA_IN)
    signalp_rows = read_tsv(SIGNALP_IN)

    manifest = index_unique(manifest_rows, "manifest")
    domains = index_unique(domain_rows, "domain")
    msa = index_unique(msa_rows, "MSA")
    signalp = index_unique(signalp_rows, "SignalP")

    expected_keys = set(manifest)
    for label, indexed in (("domain", domains), ("MSA", msa)):
        if set(indexed) != expected_keys:
            missing = sorted(expected_keys - set(indexed))
            extra = sorted(set(indexed) - expected_keys)
            raise ValueError(
                f"{label} keys differ from manifest; missing={missing}, extra={extra}"
            )

    # A promoted gene may be added to the canonical computational panel before
    # a new SignalP web/standalone run is available. Keep that gap explicit
    # rather than copying or inferring a prediction. Historical SignalP rows
    # that are no longer canonical are ignored, while missing rows become
    # honest pending records below.
    unexpected_signalp = sorted(set(signalp) - expected_keys)
    if unexpected_signalp:
        raise ValueError(
            f"SignalP contains unexpected canonical keys: {unexpected_signalp}"
        )

    rows: list[dict[str, object]] = []
    for manifest_row in manifest_rows:
        row_key = key(manifest_row)
        domain = domains[row_key]
        alignment = msa[row_key]
        raw_signalp = signalp.get(
            row_key,
            {
                "signal_peptide": "not_run",
                "signalp_prediction": "pending",
                "other_score": "",
                "sp_score": "",
                "cleavage_site": "",
                "cleavage_probability": "",
                "signalp_quality": "pending",
                "strict_check": "yes",
                "watch": "yes",
                "notes": (
                    f"Canonical {manifest_row['gene']} sequence prepared; "
                    "SignalP prediction not yet run"
                ),
                "signalp_source": f"pending {manifest_row['gene']} canonical SignalP run",
                "signalp_rerun_required": "yes",
            },
        )
        sp = canonical_signalp(row_key, raw_signalp)
        n_lrr = int(domain["n_lrr_domains"])
        status, concerns = qc_status_and_concerns(manifest_row, domain, alignment, sp)

        rows.append(
            {
                "Gene": manifest_row["gene"],
                "Species": manifest_row["species"],
                "Protein accession": manifest_row["accession"],
                "Original protein length": manifest_row["original_length"],
                "Analyzed protein length": manifest_row["analyzed_length"],
                "Sequence correction": manifest_row["sequence_correction"],
                "Forward identity percent": manifest_row["forward_pident"],
                "Forward query coverage percent": manifest_row["forward_qcovs"],
                "Reciprocal hit": manifest_row["reverse_hit"],
                "Reciprocal reference match": (
                    "yes"
                    if manifest_row["reverse_hit"]
                    == REFERENCE_ACCESSIONS[manifest_row["gene"]]
                    else "no"
                ),
                "Alignment file": alignment["alignment_file"],
                "Alignment length": alignment["alignment_length"],
                "Gap percentage": alignment["gap_percent"],
                "Mean pairwise identity percent": alignment[
                    "mean_pairwise_identity_percent"
                ],
                "Mean pairwise comparable sites": alignment[
                    "mean_pairwise_comparable_sites"
                ],
                "MSA quality": msa_quality(alignment["msa_status"]),
                "Number of significant Pfam domains": domain["n_significant_domains"],
                "Number of LRR-related domains": domain["n_lrr_domains"],
                "Detected domain families": domain["domain_names"],
                "Domain quality": domain_quality(domain["domain_status"], n_lrr),
                "Signal peptide": sp["signal_peptide"],
                "SignalP prediction": sp["signalp_prediction"],
                "SignalP SP score": sp["sp_score"],
                "SignalP cleavage site": sp["cleavage_site"],
                "SignalP cleavage probability": sp["cleavage_probability"],
                "SignalP quality": sp["signalp_quality"],
                "SignalP source": sp["signalp_source"],
                "SignalP rerun required": sp["signalp_rerun_required"],
                "Protein QC status": status,
                "QC concerns": concerns,
                "Notes": sp["notes"],
            }
        )

    protein_fields = list(rows[0])
    write_tsv(PROTEIN_OUT, rows, protein_fields)
    review_rows = [row for row in rows if row["Protein QC status"] != "Supported"]
    write_tsv(REVIEW_OUT, review_rows, protein_fields)

    by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_gene[str(row["Gene"])].append(row)

    gene_rows: list[dict[str, object]] = []
    for gene in ("BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM"):
        gene_proteins = by_gene[gene]
        lengths = [int(row["Analyzed protein length"]) for row in gene_proteins]
        identities = [float(row["Forward identity percent"]) for row in gene_proteins]
        coverages = [
            float(row["Forward query coverage percent"]) for row in gene_proteins
        ]
        msa_identities = [
            float(row["Mean pairwise identity percent"]) for row in gene_proteins
        ]
        domain_good = sum(
            row["Domain quality"] != "Needs inspection" for row in gene_proteins
        )
        signalp_positive = sum(row["Signal peptide"] == "yes" for row in gene_proteins)
        signalp_pending = sum(
            row["SignalP prediction"] == "pending" for row in gene_proteins
        )
        signalp_tested = len(gene_proteins) - signalp_pending
        needs_inspection = sum(
            row["Protein QC status"] in {"Needs inspection", "Needs SignalP rerun"}
            for row in gene_proteins
        )
        watches = sum(row["Protein QC status"] == "Watch" for row in gene_proteins)

        gene_rows.append(
            {
                "Gene": gene,
                "Protein count": len(gene_proteins),
                "Analyzed length minimum": min(lengths),
                "Analyzed length median": round(median(lengths), 1),
                "Analyzed length maximum": max(lengths),
                "Mean forward identity percent": round(mean(identities), 2),
                "Minimum forward identity percent": round(min(identities), 2),
                "Mean forward query coverage percent": round(mean(coverages), 2),
                "Minimum forward query coverage percent": round(min(coverages), 2),
                "Mean within-gene pairwise identity percent": round(
                    mean(msa_identities), 2
                ),
                "SLRP-like domain count": domain_good,
                "SignalP-positive count": signalp_positive,
                "SignalP-tested count": signalp_tested,
                "SignalP-pending count": signalp_pending,
                "Corrected sequence count": sum(
                    row["Sequence correction"] != "none" for row in gene_proteins
                ),
                "Watch count": watches,
                "Needs-inspection count": needs_inspection,
                "Data summary": (
                    f"{domain_good}/{len(gene_proteins)} SLRP-like domain support; "
                    f"{signalp_positive}/{signalp_tested} tested SignalP-positive; "
                    f"{signalp_pending} SignalP pending; "
                    f"{needs_inspection} require inspection/rerun"
                ),
            }
        )

    write_tsv(GENE_OUT, gene_rows, list(gene_rows[0]))
    print(f"Wrote {PROTEIN_OUT} ({len(rows)} proteins)")
    print(f"Wrote {GENE_OUT} ({len(gene_rows)} genes)")
    print(f"Wrote {REVIEW_OUT} ({len(review_rows)} review/watch proteins)")


if __name__ == "__main__":
    main()
