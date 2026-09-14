#!/usr/bin/env python3
"""Integrate expression, literature, and completed-pipeline evidence for SLRP selection."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analyses" / "overview" / "candidate_gene_evidence.tsv"

GENES = ("BGN", "FMOD", "PRELP", "EPYC", "OGN", "LUM", "DCN", "OMD", "ASPN")

ANNOTATION = {
    "BGN": {
        "slrp_class": "I",
        "literature_strength": "strong direct skeletal/growth-plate evidence",
        "literature_note": "Detected in growth plate; knockout mice show reduced growth and bone mass.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/9731537/",
        "decision": "main",
        "decision_note": "Strongest overall anchor: very high local P0 expression, cross-dataset growth-plate support, phenotype, and a 13-protein decontaminated canonical panel; the excluded chicken BGN-locus/asporin product remains a synteny-level annotation conflict.",
    },
    "FMOD": {
        "slrp_class": "II",
        "literature_strength": "strong direct growth-plate evidence",
        "literature_note": "Developmental knee study localizes FMOD in growth plate and late-hypertrophic chondrocytes.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/11311118/",
        "decision": "main",
        "decision_note": "Strong cross-species growth-plate expression and direct localization; technically well supported.",
    },
    "PRELP": {
        "slrp_class": "II",
        "literature_strength": "strong zone-specific growth-plate evidence",
        "literature_note": "Validated as a proliferative/columnar growth-plate-zone marker in juvenile rodents.",
        "literature_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3418671/",
        "decision": "main",
        "decision_note": "Independent zone-resolved datasets strongly support PRELP despite only moderate local P0 abundance.",
    },
    "EPYC": {
        "slrp_class": "III",
        "literature_strength": "strong direct epiphyseal/growth-plate evidence",
        "literature_note": "Expressed in growth plate; single knockout has mild long-bone/OA phenotype and interacts genetically with BGN.",
        "literature_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3013283/",
        "decision": "main",
        "decision_note": "Biologically specific class-III representative; expression is consistent but not among the most abundant.",
    },
    "OGN": {
        "slrp_class": "III",
        "literature_strength": "moderate bone/evolution evidence; weaker direct growth-plate evidence",
        "literature_note": "OGN regulates type-I collagen fibrillogenesis after proteolytic processing and has an informative vertebrate/teleost duplication history; rat growth-plate expression is weak.",
        "literature_url": "https://doi.org/10.1186/s12862-018-1310-2",
        "decision": "main, exploratory",
        "decision_note": "Keep for evolutionary/class-III scope and strong local P0 signal, but label it the least certain growth-plate member. The 14-protein vertebrate set is coherent after amphioxus exclusion; spotted gar remains tentative because SignalP is negative.",
    },
    "LUM": {
        "slrp_class": "II",
        "literature_strength": "strong cartilage/fibrillogenesis and growth-plate evidence",
        "literature_note": "Closely overlaps FMOD, is present in fetal growth plate/cartilage, and compensates in FMOD-null tissues.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/30700002/",
        "decision": "main",
        "decision_note": "Consistent rodent growth-plate expression and a complete, technically clean protein pipeline make LUM a strong class-II main candidate. Its completed 14-target SynVoy run provides locus-level candidates in all genomes, with prioritized manual confirmation still required because several genomes contain multiple or rescue-derived calls.",
    },
    "DCN": {
        "slrp_class": "I",
        "literature_strength": "strong established cartilage-matrix evidence",
        "literature_note": "Established cartilage pericellular-matrix and mechanotransduction role.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/33246102/",
        "decision": "main",
        "decision_note": "Promoted into the complete comparative package after the 2026-09-02 audit: reproducible juvenile and direct growth-plate expression, established cartilage-matrix biology, and a completed curated protein/tree/gene-structure package make DCN both a biologically relevant candidate and the essential class-I control for BGN orthology.",
    },
    "OMD": {
        "slrp_class": "II",
        "literature_strength": "bone/subchondral evidence; limited longitudinal growth-plate support",
        "literature_note": "Omd deficiency changes cortical bone width but not longitudinal bone growth.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/38722812/",
        "decision": "secondary/background; do not add to full pipeline now",
        "decision_note": "Local P0 and mouse support exist, but rat growth-plate values are absent and current evidence is more bone/subchondral than growth-plate-specific.",
    },
    "ASPN": {
        "slrp_class": "I",
        "literature_strength": "context-dependent cartilage/growth-center evidence",
        "literature_note": "Expressed in mandibular condylar cartilage; weaker support in tibial rodent growth-plate matrices.",
        "literature_url": "https://pubmed.ncbi.nlm.nih.gov/28875156/",
        "decision": "secondary/background",
        "decision_note": "Useful class-I context, but independent tibial growth-plate support is weaker than for DCN/BGN.",
    },
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    local_rows = read_tsv(
        ROOT
        / "analyses"
        / "expression"
        / "tables"
        / "mouse_expression_descriptive_metadata_corrected.tsv"
    )
    local = {
        row["gene_symbol"].upper(): row
        for row in local_rows
        if row["gene_symbol"].upper() in GENES
    }
    local_rank = {
        gene: rank
        for rank, gene in enumerate(
            sorted(
                GENES,
                key=lambda x: float(local[x]["p0_chondroprogenitor_mean_TPM"]),
                reverse=True,
            ),
            start=1,
        )
    }

    gse_rows = read_tsv(
        ROOT
        / "analyses"
        / "expression"
        / "gse114919"
        / "tables"
        / "gse114919_slrp_tibia_cross_condition_summary.tsv"
    )
    gse = {(x["species"], x["gene"]): x for x in gse_rows}

    bgee_rows = read_tsv(
        ROOT
        / "analyses"
        / "expression"
        / "tables"
        / "bgee_slrp_relevant_expression_summary.tsv"
    )
    bgee_mouse = {
        x["thesis_gene"]: x
        for x in bgee_rows
        if x["species"] == "mouse" and x["query_mode"] == "anatomy_all_data"
    }

    conservation_rows = read_tsv(
        ROOT
        / "analyses"
        / "protein_analysis"
        / "tables"
        / "protein_conservation_domain_msa_signalp_gene_summary.tsv"
    )
    conservation = {x["Gene"]: x for x in conservation_rows}
    synvoy_rows = read_tsv(
        ROOT / "analyses" / "synteny" / "tables" / "synvoy_run_summary.tsv"
    )
    synvoy = {x["gene"]: x for x in synvoy_rows}

    rows: list[dict] = []
    for gene in GENES:
        mouse = gse[("mouse", gene)]
        rat = gse[("rat", gene)]
        bgee = bgee_mouse[gene]
        protein = conservation.get(gene)
        synteny = synvoy.get(gene)
        annotation = ANNOTATION[gene]
        rows.append(
            {
                "gene": gene,
                "slrp_class": annotation["slrp_class"],
                "local_p0_mean_tpm": round(
                    float(local[gene]["p0_chondroprogenitor_mean_TPM"]), 3
                ),
                "local_p0_cv": round(float(local[gene]["p0_chondroprogenitor_CV"]), 3),
                "local_p0_rank_of_9": local_rank[gene],
                "mouse_growth_plate_mean_rank": mouse["mean_within_condition_rank"],
                "mouse_growth_plate_all_replicates_conditions": mouse[
                    "conditions_detected_in_all_replicates"
                ],
                "rat_growth_plate_mean_rank": rat["mean_within_condition_rank"],
                "rat_growth_plate_all_replicates_conditions": rat[
                    "conditions_detected_in_all_replicates"
                ],
                "bgee_mouse_direct_cartilage_call_count": bgee[
                    "direct_cartilage_or_growth_plate_call_count"
                ],
                "bgee_mouse_best_cartilage_score": bgee[
                    "best_priority_expression_score"
                ],
                "literature_strength": annotation["literature_strength"],
                "literature_note": annotation["literature_note"],
                "literature_url": annotation["literature_url"],
                "canonical_protein_panel_complete": "yes" if protein else "no",
                "canonical_protein_count": protein["Protein count"] if protein else "",
                "mean_forward_identity_percent": protein[
                    "Mean forward identity percent"
                ]
                if protein
                else "",
                "protein_needs_inspection_count": protein["Needs-inspection count"]
                if protein
                else "",
                "synvoy_14_species_complete": "yes" if synteny else "no",
                "synvoy_high_calls": synteny["final_high_goi"] if synteny else "",
                "synvoy_medium_calls": synteny["final_medium_goi"] if synteny else "",
                "recommended_role": annotation["decision"],
                "decision_note": annotation["decision_note"],
            }
        )

    fields = list(rows[0])
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT} ({len(rows)} genes)")


if __name__ == "__main__":
    main()
