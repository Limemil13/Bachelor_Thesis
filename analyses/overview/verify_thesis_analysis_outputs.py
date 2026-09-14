#!/usr/bin/env python3
"""Compact integrity checks for the audited thesis analysis outputs."""

from __future__ import annotations

import ast
import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read_tsv(
    relative: str, *, allow_empty: bool = False
) -> list[dict[str, str]]:
    path = ROOT / relative
    assert path.is_file(), f"Missing: {path}"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        assert allow_empty, f"Empty: {path}"
        return []
    width = len(rows[0])
    assert all(len(row) == width for row in rows), f"Inconsistent TSV width: {path}"
    return rows


def main() -> None:
    signalp = read_tsv(
        "analyses/protein_analysis/signal_peptides/tables/signalp_summary_clean.tsv"
    )
    assert len(signalp) == 103
    assert len({x["signalp_id"] for x in signalp}) == 103
    assert Counter(x["signal_peptide"] for x in signalp) == {"yes": 100, "no": 3}
    assert {x["signalp_job_id"] for x in signalp} == {
        "6A89C3E50004F8DA087B239F",
        "6A8AEF7700077B786803DE7F",
        "6AA7E4B90014723A58E1A86B",
    }
    assert {x["signalp_model_mode"] for x in signalp} == {"slow-sequential"}

    integrated = read_tsv(
        "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv"
    )
    assert len(integrated) == 103
    assert len({(x["Gene"], x["Species"]) for x in integrated}) == 103
    assert Counter(x["Protein QC status"] for x in integrated) == {
        "Supported": 84,
        "Watch": 13,
        "Needs inspection": 6,
    }
    assert sum(x["SignalP prediction"] == "pending" for x in integrated) == 0

    manifest = read_tsv(
        "analyses/protein_analysis/candidates/canonical_candidate_manifest.tsv"
    )
    assert len(manifest) == 103
    assert Counter(x["gene"] for x in manifest) == {
        "BGN": 13,
        "DCN": 15,
        "EPYC": 15,
        "FMOD": 14,
        "LUM": 16,
        "OGN": 14,
        "PRELP": 16,
    }
    assert all(x.get("accession") != "XP_414298.2" for x in manifest)
    assert all(x.get("accession") != "XP_066265713.1" for x in manifest)

    domain_hits = read_tsv(
        "analyses/protein_analysis/domains/domain_scan_canonical/domain_hits_canonical.tsv"
    )
    domain_summary = read_tsv(
        "analyses/protein_analysis/domains/domain_scan_canonical/domain_summary_by_protein_canonical.tsv"
    )
    assert len(domain_hits) == 693
    assert len(domain_summary) == 103

    protspace = read_tsv(
        "analyses/protein_analysis/protspace/tables/protspace_canonical_embedding_qc.tsv"
    )
    # Current projections use 103 canonical SLRPs plus four controls; the QC
    # table intentionally reports only the 103 canonical proteins.
    assert len(protspace) == 103
    assert len({x["identifier"] for x in protspace}) == 103
    assert Counter(x["embedding_qc_status"] for x in protspace) == {
        "Supported": 85,
        "Watch": 1,
        "Needs inspection": 17,
    }

    bgee_mapping = read_tsv("analyses/expression/tables/bgee_slrp_gene_mapping.tsv")
    bgee_summary = read_tsv(
        "analyses/expression/tables/bgee_slrp_relevant_expression_summary.tsv"
    )
    assert len(bgee_mapping) == 27
    assert len(bgee_summary) == 54

    gse_long = read_tsv(
        "analyses/expression/gse114919/tables/gse114919_slrp_growth_plate_long.tsv"
    )
    gse_conditions = read_tsv(
        "analyses/expression/gse114919/tables/gse114919_slrp_condition_summary.tsv"
    )
    gse_tibia = read_tsv(
        "analyses/expression/gse114919/tables/gse114919_slrp_tibia_cross_condition_summary.tsv"
    )
    assert len(gse_long) == 531
    assert Counter(x["species"] for x in gse_long) == {"mouse": 261, "rat": 270}
    assert len(gse_conditions) == 108
    assert len(gse_tibia) == 18

    synvoy_all = read_tsv("analyses/synteny/tables/synvoy_gene_species_evidence.tsv")
    synvoy_p1 = read_tsv(
        "analyses/synteny/tables/synvoy_manual_review_queue.tsv", allow_empty=True
    )
    synvoy_p2 = read_tsv(
        "analyses/synteny/tables/synvoy_batch_confirmation_queue.tsv",
        allow_empty=True,
    )
    synvoy_p3 = read_tsv(
        "analyses/synteny/tables/synvoy_spot_check_queue.tsv", allow_empty=True
    )
    synvoy_p1_done = read_tsv(
        "analyses/synteny/tables/synvoy_manual_review_completed.tsv"
    )
    synvoy_p2_done = read_tsv(
        "analyses/synteny/tables/synvoy_batch_confirmation_completed.tsv"
    )
    synvoy_p3_done = read_tsv(
        "analyses/synteny/tables/synvoy_spot_check_completed.tsv"
    )
    synvoy_genes = {x["gene"] for x in synvoy_all}
    assert synvoy_genes == {"BGN", "DCN", "EPYC", "FMOD", "LUM", "OGN", "PRELP"}
    assert len(synvoy_all) == 98
    synvoy_priorities = Counter(x["review_priority"] for x in synvoy_all)
    assert not synvoy_p1 and not synvoy_p2 and not synvoy_p3
    assert len(synvoy_p1_done) == synvoy_priorities["P1"] == 25
    assert len(synvoy_p2_done) == synvoy_priorities["P2"] == 69
    assert len(synvoy_p3_done) == synvoy_priorities["P3"] == 4
    assert all(x["review_status"] for x in synvoy_all)
    assert Counter(x["review_status"] for x in synvoy_all) == {
        "accepted": 62,
        "tentative": 20,
        "ambiguous": 12,
        "rejected": 4,
    }
    assert Counter(x["input_pair_qc"] for x in synvoy_all) == {
        "PASS": 91,
        "FAIL": 7,
    }
    assert all(
        x["input_pair_qc"] == "FAIL" and x["review_status"] == "ambiguous"
        for x in synvoy_all
        if x["species"] == "opossum"
    )
    synvoy_runs = {
        x["gene"]: x for x in read_tsv("analyses/synteny/tables/synvoy_run_summary.tsv")
    }
    assert set(synvoy_runs) == synvoy_genes
    assert synvoy_runs["LUM"]["genomes_qc_pass"] == "14"
    assert synvoy_runs["LUM"]["final_high_goi"] == "10"
    assert synvoy_runs["LUM"]["final_medium_goi"] == "39"
    assert synvoy_runs["DCN"]["genomes_qc_pass"] == "14"
    assert synvoy_runs["DCN"]["final_high_goi"] == "10"
    assert synvoy_runs["DCN"]["final_medium_goi"] == "1"
    assert synvoy_runs["BGN"]["run"] == "bgn_human_15species_dev_20260905"
    assert synvoy_runs["BGN"]["final_high_goi"] == "10"
    assert synvoy_runs["BGN"]["final_medium_goi"] == "1"

    candidates = read_tsv("analyses/overview/candidate_gene_evidence.tsv")
    assert len(candidates) == 9
    assert {x["gene"] for x in candidates} == {
        "BGN",
        "FMOD",
        "PRELP",
        "EPYC",
        "OGN",
        "LUM",
        "DCN",
        "OMD",
        "ASPN",
    }

    seven_gene = read_tsv("analyses/overview/six_gene_integrated_summary.tsv")
    assert len(seven_gene) == 7
    assert {x["gene"] for x in seven_gene} == {
        "BGN",
        "DCN",
        "FMOD",
        "PRELP",
        "EPYC",
        "OGN",
        "LUM",
    }

    tree_review = read_tsv(
        "analyses/phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv"
    )
    assert len(tree_review) == 7
    assert Counter(x["unrooted_monophyletic"] for x in tree_review) == {
        "yes": 5,
        "no": 2,
    }

    combined_tree = ROOT / (
        "analyses/phylogenetics/combined_trees/six_gene_tree/"
        "SLRP_selected_working_genes.treefile"
    )
    combined_text = combined_tree.read_text(encoding="utf-8")
    combined_tips = re.findall(r"(?<=[(,])([^(),:;]+):", combined_text)
    assert len(combined_tips) == 103
    assert "XP_066265713.1" not in combined_text

    lum_small_tree = (
        ROOT / "analyses/phylogenetics/compact_panel/trees/LUM/LUM.treefile"
    )
    lum_small_tips = re.findall(
        r"(?<=[(,])([^(),:;]+):", lum_small_tree.read_text(encoding="utf-8")
    )
    assert len(lum_small_tips) == 16
    assert (ROOT / "analyses/phylogenetics/compact_panel/itol/LUM_itol.tree").is_file()

    bgn_small_tree = (
        ROOT / "analyses/phylogenetics/compact_panel/trees/BGN/BGN.treefile"
    )
    bgn_small_tips = re.findall(
        r"(?<=[(,])([^(),:;]+):", bgn_small_tree.read_text(encoding="utf-8")
    )
    assert len(bgn_small_tips) == 13
    assert "XP_414298.2" not in bgn_small_tree.read_text(encoding="utf-8")

    fmod_small_tree = (
        ROOT / "analyses/phylogenetics/compact_panel/trees/FMOD/FMOD.treefile"
    )
    fmod_small_text = fmod_small_tree.read_text(encoding="utf-8")
    fmod_small_tips = re.findall(r"(?<=[(,])([^(),:;]+):", fmod_small_text)
    assert len(fmod_small_tips) == 14
    assert "XP_075930353.1" not in fmod_small_text

    ogn_small_tree = (
        ROOT / "analyses/phylogenetics/compact_panel/trees/OGN/OGN.treefile"
    )
    ogn_small_text = ogn_small_tree.read_text(encoding="utf-8")
    ogn_small_tips = re.findall(r"(?<=[(,])([^(),:;]+):", ogn_small_text)
    assert len(ogn_small_tips) == 14
    assert "XP_066265713.1" not in ogn_small_text

    structure = read_tsv(
        "analyses/gene_structure/tables/gene_structure_representative_5species.tsv"
    )
    assert len(structure) == 40
    assert Counter(x["query_gene"] for x in structure) == {
        "BGN": 5,
        "DCN": 5,
        "EPYC": 5,
        "FMOD": 5,
        "LUM": 5,
        "OGN": 5,
        "OMD": 5,
        "PRELP": 5,
    }
    structure_summary = read_tsv(
        "analyses/gene_structure/tables/gene_structure_summary_5species.tsv"
    )
    assert len(structure_summary) == 8
    assert {x["query_gene"] for x in structure_summary} == {
        "BGN",
        "DCN",
        "EPYC",
        "FMOD",
        "LUM",
        "OGN",
        "OMD",
        "PRELP",
    }
    structure_extended = read_tsv(
        "analyses/gene_structure/tables/gene_structure_extended_summary.tsv"
    )
    assert len(structure_extended) == 8
    assert all(float(x["protein_length_cv_percent"]) < 2.0 for x in structure_extended)
    structure_outliers = (
        ROOT
        / "analyses/gene_structure/tables/gene_structure_structural_outliers_5species.tsv"
    )
    assert structure_outliers.is_file()
    assert len(structure_outliers.read_text(encoding="utf-8").splitlines()) == 1

    full_structure_qc = read_tsv(
        "analyses/gene_structure/extended/tables/full_gene_structure_qc.tsv"
    )
    assert len(full_structure_qc) == 40
    assert all(x["cds_length_divisible_by_three"] == "yes" for x in full_structure_qc)
    assert all(
        x["translation_starts_with_methionine"] == "yes" for x in full_structure_qc
    )
    assert all(x["terminal_stop_present"] == "yes" for x in full_structure_qc)
    assert sum(int(x["internal_stop_count"]) for x in full_structure_qc) == 0
    assert all(
        x["all_splice_phases_consistent_with_gff"] == "yes" for x in full_structure_qc
    )
    splice_junctions = read_tsv(
        "analyses/gene_structure/extended/tables/splice_junction_conservation_summary.tsv"
    )
    assert len(splice_junctions) == 26
    assert all(x["phase_conserved_across_species"] == "yes" for x in splice_junctions)
    assert all(x["all_gff_phase_checks_pass"] == "yes" for x in splice_junctions)
    isoform_sensitivity = read_tsv(
        "analyses/gene_structure/extended/tables/isoform_sensitivity_5species.tsv"
    )
    assert len(isoform_sensitivity) == 40
    assert Counter(x["alternative_coding_structure"] for x in isoform_sensitivity) == {
        "no": 29,
        "yes": 11,
    }
    representative_domain_hits = read_tsv(
        "analyses/gene_structure/extended/tables/representative_pfam_domain_hits.tsv"
    )
    representative_domain_summary = read_tsv(
        "analyses/gene_structure/extended/tables/representative_domain_exon_gene_summary.tsv"
    )
    assert len(representative_domain_hits) == 294
    assert len(representative_domain_summary) == 8
    assert all(
        x["all_species_have_lrr_support"] == "yes"
        for x in representative_domain_summary
    )

    pairwise_dn_ds = read_tsv(
        "analyses/evolutionary_rates/tables/pairwise_dn_ds_human_reference.tsv"
    )
    codon_summary = read_tsv(
        "analyses/evolutionary_rates/tables/codon_constraint_gene_summary.tsv"
    )
    assert len(pairwise_dn_ds) == 32 and len(codon_summary) == 8
    assert Counter(x["estimate_status"] for x in pairwise_dn_ds) == {
        "interpretable": 22,
        "caution": 3,
        "not interpretable": 6,
        "excluded orthology conflict": 1,
    }
    finite_omega = [float(x["dN_dS_omega"]) for x in pairwise_dn_ds if x["dN_dS_omega"]]
    assert len(finite_omega) == 25 and all(value < 1 for value in finite_omega)
    assert all(x["all_interpretable_omega_below_one"] == "yes" for x in codon_summary)

    corrected_large_counts = {
        "BGN": 639,
        "EPYC": 632,
        "FMOD": 546,
        "LUM": 838,
        "OGN": 855,
        "PRELP": 848,
    }
    for gene, expected in corrected_large_counts.items():
        report = read_tsv(
            f"analyses/phylogenetics/per_gene_trees/{gene}/candidates/"
            f"{gene}_large_locus_filter_report.tsv"
        )
        assert Counter(x["status"] for x in report)["kept"] == expected

    structure_queue = read_tsv(
        "analyses/overview/gene_structure_manual_review_queue.tsv"
    )
    assert len(structure_queue) == 35
    assert Counter(x["priority"] for x in structure_queue) == {
        "P1": 1,
        "P2": 5,
        "P3": 29,
    }
    predicted_transcripts = read_tsv(
        "analyses/overview/gene_structure_predicted_transcript_ncbi_review.tsv"
    )
    assert len(predicted_transcripts) == 6
    assert Counter(x["review_status"] for x in predicted_transcripts) == {
        "database review complete; retain predicted model": 4,
        "database review complete; retain with teleost-paralog wording": 1,
        "unresolved annotation conflict": 1,
    }

    manual_review = read_tsv(
        "analyses/protein_analysis/manual_review/manual_sequence_review_5.tsv"
    )
    assert {x["Protein accession"] for x in manual_review} == {
        "XP_048463897.1",
        "XP_056660002.1",
        "NP_001013588.1",
        "XP_006631165.2",
        "XP_066265713.1",
    }
    manual_fasta = (
        ROOT / "analyses/protein_analysis/manual_review/manual_sequence_review_5.faa"
    )
    assert (
        sum(line.startswith(">") for line in manual_fasta.read_text().splitlines()) == 5
    )
    manual_decisions = read_tsv(
        "analyses/protein_analysis/manual_review/manual_sequence_review_decisions.tsv"
    )
    assert len(manual_decisions) == 9
    assert {x["accession"] for x in manual_decisions} == {
        "XP_048463897.1",
        "XP_056660002.1",
        "NP_001013588.1",
        "XP_006631165.2",
        "XP_066265713.1",
        "XP_038638277.1",
        "XP_048475905.1",
        "XP_001363160.3",
        "XP_038636759.1",
    }
    decisions_by_accession = {
        x["accession"]: x["manual_decision"] for x in manual_decisions
    }
    assert decisions_by_accession["XP_048463897.1"] == "retain as tentative"
    assert decisions_by_accession["XP_056660002.1"] == "retain"
    assert decisions_by_accession["NP_001013588.1"] == "retain"
    assert decisions_by_accession["XP_006631165.2"] == "retain as tentative"
    assert decisions_by_accession["XP_066265713.1"].startswith("exclude")
    assert decisions_by_accession["XP_038638277.1"].startswith("retain as tentative")
    assert decisions_by_accession["XP_048475905.1"].startswith("retain as tentative")
    assert decisions_by_accession["XP_001363160.3"].startswith("retain as partial")
    assert decisions_by_accession["XP_038636759.1"].startswith("retain as confident")
    manual_synteny = read_tsv(
        "analyses/protein_analysis/manual_review/"
        "manual_sequence_review_5_synteny_context.tsv"
    )
    assert len(manual_synteny) == 5
    assert {x["canonical_protein_accession"] for x in manual_synteny} == {
        "XP_048463897.1",
        "XP_056660002.1",
        "NP_001013588.1",
        "XP_006631165.2",
        "XP_066265713.1",
    }
    assert Counter(x["review_priority"] for x in manual_synteny) == {
        "P1": 1,
        "P2": 4,
    }

    protein_motif_assignments = read_tsv(
        "analyses/protein_analysis/motifs/tables/"
        "protein_candidate_motif_model_assignments.tsv"
    )
    assert len(protein_motif_assignments) == 103
    assert Counter(x["motif_model_status"] for x in protein_motif_assignments) == {
        "own gene model best": 103
    }
    protein_meme = read_tsv(
        "analyses/protein_analysis/motifs/tables/protein_meme_motif_summary.tsv"
    )
    assert len(protein_meme) == 116
    assert Counter(x["scope"] for x in protein_meme)["all_seven_genes"] == 37

    promoter_qc = read_tsv("analyses/regulatory/tables/promoter_sequence_qc.tsv")
    assert len(promoter_qc) == 70
    assert Counter(x["window"] for x in promoter_qc) == {"core": 35, "proximal": 35}
    assert {x["clipped_at_contig_edge"] for x in promoter_qc} == {"no"}
    assert all(float(x["ambiguous_percent"]) == 0 for x in promoter_qc)
    promoter_denovo = read_tsv(
        "analyses/regulatory/tables/promoter_denovo_motif_summary.tsv"
    )
    assert len(promoter_denovo) == 40
    promoter_strict = read_tsv(
        "analyses/regulatory/tables/promoter_targeted_tf_strict_scan_summary.tsv"
    )
    assert len(promoter_strict) == 1
    assert promoter_strict[0]["reported_hits"] == "0"
    promoter_ame = read_tsv(
        "analyses/regulatory/tables/promoter_ame_targeted_tf_summary.tsv"
    )
    assert len(promoter_ame) == 40
    assert {x["panel_significant"] for x in promoter_ame} == {"no"}

    mgi = read_tsv("analyses/phenotypes/tables/mgi_phenotype_evidence_summary.tsv")
    hpo = read_tsv(
        "analyses/phenotypes/tables/hpo_human_phenotype_evidence_summary.tsv"
    )
    assert len(mgi) == 7 and len(hpo) == 7
    assert {x["gene"] for x in mgi} == {
        "BGN",
        "DCN",
        "FMOD",
        "PRELP",
        "EPYC",
        "LUM",
        "OGN",
    }
    assert sum(int(x["single_gene_unique_mp_terms"]) for x in mgi) == 105
    assert {x["gene"] for x in hpo if int(x["unique_hpo_terms"]) > 0} == {"BGN", "DCN"}

    final_gene = read_tsv("analyses/overview/tables/final_gene_level_evidence.tsv")
    assert len(final_gene) == 8
    assert {x["Gene"] for x in final_gene} == {
        "BGN",
        "DCN",
        "FMOD",
        "PRELP",
        "EPYC",
        "LUM",
        "OGN",
        "OMD",
    }
    final_candidate = read_tsv(
        "analyses/overview/tables/final_candidate_level_confidence.tsv"
    )
    assert len(final_candidate) == 107
    assert {"NP_001025243.1", "XP_007897806.2"}.issubset(
        {x["accession"] for x in final_candidate}
    )
    assert sum(x["final decision"].startswith("exclude") for x in final_candidate) == 2
    assert sum(x["final decision"].startswith("pending") for x in final_candidate) == 0
    elephant_fmod = next(
        x for x in final_candidate if x["accession"] == "XP_007897806.2"
    )
    assert elephant_fmod["protein length"].startswith("684")
    assert "compound" in elephant_fmod["protein length"]
    assert "exclude the compound protein model" in elephant_fmod["final decision"]

    elephant_qc = read_tsv(
        "analyses/synteny/diagnostics/priority_locus_reviews/"
        "FMOD_elephant_shark/XP_007897806.2_qc.tsv"
    )
    assert len(elephant_qc) == 1
    assert elephant_qc[0]["protein_length_aa"] == "684"
    assert elephant_qc[0]["internal_stop_count"] == "0"
    assert elephant_qc[0]["terminal_stop_present"] == "yes"
    elephant_diagnostic = read_tsv(
        "analyses/synteny/diagnostics/priority_locus_reviews/"
        "FMOD_elephant_shark/XP_007897806.2_diagnostic_summary.tsv"
    )
    assert len(elephant_diagnostic) == 1
    assert elephant_diagnostic[0]["n_terminal_best_reference"] == "DCN"
    assert elephant_diagnostic[0]["c_terminal_best_reference"] == "FMOD"
    assert "54-76" in elephant_diagnostic[0]["pfam_architecture"]
    assert "383-412" in elephant_diagnostic[0]["pfam_architecture"]

    species_summary = read_tsv(
        "analyses/overview/tables/species_lineage_evidence_summary.tsv"
    )
    assert len(species_summary) == 17
    species_by_name = {row["display_species"]: row for row in species_summary}
    assert species_by_name["Sea lamprey"]["canonical_proteins"] == "3"
    assert species_by_name["Spotted gar"]["synteny_accepted"] == "6"
    assert species_by_name["Coelacanth"]["synteny_accepted"] == "6"
    assert species_by_name["Opossum†"]["synvoy_input_qc"] == "FAIL"
    assert species_by_name["Amphioxus‡"]["canonical_proteins"] == "0"

    for document in (
        "analyses/README.md",
        "analyses/overview/AUTHORITATIVE_RESULT_LOCATIONS.md",
        "analyses/overview/RESULTS_INTERPRETATION_GUIDE.md",
        "analyses/overview/WORKFLOW_OVERVIEW.md",
        "analyses/overview/THESIS_CONCLUSION.md",
        "analyses/overview/METHODS_ALGORITHM_NOTES.md",
        "analyses/overview/SIX_GENE_PROFILES.md",
        "analyses/overview/FIGURE_CAPTIONS.md",
        "analyses/overview/THESIS_NEXT_STEPS.md",
        "analyses/overview/THESIS_SELF_WRITING_GUIDE.md",
        "analyses/overview/THESIS_STRUCTURE_AND_WRITING_PLAN.md",
        "analyses/overview/GITHUB_PUSH_MANIFEST.md",
        "analyses/overview/history/THESIS_ANALYSIS_LOG_2026-08-25.md",
        "analyses/overview/history/THESIS_ANALYSIS_LOG_2026-08-29.md",
        "analyses/overview/history/THESIS_ANALYSIS_LOG_2026-09-07.md",
        "analyses/overview/history/THESIS_ANALYSIS_LOG_2026-09-10.md",
        "analyses/overview/history/THESIS_ANALYSIS_LOG_2026-09-14.md",
        "analyses/literature/README.md",
        "thesis/chapters/results.tex",
        "thesis/chapters/discussion.tex",
        "thesis/bibliography.bib",
        "analyses/synteny/SYNTENY_ANALYSIS_STATUS.md",
        "analyses/synteny/diagnostics/priority_locus_reviews/"
        "FMOD_elephant_shark/README.md",
        "analyses/gene_structure/README.md",
        "analyses/evolutionary_rates/README.md",
        "analyses/regulatory/README.md",
        "analyses/protein_analysis/motifs/README.md",
        "analyses/phenotypes/README.md",
        "thesis/chapters/introduction.tex",
        "thesis/chapters/methods.tex",
        "thesis/chapters/conclusion.tex",
        "thesis/chapters/appendix.tex",
    ):
        assert (ROOT / document).is_file(), f"Missing: {document}"
    assert (
        ROOT / "analyses/expression/figures/six_gene_expression_evidence.png"
    ).is_file()
    assert (
        ROOT / "analyses/expression/figures/six_gene_expression_evidence.svg"
    ).is_file()
    for gene in ("BGN", "DCN", "EPYC", "FMOD", "LUM", "OGN", "PRELP"):
        for suffix in ("png", "svg"):
            assert (
                ROOT
                / f"analyses/protein_analysis/figures/sequence_logos/{gene}_sequence_logo.{suffix}"
            ).is_file()

    overview_figures = (
        "analyses/overview/figures/six_gene_evidence_dashboard",
        "analyses/overview/figures/thesis_workflow_overview",
        "analyses/overview/figures/species_lineage_evidence_summary",
        "analyses/protein_analysis/figures/protein_conservation_overview",
        "analyses/expression/figures/expression_age_and_bgee_overview",
        "analyses/synteny/figures/synvoy_species_confidence_heatmap",
        "analyses/phylogenetics/figures/phylogenetic_result_overview",
        "analyses/phylogenetics/figures/compact_per_gene_trees",
        "analyses/phylogenetics/figures/six_gene_combined_tree",
        "analyses/gene_structure/figures/gene_structure_conservation_overview",
        "analyses/gene_structure/extended/figures/full_gene_structure_5species",
        "analyses/gene_structure/extended/figures/coding_exon_splice_phase_map",
        "analyses/gene_structure/extended/figures/human_domain_exon_architecture",
        "analyses/evolutionary_rates/figures/pairwise_dn_ds_human_reference",
    )
    for stem in overview_figures:
        for suffix in ("png", "svg"):
            assert (ROOT / f"{stem}.{suffix}").is_file(), f"Missing {stem}.{suffix}"

    supplementary_figures = (
        "analyses/protein_analysis/motifs/figures/protein_gene_motif_model_matrix.png",
        "analyses/protein_analysis/motifs/figures/protein_candidate_motif_model_margins.png",
        "analyses/regulatory/figures/promoter_proximal_gc_heatmap.png",
        "analyses/regulatory/figures/promoter_targeted_tf_ame_enrichment.png",
        "analyses/regulatory/figures/promoter_targeted_tf_species_recurrence.png",
        "analyses/phenotypes/figures/mgi_single_gene_phenotype_heatmap.png",
    )
    for relative in supplementary_figures:
        assert (ROOT / relative).is_file(), f"Missing {relative}"

    lum_large_tree = (
        ROOT
        / "analyses/phylogenetics/per_gene_trees/LUM/final_tree/LUM_large_final.treefile"
    )
    lum_large_report = (
        ROOT
        / "analyses/phylogenetics/per_gene_trees/LUM/final_tree/LUM_large_final.iqtree"
    )
    # IQ-TREE can write a temporary unique-sequence treefile during a live run;
    # the .iqtree report appears only when the expanded final result is complete.
    if lum_large_report.exists():
        assert lum_large_tree.is_file()
        text = lum_large_tree.read_text(encoding="utf-8")
        tips = re.findall(r"(?<=[(,])([^(),:;]+):", text)
        assert len(tips) == 781
        assert len(set(tips)) == 781
        lum_itol_base = ROOT / "analyses/phylogenetics/per_gene_trees/itol"
        for suffix in ("tree", "labels.txt", "taxgroups.txt", "taxonomy_mapping.tsv"):
            assert (lum_itol_base / f"LUM_large_final_itol.{suffix}").is_file()
        lum_taxonomy = read_tsv(
            "analyses/phylogenetics/per_gene_trees/itol/"
            "LUM_large_final_itol.taxonomy_mapping.tsv"
        )
        assert len(lum_taxonomy) == 781
        assert Counter(x["group"] for x in lum_taxonomy).get("Unknown", 0) == 0

    scripts = [
        "analyses/protein_analysis/signal_peptides/scripts/parse_signalp_table.py",
        "analyses/protein_analysis/protspace/scripts/analyze_protspace_canonical.py",
        "analyses/expression/fetch_bgee_slrp_expression.py",
        "analyses/expression/summarize_bgee_slrp_expression.py",
        "analyses/expression/gse114919/summarize_gse114919_growth_plate.py",
        "analyses/expression/plot_six_gene_expression_evidence.py",
        "analyses/synteny/scripts/build_synvoy_review.py",
        "analyses/overview/build_candidate_gene_evidence.py",
        "analyses/protein_analysis/scripts/build_manual_sequence_review_packet.py",
        "analyses/protein_analysis/scripts/make_sequence_logos.py",
        "analyses/phylogenetics/per_gene_trees/scripts/make_large_taxgroup_annotations.py",
        "analyses/phylogenetics/per_gene_trees/scripts/filter_large_slrp_fasta.py",
        "analyses/phylogenetics/compact_panel/scripts/plot_compact_phylogenies.py",
        "analyses/overview/scripts/build_results_synthesis.py",
        "analyses/overview/scripts/build_final_evidence_tables.py",
        "analyses/overview/scripts/build_species_lineage_summary.py",
        "analyses/protein_analysis/motifs/scripts/prepare_protein_motif_inputs.py",
        "analyses/protein_analysis/motifs/scripts/summarize_protein_motifs.py",
        "analyses/regulatory/scripts/extract_promoters.py",
        "analyses/regulatory/scripts/subset_meme_database.py",
        "analyses/regulatory/scripts/summarize_regulatory_motifs.py",
        "analyses/phenotypes/scripts/build_mgi_phenotype_summary.py",
        "analyses/phenotypes/scripts/summarize_hpo_focal_genes.py",
        "analyses/gene_structure/scripts/build_extended_gene_structure.py",
        "analyses/gene_structure/scripts/map_domains_to_exons.py",
        "analyses/evolutionary_rates/scripts/analyze_codon_constraint.py",
        "analyses/synteny/scripts/audit_synvoy_canonical_loci.py",
        "analyses/synteny/scripts/finalize_all_synvoy_reviews.py",
        "analyses/synteny/scripts/validate_synvoy_target_gffs.py",
    ]
    for relative in scripts:
        path = ROOT / relative
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    thesis_tex_files = [ROOT / "thesis/main.tex"]
    thesis_tex_files.extend(sorted((ROOT / "thesis/chapters").glob("*.tex")))
    include_pattern = re.compile(
        r"\\includegraphics(?:\[[^\]]*\])?\s*\{([^}]+)\}", re.MULTILINE
    )
    supplement_pattern = re.compile(
        r"\\thesissuppfigure\s*\{([^}]+)\}", re.MULTILINE
    )
    figure_paths: list[Path] = []
    labels: list[str] = []
    references: list[str] = []
    for tex_path in thesis_tex_files:
        tex = tex_path.read_text(encoding="utf-8")
        relative_figures = [
            relative for relative in include_pattern.findall(tex) if "#" not in relative
        ]
        relative_figures.extend(supplement_pattern.findall(tex))
        figure_paths.extend(ROOT / "thesis" / relative for relative in relative_figures)
        labels.extend(re.findall(r"\\label\{([^}]+)\}", tex))
        references.extend(re.findall(r"\\(?:ref|pageref)\{([^}]+)\}", tex))
        assert tex.count(r"\begin{figure}") == tex.count(r"\end{figure}"), (
            f"Unbalanced figure environment: {tex_path}"
        )
    assert figure_paths, "No LaTeX figure references found"
    for figure_path in figure_paths:
        assert figure_path.is_file(), f"Missing LaTeX figure: {figure_path}"
    duplicate_labels = [label for label, count in Counter(labels).items() if count > 1]
    assert not duplicate_labels, f"Duplicate LaTeX labels: {duplicate_labels}"
    unresolved_references = sorted(set(references) - set(labels))
    assert not unresolved_references, (
        f"Unresolved LaTeX references: {unresolved_references}"
    )

    print("PASS: canonical manifest 103 and Pfam 693 hits / 103 summaries")
    print(
        "PASS: SignalP 103 calls (100 positive, 3 negative) and no pending calls"
    )
    print("PASS: integrated protein table 103 and ProtSpace canonical-QC table 103")
    print("PASS: Bgee 27 mappings / 54 summaries")
    print("PASS: GSE114919 531 long rows / 108 conditions / 18 tibia summaries")
    print(
        f"PASS: SynVoy {len(synvoy_all)} rows across {len(synvoy_genes)} genes -> "
        + " / ".join(f"{synvoy_priorities[p]} {p}" for p in ("P1", "P2", "P3"))
    )
    print("PASS: all SynVoy reviews complete; seven opossum input-pair failures quarantined")
    print("PASS: candidate evidence 9 genes")
    print(
        "PASS: integrated seven-gene evidence, 103-tip review, and canonical compact trees"
    )
    if lum_large_report.exists():
        print("PASS: 781-tip large LUM tree and complete iTOL taxonomy bundle")
    print("PASS: local gene structure 40 rows / 8 summaries / no structural outliers")
    print(
        "PASS: extended structure 40 complete CDSs / 26 conserved splice phases / Pfam support in all 40"
    )
    print(
        "PASS: coding constraint 25 finite NG86 estimates, all below one; deep comparisons explicitly limited"
    )
    print(
        "PASS: corrected large-tree locus-filter reports and thesis-wide overview figures"
    )
    print(
        "PASS: gene-structure queue 35 / predicted-transcript audit 6 / five-protein manual packet and explicit decisions"
    )
    print("PASS: final 8-row gene evidence table and 107-row candidate audit")
    print("PASS: protein motifs 116 motifs / 103 own-gene-best assignments")
    print(
        "PASS: promoter motifs 70 QC windows / 40 de-novo motifs / no corrected targeted hit"
    )
    print("PASS: curated MGI/HPO summaries for all seven comparative genes")
    print("PASS: analysis scripts parse successfully")
    print(
        f"PASS: {len(figure_paths)} LaTeX figure paths and "
        f"{len(references)} cross-references resolve; labels are unique"
    )


if __name__ == "__main__":
    main()
