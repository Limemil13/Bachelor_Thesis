# Authoritative thesis result locations

Last audited: 2026-10-01

This is the current result map. Some historical directory names still contain
`six_gene`; those stable paths now hold the rebuilt seven-gene/103-protein
results and should not be interpreted from their names alone.

## Final integrated tables

- Gene-level evidence (seven comparative genes plus OMD context):
  `analyses/overview/tables/final_gene_level_evidence.tsv`
- Candidate-level confidence (103 canonical rows plus two exclusion-provenance rows):
  `analyses/overview/tables/final_candidate_level_confidence.tsv`
- Integrated seven-gene measurements (legacy filename):
  `analyses/overview/six_gene_integrated_summary.tsv`
- Candidate-selection evidence for nine screened genes:
  `analyses/overview/candidate_gene_evidence.tsv`

The protein conservation/domain/MSA/signal-peptide summary belongs **only** at:

`analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv`

Its one-row-per-gene aggregate is:

`analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_gene_summary.tsv`

Do not maintain a second editable copy in a thesis, ProtSpace, or tree folder.

## Candidate proteins

- Canonical 103-protein FASTA:
  `analyses/protein_analysis/candidates/slrp_candidates_canonical.faa`
- Canonical manifest:
  `analyses/protein_analysis/candidates/canonical_candidate_manifest.tsv`
- Per-gene canonical FASTAs:
  `analyses/protein_analysis/candidates/canonical_by_gene/`
- BGN/DCN correction audit:
  `analyses/protein_analysis/candidates/BGN_deep_lineage_correction.tsv`
- Manual decisions:
  `analyses/protein_analysis/manual_review/manual_sequence_review_decisions.tsv`

## Domains, MSA, conservation, and motifs

- Pfam hits and per-protein summaries:
  `analyses/protein_analysis/domains/domain_scan_canonical/`
- Canonical alignments:
  `analyses/protein_analysis/alignments/canonical/`
- Consensus, identity matrices, and conservation:
  `analyses/protein_analysis/alignments/canonical/conservation/`
- Sequence logos:
  `analyses/protein_analysis/figures/sequence_logos/`
- Protein motif tables and figures:
  `analyses/protein_analysis/motifs/tables/` and
  `analyses/protein_analysis/motifs/figures/`

## Signal peptides

- Current merged calls for all 103 proteins:
  `analyses/protein_analysis/signal_peptides/tables/signalp_summary_clean.tsv`
- Completed final 17-protein input:
  `analyses/protein_analysis/signal_peptides/input/pending_17_canonical_for_signalp.faa`
- Final 17-protein ID mapping:
  `analyses/protein_analysis/signal_peptides/input/pending_17_signalp_id_mapping.tsv`
- Raw final-job bundle:
  `analyses/protein_analysis/signal_peptides/raw_outputs/pending_17_slow_20260914_job_6AA7E4B90014723A58E1A86B/`
- Integrated 103-row table:
  `analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv`

## Phylogenies and iTOL

- Curated per-gene trees:
  `analyses/phylogenetics/compact_panel/trees/`
- Combined 103-tip tree (legacy directory name):
  `analyses/phylogenetics/combined_trees/six_gene_tree/SLRP_selected_working_genes.treefile`
- Automated clade review:
  `analyses/phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv`
- BGN/DCN discrimination diagnostic:
  `analyses/phylogenetics/compact_panel/tables/class_I_BGN_DCN_tree_discrimination.tsv`
- iTOL bundles:
  `analyses/phylogenetics/compact_panel/itol/`
- Thesis tree figures:
  `analyses/phylogenetics/figures/`
- Corrected large-tree input/filter reports:
  `analyses/phylogenetics/per_gene_trees/<GENE>/candidates/`

The compact trees are the central thesis trees. Corrected large-tree reruns are
optional unless those expanded trees will be presented as final evidence.

## ProtSpace

- Canonical-with-controls input:
  `analyses/protein_analysis/protspace/input/slrp_embedding_canonical_with_controls.faa`
- Current dated output bundle:
  `analyses/protein_analysis/protspace/output/slrp_protspace_seven_gene_curated_with_controls_20260903/`
- Canonical QC tables:
  `analyses/protein_analysis/protspace/tables/`
- Projection figure:
  `analyses/protein_analysis/protspace/figures/protspace_canonical_projections.png`

## Gene structure and coding constraint

- Five-species representative table:
  `analyses/gene_structure/tables/gene_structure_representative_5species.tsv`
- Current eight-gene structure summary:
  `analyses/gene_structure/tables/gene_structure_extended_summary.tsv`
- Full transcript/CDS/splice/domain outputs:
  `analyses/gene_structure/extended/`
- Gene-structure figures:
  `analyses/gene_structure/figures/` and
  `analyses/gene_structure/extended/figures/`
- Human-reference NG86 tables and figure:
  `analyses/evolutionary_rates/tables/` and
  `analyses/evolutionary_rates/figures/`

## Expression, promoter, and phenotype context

- Corrected GSE305415 WT P0 expression summary (from stored Salmon files):
  `analyses/expression/tables/mouse_expression_descriptive_metadata_corrected.tsv`
- Mouse/rat microdissected growth-plate results:
  `analyses/expression/gse114919/tables/`
- Bgee summaries:
  `analyses/expression/tables/bgee_slrp_relevant_expression_summary.tsv`
- Expression figures:
  `analyses/expression/figures/`
- Seven-gene promoter tables and figures:
  `analyses/regulatory/tables/` and `analyses/regulatory/figures/`
- Seven-gene MGI/HPO tables and phenotype heatmap:
  `analyses/phenotypes/tables/` and `analyses/phenotypes/figures/`

## Synteny

- Current completed seven-gene automated summary:
  `analyses/synteny/tables/synvoy_run_summary.tsv`
- One-row-per-gene/species evidence table:
  `analyses/synteny/tables/synvoy_gene_species_evidence.tsv`
- Complete 98-row review decisions:
  `analyses/synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`
- Completed P1/P2/P3 views:
  `analyses/synteny/tables/synvoy_manual_review_completed.tsv`,
  `synvoy_batch_confirmation_completed.tsv`, and
  `synvoy_spot_check_completed.tsv`
- Manual-review procedure:
  `analyses/synteny/SYNTENY_MANUAL_REVIEW_GUIDE.md`
- Elephant-shark compound FMOD-like model diagnostic:
  `analyses/synteny/diagnostics/priority_locus_reviews/FMOD_elephant_shark/`

## Species-lineage synthesis

- Per-species quantitative table:
  `analyses/overview/tables/species_lineage_evidence_summary.tsv`
- Thesis-ready figure:
  `analyses/overview/figures/species_lineage_evidence_summary.png`
- Evolutionary interpretation and literature context:
  `analyses/overview/SLRP_EVOLUTIONARY_ORIGIN_AND_SPECIES.md`

DCN completed at `results/dcn_human_15species_dev_20260903`, BGN at
`results/bgn_human_15species_dev_20260905`, and FMOD at
`results/fmod_human_15species_dev_20260906`; all three updated reports are
included in the authoritative 98-row table. All current locus reviews are
complete. A matched opossum FASTA/GFF rerun is required; updated PRELP and EPYC
runs are optional sensitivity analyses. The staged WSL procedure is in
`analyses/synteny/SYNVOY_UPDATED_RERUN_RUNBOOK.md`.

## Interpretation and writing

- Detailed gene biology: `analyses/overview/GENE_FUNCTION_REFERENCE.md`
- Gene-to-thesis profiles: `analyses/overview/SIX_GENE_PROFILES.md`
- Evolutionary origin/species context:
  `analyses/overview/SLRP_EVOLUTIONARY_ORIGIN_AND_SPECIES.md`
- Figure captions: `analyses/overview/FIGURE_CAPTIONS.md`
- Reproducibility status: `analyses/overview/THESIS_REPRODUCIBILITY_STATUS.md`
- Thesis drafts: `thesis/chapters/`
- Bibliography: `thesis/bibliography.bib`

Personal writing aids, dated work logs, and the current task checklist are kept
locally and excluded from the public repository.
