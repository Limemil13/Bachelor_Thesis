# Analysis and result overview

The repository is organized by analysis step. The files below provide the
shortest route from the thesis statements to the supporting results.

## Final evidence files

| File | Content |
| --- | --- |
| [`tables/final_gene_level_evidence.tsv`](tables/final_gene_level_evidence.tsv) | One summary row for each of the seven comparative genes and the OMD context gene. |
| [`tables/final_candidate_level_confidence.tsv`](tables/final_candidate_level_confidence.tsv) | Candidate-level evidence, final decisions and reasons for exclusions or cautions. |
| [`tables/results_interpretation_summary.tsv`](tables/results_interpretation_summary.tsv) | Short summary of what each analysis supports and what it cannot prove. |
| [`tables/species_lineage_evidence_summary.tsv`](tables/species_lineage_evidence_summary.tsv) | Coverage and evidence status across the species panel. |
| [`verify_thesis_analysis_outputs.py`](verify_thesis_analysis_outputs.py) | Checks expected row counts, identifiers, key decisions, figures and script syntax. |
| [`../synteny/tables/synvoy_gene_species_evidence.tsv`](../synteny/tables/synvoy_gene_species_evidence.tsv) | Gene-by-species SynVoy evidence used for the final synteny summary. |
| [`../synteny/diagnostics/priority_locus_reviews/manual_review_decisions.tsv`](../synteny/diagnostics/priority_locus_reviews/manual_review_decisions.tsv) | Manual decisions for loci that required additional review. |
| [`../gene_structure/extended/tables/structure_protein_crosswalk.tsv`](../gene_structure/extended/tables/structure_protein_crosswalk.tsv) | Links each gene-structure representative to the canonical protein set. |
| [`../protein_analysis/motifs/validation/results/tables/loso_summary.tsv`](../protein_analysis/motifs/validation/results/tables/loso_summary.tsv) | Independent leave-one-species-out protein-pattern validation summary. |
| [`../regulatory/cross_species_conservation/results/tables/promoter_alignment_summary.tsv`](../regulatory/cross_species_conservation/results/tables/promoter_alignment_summary.tsv) | Cross-species promoter-alignment summary for each gene and window. |

## Analysis map

| Analysis | Main scripts | Main result files |
| --- | --- | --- |
| Expression | `expression/*.py`, `expression/gse114919/*.py` | corrected mouse-expression tables, full-SLRP sensitivity tables and expression figures |
| Protein candidates | `protein_analysis/scripts/` | canonical 103-protein manifest, final FASTA files, conservation and manual-review tables |
| Domains and SignalP | `protein_analysis/scripts/parse_domain_scan.py`, `protein_analysis/signal_peptides/scripts/` | Pfam summary and final SignalP summary |
| Protein motifs and ProtSpace | `protein_analysis/motifs/scripts/`, `protein_analysis/motifs/validation/scripts/`, `protein_analysis/protspace/scripts/` | motif assignments, leave-one-species-out validation, canonical projection tables and figures |
| Phylogeny | `phylogenetics/compact_panel/scripts/`, `phylogenetics/combined_trees/`, `phylogenetics/targeted_extensions/scripts/` | compact per-gene trees, the selected-gene tree, and the expanded BGN--DCN and FMOD--LUM family trees |
| Synteny | `synteny/scripts/` | reviewed SynVoy table, confidence table and focused locus diagnostics |
| Gene structure | `gene_structure/scripts/`, `overview/scripts/select_gene_structure_representatives.py` | orthology-supported structure, splice, protein-crosswalk and domain-exon tables |
| Coding constraint | `evolutionary_rates/scripts/analyze_codon_constraint.py` | pairwise NG86 estimates and gene summary |
| Promoter motifs | `regulatory/scripts/`, `regulatory/cross_species_conservation/scripts/` | promoter QC, motif summaries, cross-species identity tables and exploratory figures |
| Phenotypes | `phenotypes/scripts/` | curated MGI and HPO summaries |

## Scope and interpretation

The final protein panel contains 103 proteins from seven SLRP genes. The
different analysis steps were used to check different parts of the candidate
assignments. For example, the protein sequence, phylogenetic tree, genomic
location and gene structure can each support an assignment.

Expression and phenotype data were mainly used to explain why the genes are
relevant to cartilage and the growth plate. They do not prove that an ortholog
has exactly the same growth-plate function in every species.

The original opossum SynVoy rows used a genome FASTA and GFF annotation from
different assembly versions and were therefore replaced. Seven matched-input
reruns were checked at the accession, locus and neighbouring-gene level. DCN,
EPYC, FMOD, LUM, OGN and PRELP were accepted. BGN remains unresolved because
the expected ATP2B3--HAUS7 interval is mostly ambiguous sequence, so this result
cannot be used to claim that BGN is absent. Other suspicious cases remain in
the candidate-level table, including the ASPN-like chicken BGN record, the
uncertain amphioxus OGN record, the partial deep-lineage BGN proteins and the
compound elephant-shark FMOD-like model.

The opossum BGN result was also checked against the older MonDom5 annotation
and Ensembl Compara. Neither source supplied a confident BGN replacement, and
the strongest old-interval SLRP control match was ECM2-like rather than
BGN-specific. Because the current expected interval is mostly unresolved and
published antibody evidence detected biglycan in opossum cartilage, the final
interpretation remains unresolved rather than absent.
    
