# Analysis and result overview

The repository is arranged by analysis step
## Final evidence files

| File | Content |
| --- | --- |
| `tables/final_gene_level_evidence.tsv` | One summary row for each of the seven comparative genes and the OMD context gene. |
| `tables/final_candidate_level_confidence.tsv` | Candidate-level evidence, final decisions and reasons for exclusions or cautions. |
| `tables/results_interpretation_summary.tsv` | Short summary of what each analysis supports and what it cannot prove. |
| `tables/species_lineage_evidence_summary.tsv` | Coverage and evidence status across the species panel. |
| `verify_thesis_analysis_outputs.py` | Checks expected row counts, identifiers, key decisions, figures and script syntax. |

## Analysis map

| Analysis | Main scripts | Main result files |
| --- | --- | --- |
| Expression | `expression/*.py`, `expression/gse114919/*.py` | corrected mouse-expression tables, full-SLRP sensitivity tables and expression figures |
| Protein candidates | `protein_analysis/scripts/` | canonical 103-protein manifest, final FASTA files, conservation and manual-review tables |
| Domains and SignalP | `protein_analysis/scripts/parse_domain_scan.py`, `protein_analysis/signal_peptides/scripts/` | Pfam summary and final SignalP summary |
| Protein motifs and ProtSpace | `protein_analysis/motifs/scripts/`, `protein_analysis/protspace/scripts/` | motif assignments, canonical projection tables and figures |
| Phylogeny | `phylogenetics/compact_panel/scripts/`, `phylogenetics/combined_trees/` | compact per-gene trees, selected-gene tree and family-reference tree |
| Synteny | `synteny/scripts/` | reviewed SynVoy table, confidence table and focused locus diagnostics |
| Gene structure | `gene_structure/scripts/`, `overview/scripts/select_gene_structure_representatives.py` | orthology-supported structure, splice and domain-exon tables |
| Coding constraint | `evolutionary_rates/scripts/analyze_codon_constraint.py` | pairwise NG86 estimates and gene summary |
| Promoter motifs | `regulatory/scripts/` | promoter QC, motif summaries and exploratory figures |
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
    
