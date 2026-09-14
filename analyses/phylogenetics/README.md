# Phylogenetic analyses

All thesis tree work is grouped here and separated by biological purpose.

- `compact_panel/`: one small cross-species tree per gene using the canonical
  103-protein panel. `trees/` contains the seven final per-gene trees and `itol/`
  contains their matching iTOL annotations.
- `combined_trees/six_gene_tree/`: one 103-tip tree containing all seven genes,
  used to inspect gene-clade separation and flagged proteins. The directory name
  is historical and retained for stable links.
- `combined_trees/slrp_family_reference_tree/`: broader SLRP-family context,
  used for annotation diagnostics such as the chicken BGN/ASPN conflict.
- `per_gene_trees/`: large NCBI ortholog datasets and historical tree outputs.
  An audit found a header-parsing problem in the legacy filter; corrected
  per-locus candidates are ready, but corrected IQ-TREE reruns are pending.

Directories named `archive/` contain superseded or diagnostic runs, not the
tree to use in the thesis. For exact final filenames, use
`../overview/AUTHORITATIVE_RESULT_LOCATIONS.md`.

For current thesis figures, use `figures/compact_per_gene_trees.png` and
`figures/six_gene_combined_tree.png`. The large legacy trees are exploratory
and must not be used for final topology claims until their corrected reruns are
complete.
