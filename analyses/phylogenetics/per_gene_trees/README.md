# Large per-gene trees

Each gene directory uses the same explicit layout:

- `ncbi_download/`: raw NCBI ortholog download and extracted dataset.
- `candidates/`: raw and filtered candidate FASTAs plus reports. Files ending
  `_locus_filtered.faa` and `_locus_filter_report.tsv` are the corrected inputs.
- `alignments/`: MAFFT alignment.
- `trimmed_alignments/`: trimAl alignment used for IQ-TREE.
- `final_tree/`: supported final IQ-TREE output.
- `logs/`: retained run logs.

## Important audit status

The tree files currently under `final_tree/` and `itol/` are legacy exploratory
results. The original filter used the last bracketed NCBI header field as a
species value; in many headers that field was an isoform such as `X1`. This
could collapse unrelated loci globally. The filter now parses organism,
GeneID/locus, and isoform separately and supports `--one-per-locus`.

Corrected candidate counts are BGN 639, EPYC 632, FMOD 546, LUM 838, OGN 855,
and PRELP 848. To regenerate corrected large trees without overwriting legacy
provenance, run in WSL:

```bash
cd /path/to/Bachelor_Thesis
bash analyses/phylogenetics/per_gene_trees/scripts/run_corrected_large_per_gene_trees.sh
```

This is an overnight-scale job. Corrected outputs are written under each gene's
`corrected_final_tree/` directory. Until they exist and pass review, use the
curated compact trees for thesis conclusions.
