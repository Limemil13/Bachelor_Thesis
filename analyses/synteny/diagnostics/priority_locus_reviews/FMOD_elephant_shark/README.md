# Elephant-shark FMOD-like model diagnostic

`XP_007897806.2` was reconstructed directly from the nine CDS records in the
exact elephant-shark FASTA/GFF3 pair used by SynVoy. The reconstruction is
684 aa, begins with methionine, ends at a terminal stop, and contains no
internal stop codons.

This is not a routine long FMOD ortholog. Local BLASTP splits the model into
two SLRP-like units: residues 1--330 align best to human DCN (47.0% identity),
whereas residues 376--684 align best to human FMOD (57.0% identity and 82.2%
human-reference coverage). Pfam detects LRRs across both units and significant
LRR N-terminal caps at residues 54--76 and 383--412. The protein also has a
hydrophobic N terminus and a second hydrophobic segment at residues 330--347.

The conservative interpretation is a likely compound/fused gene annotation.
The syntenic `LOC103182549` locus remains tentative FMOD evidence, but this
protein is not part of the confident FMOD protein/MSA/tree panel. Resolution
would require transcript evidence or a revised gene model.

## Reproducible outputs

- `XP_007897806.2.faa`: CDS-derived protein.
- `XP_007897806.2_qc.tsv`: translation and phase QC.
- `XP_007897806.2_reference_blastp.tsv`: comparisons with eight human SLRPs.
- `XP_007897806.2_pfam.domtblout`: raw Pfam domain table.
- `XP_007897806.2_pfam.log`: complete HMMER log.
- `XP_007897806.2_diagnostic_summary.tsv`: one-row interpretation summary.
- `pairwise_comparison.tsv`: human--elephant-shark neighbourhood comparison.
- `neighborhoods.tsv`: oriented exact-GFF gene windows.

The scripts are `../../../scripts/extract_gff_protein.py` and
`../../../scripts/compare_protein_to_reference_queries.py`. The Pfam scan used
HMMER 3.4 and `data/domains/Pfam-A.hmm`.
