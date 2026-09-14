# Canonical compact SLRP protein dataset

Last rebuilt and verified: 2026-09-14

## Scope

The current canonical comparative panel contains **103 nonredundant proteins**
from seven genes:

| Gene | Proteins |
|---|---:|
| BGN | 13 |
| DCN | 15 |
| EPYC | 15 |
| FMOD | 14 |
| LUM | 16 |
| OGN | 14 |
| PRELP | 16 |

OMD is retained as an expression, literature, gene-structure, and coding-
constraint context gene. It is not part of this canonical protein/tree panel.

The authoritative concatenated FASTA and manifest are:

- `candidates/slrp_candidates_canonical.faa`
- `candidates/canonical_candidate_manifest.tsv`

## Per-gene source FASTAs

- BGN: `candidates/BGN_domain_input_canonical_decontaminated.faa`
- DCN: `candidates/DCN_domain_input.faa`
- EPYC: `candidates/EPYC_domain_input.faa`
- FMOD: `candidates/FMOD_domain_input_lamprey_deduplicated.faa`
- LUM: `candidates/LUM_domain_input.faa`
- OGN: `candidates/OGN_domain_input_dog_corrected.faa`
- PRELP: `candidates/PRELP_domain_input.faa`

## Applied curation decisions

- Chicken `XP_414298.2` is excluded from canonical BGN because its product and
  reciprocal/tree evidence are ASPN-like despite the BGN locus label.
- Historical opossum `XP_001363160.2`, elephant-shark `NP_001279248.1`,
  catshark `XP_038636760.1`, and whale-shark `XP_048463588.1` must not be
  interpreted as BGN; they are DCN records or DCN-like models.
- Catshark BGN `XP_038638277.1` and whale-shark BGN `XP_048475905.1` are
  retained as tentative partial models and require focused manual review.
- Dog OGN `XP_038383338.1` is analyzed after removing an erroneous 64-aa
  N-terminal extension; the original sequence remains provenance.
- Spotted-gar OGN `XP_006631165.2` is retained as tentative and is the only
  known SignalP-negative canonical protein.
- Amphioxus `XP_066265713.1` is excluded from the confident OGN set and retained
  only as uncertain SLRP-like provenance.
- Lamprey `XP_075930353.1` is retained once under LUM, not FMOD, and remains a
  tentative deep-lineage LUM-like assignment.
- Whale-shark EPYC `XP_048463897.1` is retained as tentative following manual
  alignment review.

## Verified current outputs

- Pfam/HMMER: 693 significant hits and 103 per-protein summaries.
- MAFFT/conservation: seven alignments, 103 QC rows, consensus tables,
  identity matrices, conservation profiles, and seven sequence logos.
- Integrated protein QC: 84 Supported, 13 Watch, and 6 Needs inspection; no
  protein remains in the SignalP-rerun category.
- SignalP: retained calls are available for all 103 canonical proteins (100
  positive and three negative). The former 17-protein upload is preserved at
  `signal_peptides/input/pending_17_canonical_for_signalp.faa`, and its completed
  raw job bundle is archived under `signal_peptides/raw_outputs/`.
- ProtSpace: 103 canonical proteins plus four controls were embedded; canonical
  QC contains 85 Supported, 1 Watch, and 17 Needs-inspection calls.
- Protein motifs: 116 MEME motifs, 35 STREME motifs, and 103/103 proteins with
  their assigned gene model ranked best by MAST. The 2026-09-14 rerun removed
  the signal peptide from all 100 final SignalP-positive proteins and left only
  the three SignalP-negative proteins untrimmed. This is supportive
  classification evidence, not validation of biochemical binding sites.
- Compact combined tree: 103 tips, 172 trimmed columns, LG+I+G4, 1,000
  SH-aLRT and 1,000 ultrafast-bootstrap replicates. DCN, EPYC, FMOD, OGN, and
  PRELP are unrooted-monophyletic. LUM is 15/16 because lamprey lies outside
  the main split. BGN is not a single pure split after DCN is included, but a
  direct BGN/DCN nearest-neighbour diagnostic assigns every BGN and DCN tip to
  a closer same-gene tip.

## Authoritative result locations

- Protein-level integrated table:
  `tables/protein_conservation_domain_msa_signalp_summary.tsv`
- Gene-level protein summary:
  `tables/protein_conservation_domain_msa_signalp_gene_summary.tsv`
- Final candidate decisions:
  `../overview/tables/final_candidate_level_confidence.tsv`
- Final gene-level evidence:
  `../overview/tables/final_gene_level_evidence.tsv`

Legacy folders and filenames containing `six_gene` or `90` are retained for
provenance or stable links. Their names do not define the current dataset.

## Rebuild policy

`scripts/run_canonical_refresh.sh` rebuilds the manifest, Pfam parsing/scan,
MAFFT alignments, conservation outputs, SignalP merge, and integrated tables
from the source FASTAs above. SignalP is an external step and missing calls must
not be invented or inferred. Do not hand-edit generated summary
tables; change the curated source/decision record and rebuild the affected
downstream products.
