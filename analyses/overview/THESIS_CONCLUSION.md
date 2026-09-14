# Defensible overall thesis conclusion

Last audited: 2026-09-10

## Recommended central conclusion

Expression and literature identified seven growth-plate/cartilage-relevant
SLRPs whose curated vertebrate candidates generally conserve an LRR-rich SLRP
core, coding-exon organization, gene-associated sequence patterns, and
phylogenetic identity. Available SignalP calls support classical secretion for
85 of 86 tested current proteins; 17 new DCN/replacement-BGN calls remain
pending and must not be treated as either positive or negative. Together, the
evidence is consistent with long-term constraint on a secreted vertebrate ECM
toolkit, while documented exceptions reveal paralogy, incomplete gene models,
lineage-specific history, and annotation uncertainty.

The evolutionary-origin result is deliberately qualified. Published evidence
places an SLRP-like precursor early in chordates and major family expansion in
early vertebrates. The present amphioxus, lamprey, cartilaginous-fish, gar,
coelacanth, teleost, and tetrapod sampling brackets that history; it does not
treat any living species as an ancestor or reconstruct the ancestral sequence
of an individual gene. Gene-level assignments become most uncertain at the
amphioxus/lamprey boundary and in partial shark annotations, while a broadly
recognizable repertoire is supported across jawed vertebrates.

## Why this is more informative than “everything is conserved”

- Four historical BGN assignments were recognized as DCN or DCN-like records,
  demonstrating why close class-I paralogs require reciprocal, tree, and locus
  evidence rather than a single similarity hit.
- DCN, EPYC, FMOD, OGN, and PRELP are unrooted-monophyletic in the 103-tip tree.
- BGN is not one pure combined-tree split, but every BGN and DCN sequence has a
  closer same-gene than opposite-class-I neighbour. Partial shark BGN models
  remain appropriately tentative.
- LUM has a 15/16 main split; the lamprey record lies outside and remains a
  tentative deep-lineage assignment.
- OGN becomes coherent after excluding amphioxus as uncertain SLRP-like
  provenance. Spotted-gar OGN remains structurally plausible but SignalP
  negative.
- Coding-exon number, complete translation, and intron phase are stable, while
  gene span and UTR annotation vary more strongly.
- Expression patterns differ by dataset, species, age, and growth-plate zone;
  conservation of protein structure does not imply identical expression.
- Completed review of all 98 SynVoy rows yielded 62 accepted, 20 tentative, 12
  ambiguous, and four rejected locus decisions. Alternative retained loci
  corrected misleading top fragments, while unrelated amphioxus hits and the
  chicken BGN/ASPN conflict were rejected.
- The species summary showed a broad rise in median human-reference identity
  from lamprey/sharks toward mammals, but not a linear ladder: spotted gar
  exceeded zebrafish and coelacanth and retained five of seven accepted loci.
  This separates evolutionary divergence from gene-specific rate and annotation
  effects.
- All opossum synteny calls are excluded pending a matched FASTA/GFF rerun;
  independently supported opossum proteins remain usable.
- The deep-lineage species are evolutionary representatives, not ancestors;
  their value is in showing where modern paralog assignments become uncertain.

## Gene-level interpretation

| Gene | Thesis role | Integrated interpretation |
|---|---|---|
| BGN | principal class-I anchor | Highest local expression and strong direct growth-plate support; conserved domain/MSA/coding evidence, but partial shark models, chicken ASPN-like provenance, and non-monophyletic combined-tree topology require qualified wording. |
| DCN | principal class-I candidate/paralog control | Direct cartilage/ECM function, expression, phenotype, domain/MSA, 15/15 pure tree split, conserved seven-CDS-exon structure, and conserved reviewed loci; SignalP is still pending for the current panel. |
| FMOD | principal class-II anchor | Strong direct developmental growth-plate localization plus coherent protein/tree/coding evidence; duplicated zebrafish co-orthologs and the alternative spotted-gar locus are informative evolutionary findings. |
| PRELP | principal class-II candidate | Direct growth-plate expression and experimental bone-interface roles align with strong comparative conservation; zebrafish locus accepted, catshark tentative, amphioxus unresolved. |
| EPYC | principal class-III candidate | Direct growth-plate ECM localization and coherent comparative evidence; whale-shark/zebrafish are tentative and the amphioxus SynVoy call is rejected. |
| LUM | principal class-II candidate | Strong protein/coding evidence and accepted canonical shark/gar/coelacanth loci; lamprey placement and unresolved amphioxus require caution. |
| OGN | exploratory class-III candidate | Conserved vertebrate set but less consistent rat growth-plate expression and explicit deep-lineage/SignalP caveats. |
| OMD | context gene | Relevant to osteoblast/mineral/primary-spongiosa biology and conserved coding structure, but not run through the full protein/tree/synteny panel. |

## Claims to avoid

Do not state that all genes are equally conserved, that all analyses used the
same 15 species, or that conserved sequence proves conserved growth-plate
function in every species. Do not interpret pending SignalP as negative,
missing Bgee/SynVoy output as absence, MGI/HPO counts as effect sizes, motif
matches as binding sites, or an unrooted tree as evolutionary direction.

## Thesis-ready one-sentence version

“A multi-evidence comparison of seven growth-plate-associated SLRPs revealed a
strongly conserved vertebrate LRR/coding framework compatible with secreted ECM
function, together with gene-specific expression, paralog-resolution, and a
small set of informative deep-lineage and annotation uncertainties.”
