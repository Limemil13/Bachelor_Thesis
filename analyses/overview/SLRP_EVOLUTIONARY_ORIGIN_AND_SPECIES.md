# SLRP evolutionary origin, species roles, and outgroups

## Thesis-safe conclusion

No living species in the project is the species "from which" SLRPs originated.
Living tunicates, amphioxus, lamprey, sharks, and bony vertebrates are sampled
representatives whose genomes bracket stages of chordate and vertebrate
evolution. The defensible synthesis is:

> An ancestral SLRP-like protein was probably present early in chordate
> evolution. The family then expanded during early vertebrate evolution, and a
> broadly recognizable modern SLRP repertoire was established before the
> divergence of cartilaginous and bony vertebrates.

This wording does not claim that sharks are ancestral to other living
vertebrates and does not claim that the project reconstructed the ancestral
sequence of each gene.

## What the origin literature sampled

| Evolutionary position | Representative species | Evidence and interpretation | Role in this thesis |
|---|---|---|---|
| Tunicates (urochordates) | *Ciona intestinalis*, *Ciona savignyi* | Park et al. reported one LRRCE-containing SLRP-like gene in each *Ciona* genome and proposed it as a possible representative of the precursor of canonical vertebrate SLRPs. | Literature-only family-origin context. It cannot be assigned confidently to BGN, DCN, FMOD, PRELP, LUM, EPYC, or OGN and therefore should not be inserted into every gene-specific analysis. |
| Cephalochordates | *Branchiostoma floridae* | The 2008 genome-era search did not recover a clear canonical LRRCE-containing SLRP. The authors discussed either emergence after the cephalochordate split or secondary loss; annotation incompleteness is an additional modern caveat. | Provides an uncertainty boundary, not evidence that amphioxus is the ancestor or that all amphioxus SLRP-like proteins are orthologs of modern genes. |
| Jawless vertebrates | sea lamprey, *Petromyzon marinus* | Several distinguishable SLRP-like genes were reported, including biglycan-, decorin-, keratocan-, and epiphycan-like forms, indicating substantial early vertebrate diversification. | Deep-lineage context. The project retains lamprey BGN, PRELP/keratocan-like, and LUM-like assignments with gene-specific caveats. |
| Cartilaginous fishes | elephant shark/chimaera, *Callorhinchus milii*; whale shark, *Rhincodon typus* | Costa et al. recovered most human SLRP family members and inferred conservation of the family at least since the Chondrichthyes–Osteichthyes split. Failure to find BGN, OPTC, and ECMX in the 2018 shark annotations was a search/annotation result, not proof of absolute biological absence. | Earliest main-panel jawed-vertebrate branch. Partial or conflicting shark BGN models must remain tentative, which is consistent with the older literature's uncertainty. |
| Non-teleost bony fishes | spotted gar, *Lepisosteus oculatus*; coelacanth, *Latimeria chalumnae* | These lineages help separate ancient vertebrate conservation from teleost-specific duplication. | Core evolutionary bridge in the project panel. |
| Teleosts | zebrafish, *Danio rerio*, plus other teleosts in Costa et al. | Costa et al. supported teleost-specific duplication of OGN, FMOD, and ECM2; `ogn1` is syntenic with human OGN, whereas `ogn2` is fish-specific. | Explains why teleost paralogy must be checked rather than treating the top protein hit automatically as the one-to-one ortholog. |
| Tetrapods | frog, anole, chicken, mouse, human and other mammals | Provide the best-annotated reference framework for gene-specific orthology and conserved structure/function. | Main reference and ingroup comparisons. Human is the query/locus anchor, not an evolutionary ancestor of the other species. |

Primary sources: [Park et al. 2008](https://doi.org/10.1186/1471-2164-9-599)
and [Costa et al. 2018](https://doi.org/10.1186/s12862-018-1310-2).

## Quantitative findings from this project

The reproducible species summary is
`tables/species_lineage_evidence_summary.tsv`; its figure is
`figures/species_lineage_evidence_summary.png`.

- Median identity to the human query generally increased from the deepest
  vertebrate comparisons toward mammals: 48.5% across the three retained
  lamprey proteins, 51.2% in elephant shark, 53.2% in catshark, and
  85.0--93.5% in non-human mammals.
- The pattern was not a simple evolutionary ladder. Spotted gar had a median
  of 68.1%, compared with 61.2% in zebrafish and 59.8% in coelacanth. These
  values also reflect gene-specific rates, different retained gene subsets,
  and annotation/model quality.
- Spotted gar and coelacanth each had five accepted synteny loci out of seven,
  making them especially useful bridges between the deepest comparisons and
  tetrapods. The three shark species each had three accepted loci, while their
  remaining calls exposed partial or compound gene models and annotation
  uncertainty.
- Zebrafish retained all seven canonical genes, but `fmoda` and `fmodb` provide
  complementary sequence and locus evidence. This is a direct project example
  of why teleost duplication cannot be reduced to one top BLAST hit.
- Chicken retained six canonical proteins. The rejected BGN-locus product is
  ASPN-like; this is an annotation/paralog finding and has nothing to do with a
  supposed absence of bone in birds.
- Opossum demonstrates a technical rather than evolutionary failure: its six
  canonical proteins have a median identity of 85.0%, but all seven synteny
  rows are unusable because the FASTA and GFF sequence versions are
  incompatible.
- Amphioxus produced no confident gene-specific protein in the canonical panel
  and no accepted synteny row. Because the protein and synteny work used
  different *Branchiostoma* species and deep annotations are incomplete, this
  defines an uncertainty boundary rather than proving absence of these genes.

The thesis-level finding is therefore not only “deep species are more
different.” The comparison separates genuine lineage-specific events, such as
teleost duplication, from paralog confusion, compound gene models, incomplete
annotation, and incompatible genome inputs.

## Mapping to the current project species

The clean protein manifest contains 103 proteins from seven genes and up to 16
species labels overall. Not every gene has a defensible sequence in every
species. Direct overlap with the Costa sampling includes human, mouse, chicken,
anole, frog, zebrafish, spotted gar, coelacanth, elephant shark, and whale
shark. Opossum supplies marsupial coverage analogous to the Tasmanian devil in
Costa et al.; catshark supplies an additional cartilaginous-fish comparison.

The project's amphioxus proteome is *Branchiostoma lanceolatum*, whereas Park et
al. examined *B. floridae*. The excluded putative OGN sequence
`XP_066265713.1` is annotated as chondroadherin-like, is divergent in the MSA,
and lacks convincing OGN tree/synteny support. It is therefore an informative
uncertain SLRP-like/other-SLRP case, not a confident OGN ortholog.

## Outgroup versus outlier

- An **outgroup** is deliberately chosen to root a tree or polarize character
  change. It should be homologous but outside the focal ingroup.
- An **outlier** is an unusually divergent observation. A thesis dataset does
  not need to contain one, and an outlier must never be manufactured merely to
  make the results look more interesting.
- The compact gene-specific and combined trees are valid as unrooted evidence
  of clustering. Their left-to-right display does not show ancestry.
- If a rooted OGN tree is required, one or more related class-III paralogs such
  as EPYC/OPTC can serve as the gene-level outgroup. Costa et al. used other
  SLRPs as the OGN outgroup.
- If a broad family-origin tree is required, a separately curated external LRR
  homolog such as LINGO3 can root it; Costa et al. rooted their broad family
  analysis using whale-shark LINGO3. This should be a separate context figure,
  not mixed into the per-gene conservation statistics.

## Real uncertainty examples already available

The project already contains genuine edge cases, so no artificial outlier is
needed:

- amphioxus `XP_066265713.1`: excluded from confident OGN; retain only as
  uncertain SLRP-like/chondroadherin-like provenance;
- spotted-gar OGN `XP_006631165.2`: structurally plausible but tentative because
  SignalP is negative;
- lamprey LUM-like `XP_075930353.1`: outside the main LUM split and therefore
  retained as a tentative deep-lineage assignment;
- catshark BGN `XP_038638277.1` and whale-shark BGN `XP_048475905.1`: reciprocal
  BGN support but partial/reduced-coverage gene models requiring cautious
  interpretation.

These cases demonstrate the limit of sequence-only orthology assignment and
justify combining domains, MSA, phylogeny, synteny, gene structure, and manual
review.
