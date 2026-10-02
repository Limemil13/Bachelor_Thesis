# Functional reference for the SLRP genes

Last updated: 2026-09-03

## Purpose and reading order

This is the detailed biological reference for understanding the genes in the
thesis. It is broader than `SIX_GENE_PROFILES.md`: this file explains what the
proteins are thought to do, which experiments support those functions, and what
remains unknown. `SIX_GENE_PROFILES.md` is the shorter integration layer that
connects this biology to the project's expression, protein, MSA, tree, synteny,
and gene-structure results.

The expanded comparative panel is **BGN, DCN, FMOD, PRELP, EPYC, LUM, and
OGN**. **OMD** is included as a relevant cartilage--bone-interface context gene
but has not undergone the full protein/tree/synteny workflow. DCN was promoted
on 2026-09-02 after its expression and established cartilage biology were
integrated with a completed protein, tree, and gene-structure package.

This is a thesis-focused synthesis, not a claim to cover every publication in
every organ system. Evidence is labelled by experiment type because a binding
assay, cultured-cell experiment, expression pattern, knockout phenotype, and
cross-species computational result answer different questions.

The "local P0" shorthand below means the project's stored Salmon
quantifications of three WT P0 mouse samples from public NCBI GEO GSE305415
(BioProject PRJNA1305489), not an original local sequencing study. The
independent direct growth-plate dataset is GSE114919.

## Shared molecular framework

Small leucine-rich proteoglycans (SLRPs) are secreted extracellular-matrix
(ECM) proteins with a central series of leucine-rich repeats (LRRs), usually
flanked by cysteine-rich terminal regions. The LRR domain creates a curved
protein-interaction surface. Different family members carry different
glycosaminoglycan (GAG) chains or N-linked oligosaccharides, and some occur as a
mixture of glycoprotein and proteoglycan forms. The modern five-class system is
based on sequence homology, domain organization, genomic organization, and
chromosomal relationships rather than on GAG chemistry alone
([Iozzo and Schaefer, 2015](https://doi.org/10.1016/j.matbio.2015.02.003)).

Three shared ideas are important for this thesis:

1. **Secretion is necessary but not sufficient evidence of ECM function.** A
   conserved N-terminal signal peptide supports entry into the secretory
   pathway. SignalP does not prove the final tissue localization, binding
   partner, or physiological effect of the mature protein.
2. **LRR conservation supports an interaction scaffold.** SLRPs can organize
   collagen fibrils and can modulate growth-factor or receptor pathways, but a
   conserved LRR domain does not prove that all orthologs bind the same partners
   with the same affinity
   ([Schaefer and Iozzo, 2008](https://doi.org/10.1074/jbc.R800020200)).
3. **Family members can compensate for one another.** Mild single-knockout
   phenotypes and stronger double-knockout phenotypes occur in this family.
   Therefore, a mild knockout does not necessarily mean that a protein is
   unimportant; it may reflect redundancy, developmental compensation, or an
   untested physiological challenge. Conversely, a double knockout does not
   identify the separate contribution of each gene.

The project tests whether these molecular systems are evolutionarily conserved.
It does **not** experimentally test collagen binding, GAG modification, receptor
signalling, or growth-plate function in each species.

## BGN — biglycan

### What the molecule is

BGN is a class-I SLRP. Its mature extracellular core protein commonly carries
chondroitin-sulfate or dermatan-sulfate chains near the N-terminus. It is
abundant in skeletal and other collagen-rich connective tissues. BGN is closely
related to decorin (DCN), but the two proteins are not interchangeable in every
tissue or pathway.

### Functions supported by experiments

**Collagen matrix and tissue mechanics.** BGN contributes to the organization
of collagen fibrils. Bgn-deficient mice show abnormal fibril morphology in bone
and other connective tissues, reduced skeletal growth and bone mass, while
combined Bgn/Dcn loss produces a much stronger bone and connective-tissue
phenotype. This supports both a direct matrix role and partial compensation
among class-I SLRPs
([Corsi et al., 2002](https://pubmed.ncbi.nlm.nih.gov/12102052/)).

**Bone formation and remodelling.** BGN deficiency delays new bone formation
after marrow ablation, consistent with impaired osteogenic capacity
([Chen et al., 2003](https://doi.org/10.1007/s00223-002-1101-y)). BGN can also
influence BMP/TGF-beta-associated osteoblast behaviour and the balance between
osteoblast and osteoclast activity, but the direction and importance of these
effects depend on age, sex, tissue, and experimental model. Bgn/Fmod
double-deficient mice have low bone mass and increased osteoclastogenesis; BGN
and FMOD can bind TNF-alpha and RANKL and suppress osteoclast formation in the
experimental system used by Kram and colleagues
([Kram et al., 2017](https://doi.org/10.1038/s41598-017-12651-6)). OPG-Fc only
partially rescued the mature double-knockout skeleton and was harmful in young
mice, showing that the phenotype cannot be reduced to a single RANKL pathway
([Kram et al., 2020](https://doi.org/10.1016/j.jsb.2020.107627)).

**Developmental distribution.** Developmental expression studies place Bgn in
multiple mesenchymal and skeletal compartments, including sites of cartilage
and bone formation, but also at non-skeletal sites. This argues against treating
BGN as a uniquely growth-plate-specific marker
([Wilda et al., 2000](https://doi.org/10.1359/jbmr.2000.15.11.2187)).

### What is not fully known

BGN has both structural ECM effects and signalling-associated effects. Their
relative importance in a particular growth-plate zone is not fully resolved.
Many strong skeletal observations come from whole-body single or compound
knockouts, so primary growth-plate effects can be difficult to separate from
changes in bone cells, tendons, joint mechanics, inflammation, or compensation
by other SLRPs. The literature does not demonstrate that BGN performs every
reported mammalian function in all vertebrate orthologs.

### Relation to this thesis

BGN is the strongest expression and comparative anchor in the project. It
ranked first in the local P0 expression screen. Its decontaminated 13-protein
set has SLRP-domain support throughout, 12 SignalP-positive calls and one
negative fragmentary whale-shark model, a continuous LRR-rich alignment core, and
stable seven-CDS-exon organization. The seven-gene tree does not recover BGN as
one pure split, but a focused BGN/DCN nearest-neighbour analysis assigns every
BGN and DCN tip to a closer same-gene tip. This convergence supports a preserved
BGN-like class-I matrix-interaction scaffold while retaining caution for partial
deep-lineage models. It does not by itself prove conserved
osteoclast or BMP signalling. The chicken `XP_414298.2` record was excluded
because protein and tree evidence supported ASPN rather than BGN, illustrating
why locus labels alone are insufficient.

An accession-level audit removed historical opossum, elephant-shark, catshark,
and whale-shark DCN/DCN-like records from BGN. Catshark `XP_038638277.1` and
whale-shark `XP_048475905.1` are tentative BGN replacements. Catshark BGN has a
strong signal peptide but an incomplete C terminus; whale-shark BGN is a
discontinuous fragment without a detected signal peptide. This correction is an
informative result: protein similarity within SLRP class I is insufficient
without reciprocal, tree, and locus evidence.

## DCN — decorin

### What the molecule is

DCN is a class-I SLRP and one of BGN's closest paralogs. It is synthesized with
an N-terminal signal peptide, secreted, and commonly carries one dermatan- or
chondroitin-sulfate chain near the N terminus. Its curved LRR-rich core binds
components of collagenous ECM. Similarity to BGN makes DCN a biologically
important comparison and an orthology challenge: sequence similarity alone can
confuse the two genes in incompletely annotated species.

### Functions supported by experiments

**Collagen organization.** Dcn-null mice have abnormal collagen fibril shape,
diameter distributions, and tissue mechanics, establishing an in-vivo role in
fibril organization rather than only an in-vitro binding interaction
([Danielson et al., 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC2134287/)).
The strength of the phenotype varies by tissue, and other SLRPs can occupy some
of the same collagen-associated sites.

**Cartilage pericellular matrix and mechanotransduction.** In decorin-null mouse
articular cartilage, the chondrocyte pericellular matrix contained less
aggrecan, showed a mild increase in collagen-II fibril diameter, and had a lower
micromodulus. Chondrocyte calcium responses were also impaired. These results
link DCN to the local mechanical/charge environment around chondrocytes and to
mechanotransduction
([Chery et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7902451/)).

**Developmental skeletal distribution.** Human developmental localization
showed DCN in collagen-rich connective tissues and in spatially restricted
skeletal compartments, including epiphyseal/articular regions. DCN and BGN
patterns were sometimes divergent or mutually exclusive, supporting related
but non-identical functions
([Bianco et al., 1990](https://pubmed.ncbi.nlm.nih.gov/2212616/)).

### What is not fully known

DCN participates in many ECM and signalling contexts outside the growth plate,
so it is not a growth-plate-specific marker. Articular-cartilage mechanics,
developing epiphyseal localization, and growth-plate RNA expression support
relevance, but they do not demonstrate that DCN is a unique regulator of
longitudinal growth. The relative contributions of the LRR core, its GAG chain,
and compensation by BGN or other SLRPs depend on tissue and developmental stage.

### Relation to this thesis

DCN is reproducibly expressed in the local P0 Prx1-lineage samples and in the
independent mouse/rat growth-plate dataset. Its 15-species curated set has
SLRP/LRR Pfam support throughout and a globally conserved MSA. A compact
gene-specific tree is complete. The five-species structure comparison retains
seven CDS blocks and conserved intron phase at all six homologous junctions;
finite human-target omega estimates are below one. Together these results
support a conserved secreted ECM scaffold and stable ortholog models, while not
proving identical mechanotransduction effects across vertebrates. SignalP is
positive for 14 of 15 DCN proteins; the negative opossum record is a conserved
C-terminal fragment that lacks the signal-peptide-containing N terminus. The updated SynVoy run is
complete, with ten HIGH and one MEDIUM retained candidates; amphioxus,
coelacanth, spotted gar, and zebrafish received detailed locus review. Canonical
ordered loci support coelacanth, spotted gar, and zebrafish, while amphioxus
remains unresolved. The current opossum synteny score is excluded because of
the FASTA/GFF mismatch. Focused review retained the partial opossum protein and
accepted the full-span catshark isoform, with its N-terminal extension and later
SignalP cleavage recorded as a model caveat.

## FMOD — fibromodulin

### What the molecule is

FMOD is a class-II SLRP and collagen-binding extracellular proteoglycan. Its
core can carry keratan-sulfate or related N-linked oligosaccharides, and its
N-terminal region contains additional modifications that may affect matrix
interactions. It is prominent in tendon, cartilage, bone, and other fibrillar
collagen matrices.

### Functions supported by experiments

**Collagen fibrillogenesis.** Fmod-null mice have disorganized tendon collagen,
irregular fibril outlines, a shift toward thin fibrils, and altered LUM
deposition. This is direct in-vivo evidence that FMOD participates in collagen
fibril assembly and that LUM can respond to FMOD loss
([Svensson et al., 1999](https://doi.org/10.1074/jbc.274.14.9636)).

**Growth plate and developing joint.** In the developing mouse knee, Fmod
expression was highest in growing epiphyses, and FMOD protein localized around
late-hypertrophic chondrocytes, the growth plate, and the secondary ossification
centre. Expression and distribution changed with age
([Säämänen et al., 2001](https://pmc.ncbi.nlm.nih.gov/articles/PMC1221771/)).
This is more direct growth-plate evidence than is available for most genes in
the panel.

**Cartilage remodelling and disease.** FMOD is cleaved during cartilage
development and pathology. MMP-13-cleaved FMOD occurs near hypertrophic
chondrocytes in fetal growth plates; in biochemical experiments FMOD was
processed by MMP-13 and ADAMTS-4. These findings connect FMOD to the remodelling
stage of endochondral ossification, but cleavage is not automatically evidence
that the fragments have a specific signalling function
([Shu et al., 2019](https://doi.org/10.3390/ijms20030579)). Fmod-null and
Bgn/Fmod compound-null models also develop joint, tendon, and bone phenotypes,
including increased osteoarthritis susceptibility and altered osteoclast
activity.

### What is not fully known

FMOD and LUM have overlapping collagen-regulatory roles, so the isolated effect
of FMOD depends on tissue and developmental stage. Some articular-cartilage
degeneration in Fmod-null mice may be secondary to altered ligaments or joint
mechanics rather than a primary defect in chondrocytes. The precise function of
each FMOD proteolytic fragment in the growth plate is not established.

### Relation to this thesis

FMOD combines high juvenile/growth-plate expression with unusually direct
developmental localization. All 14 curated proteins retain SLRP-family and
SignalP support, form a strongly supported gene-specific tree split, and share
a representative two-CDS-exon organization. The results support a conserved
secreted collagen-associated FMOD scaffold. They cannot establish identical
proteolysis or compensatory interaction with LUM in each species.

## PRELP — proline/arginine-rich end leucine-rich repeat protein

### What the molecule is

PRELP is a class-II extracellular LRR protein distinguished by a strongly basic,
proline/arginine-rich N-terminal region. It is found in cartilage and several
other connective tissues. Depending on tissue and terminology, it is described
as a glycoprotein or proteoglycan-family member rather than as a uniformly
GAG-substituted molecule.

### Functions supported by experiments

**Matrix-to-cell-surface linkage.** The PRELP N-terminus binds heparin and
heparan sulfate with nanomolar affinity. Fibroblast binding to PRELP is inhibited
by heparin, supporting a model in which PRELP can connect extracellular matrix
to cell-surface heparan-sulfate proteoglycans
([Bengtsson et al., 2000](https://doi.org/10.1074/jbc.M007917200)). This is a
biochemical and cell-binding model; its importance in a living growth plate has
not been directly quantified.

**Cartilage distribution and conservation.** Human PRELP was characterized as
an articular-cartilage matrix gene, and mouse PRELP retains the LRR, cysteine,
and glycosylation-related features while being expressed in fetal and postnatal
cartilage
([Grover et al., 1996](https://doi.org/10.1006/geno.1996.0605);
[Grover and Roughley, 2002](https://doi.org/10.1016/S0945-053X(01)00165-2)).
The timing of expression differs between reports and species, so the human age
pattern should not be generalized to all vertebrates.

**Bone-cell regulation.** An N-terminal PRELP-derived peptide can enter
osteoclast-lineage cells, reduce NF-kappaB activity, inhibit osteoclastogenesis,
and counteract bone loss in several mouse models
([Rucci et al., 2009](https://doi.org/10.1083/jcb.200906014);
[Rucci et al., 2013](https://doi.org/10.1002/jbmr.1951)). PRELP knockdown in an
osteoblast-lineage cell system also reduced osteogenic differentiation and
mineralization
([Li et al., 2016](https://doi.org/10.1016/j.bbrc.2016.01.106)).

### What is not fully known

The pharmacological activity of an isolated N-terminal peptide does not prove
that endogenous full-length PRELP is released, internalized, or active at the
same concentration in normal growth plates. The proposed heparan-sulfate linker
function, osteoclast inhibition, and osteoblast effects have not yet been
assembled into one definitive in-vivo mechanism. PRELP therefore has strong
biochemical plausibility and skeletal relevance, but its exact physiological
role in longitudinal bone growth remains incompletely resolved.

### Relation to this thesis

PRELP had moderate abundance in the broad local P0 sample but strong direct
mouse/rat growth-plate ranks. All 16 proteins show SLRP-family, SignalP,
alignment, tree, and two-CDS-exon support. Conservation of the mature LRR core
supports a preserved interaction scaffold; conservation of a secretion signal
is compatible with ECM localization. The present analysis did not separately
measure conservation of every basic N-terminal residue or demonstrate
heparan-sulfate/osteoclast activity across species.

## EPYC — epiphycan (PG-Lb)

### What the molecule is

EPYC is a class-III SLRP and a chondroitin/dermatan-sulfate proteoglycan enriched
in cartilage. Its historical name PG-Lb reflects its original biochemical
description. Compared with BGN, FMOD, and LUM, fewer molecular binding partners
and signalling mechanisms have been established for EPYC.

### Functions supported by experiments

**Developmental cartilage localization.** During mouse embryonic limb
development, Epyc transcription is associated mainly with resting and
proliferating chondrocytes, whereas EPYC protein is distributed throughout the
growth-plate ECM, including around hypertrophic chondrocytes. This difference is
biologically reasonable because a secreted matrix protein can persist and
spread beyond the cells currently producing its RNA
([Johnson et al., 1999](https://pubmed.ncbi.nlm.nih.gov/10633869/)). The study
directly supports growth-plate relevance but proposed, rather than demonstrated,
a role in regulating collagen or other matrix components.

**Joint maintenance and genetic interaction.** Epyc-null mice are viable and
grossly normal at birth. With age, male single-knockout mice showed shorter
femurs and many deficient animals developed osteoarthritis. The phenotype was
earlier or stronger in Epyc/Bgn double-deficient mice, with compensatory changes
in other SLRP transcripts. These data support a role in joint integrity and an
in-vivo interaction with BGN, but also show that loss of EPYC alone is not a
severe early skeletal-development defect
([Nuka et al., 2010](https://doi.org/10.1016/j.joca.2009.11.006)).

### What is not fully known

EPYC's precise molecular mechanism in the growth plate remains one of the least
resolved in the panel. Localization and knockout data support cartilage-matrix
relevance, but they do not identify a single indispensable receptor or binding
pathway. The mild early phenotype can reflect redundancy, a maintenance role
that emerges with age, or effects that appear only under mechanical or disease
challenge.

### Relation to this thesis

EPYC is justified by direct anatomical evidence rather than by being among the
highest RNA hits. Its 15 proteins share SignalP and SLRP-family support and form
a gene-specific tree split; six CDS exons are conserved in the five-species
structure screen. These results strengthen ortholog identity and support a
conserved secreted cartilage-proteoglycan architecture. Whale-shark
`XP_048463897.1` is retained as tentative: its N-terminal beginning and
middle-to-C-terminal core are compatible with EPYC, but numerous earlier
mismatches, reduced coverage, and an alternative reciprocal hit prevent an
unqualified assignment.

## LUM — lumican

### What the molecule is

LUM is a class-II collagen-binding SLRP. It can occur with keratan-sulfate or
other N-linked carbohydrate substitutions, depending on tissue. It is widely
distributed in collagen-rich tissues including cornea, skin, tendon, cartilage,
and bone, which means it is biologically relevant but not growth-plate-specific.

### Functions supported by experiments

**Collagen fibril regulation.** Lum-null mice and Lum/Fmod compound-null models
show abnormal collagen fibril diameter or shape in several tissues. FMOD and LUM
have overlapping but non-identical functions: FMOD is thought to act strongly
at earlier fibril-assembly stages, while LUM limits later lateral fibril growth
in some tissues. Compensation is tissue dependent, so this model should not be
treated as a universal sequence of events.

**Cartilage and skeletal development.** LUM changes during osteoblast
differentiation and mineralization and is present in developing cartilaginous
and bone matrices
([Raouf et al., 2002](https://doi.org/10.1016/S0945-053X(02)00027-6)). In
tissue-engineered cartilage, LUM knockdown increased type-II collagen
accumulation and fibril diameter, supporting a context-dependent restraint on
collagen deposition
([Kafienah et al., 2008](https://doi.org/10.1016/j.matbio.2008.04.002)). LUM is
also present in fetal cartilage rudiments and growth plates and is less readily
cleaved than FMOD in the protease systems examined by Shu and colleagues
([Shu et al., 2019](https://doi.org/10.3390/ijms20030579)).

**Bone resorption.** Recombinant LUM inhibited osteoclast differentiation,
migration/fusion, and resorption in vitro by suppressing Akt activity
([Lee et al., 2021](https://doi.org/10.3390/ijms22094717)). This supports a
possible osteoprotective role, but it is not yet equivalent to demonstrating a
normal growth-plate function in vivo.

### What is not fully known

LUM has broad and context-dependent effects: the same collagen-limiting activity
can be beneficial for controlled fibril assembly yet restrict collagen
accumulation in engineered cartilage. Redundancy with FMOD complicates
single-gene interpretation. Its direct contribution to juvenile growth-plate
elongation, as distinct from cartilage, tendon, bone, or osteoclast biology, is
not fully established.

### Relation to this thesis

LUM has consistent expression support and 16 curated proteins with conserved
SignalP and SLRP-family architecture. Fifteen of 16 lie in the main LUM tree
split; lamprey `XP_075930353.1` remains a tentative deep-lineage LUM-like
assignment because it groups near FMOD-like sequences. The conserved two-CDS-
exon organization is strong. The completed 14-target LUM SynVoy analysis passed
QC for every target and yielded best HIGH calls in 10 genomes and MEDIUM calls
in four. Manual review accepted the canonical catshark, coelacanth, and
spotted-gar loci using ordered flanking blocks together with protein and tree
support, even though the automated summary selected other fragments. The
amphioxus interval was intergenic and remained ambiguous. LUM should therefore
be described as a strong vertebrate multi-evidence candidate whose exceptions
also show why automated candidate rank is not equivalent to orthology.

## OGN — osteoglycin (mimecan)

### What the molecule is

OGN is a class-III secreted SLRP, historically also called mimecan. It is found
in several connective tissues and can be proteolytically processed. It is not a
growth-plate-specific protein, and the name “osteoglycin” should not be taken as
proof of a unique or exclusive bone function.

### Functions supported by experiments

**Collagen fibrillogenesis.** Ogn-null mice are viable and fertile and lack an
obvious gross developmental phenotype under standard conditions, but their skin
and corneal collagen fibrils are thicker and more variable and skin tensile
strength is moderately reduced. This is direct evidence that OGN regulates
collagen fibrillogenesis in vivo while also showing that compensation limits the
gross phenotype
([Tasheva et al., 2002](https://pubmed.ncbi.nlm.nih.gov/12432342/)).

**Proteolytic activation.** BMP1/Tolloid-family metalloproteinases process
pro-OGN. In biochemical collagen assays, processed OGN more strongly regulated
type-I collagen fibrillogenesis than the precursor
([Ge et al., 2004](https://doi.org/10.1074/jbc.M406630200)). This provides a
specific reason to take terminal sequence integrity seriously, but it does not
demonstrate growth-plate-specific activity or establish that processing is
identical in every species.

**Evolutionary history.** OGN-family evolution includes teleost gene
duplication and subfunctionalization. Consequently, a divergent fish sequence
can be a real OGN-family ortholog or co-ortholog, while a similar sequence from
a deep lineage is not automatically a gene-specific OGN ortholog
([Costa et al., 2018](https://doi.org/10.1186/s12862-018-1310-2)).

### What is not fully known

OGN has convincing general collagen-matrix functions, but the direct evidence
for a specific role in the juvenile growth plate is weaker than for FMOD or
EPYC. Other studies connect OGN to bone formation, metabolism, cardiovascular
remodelling, and inflammatory contexts, but these broader effects do not yet
form one universally accepted mechanism and are not required to justify the
growth-plate thesis. It remains uncertain which reported non-matrix effects are
primary, tissue-specific, or secondary to ECM changes.

### Relation to this thesis

OGN ranked highly in the broad local P0 screen but had weaker direct rodent
growth-plate ranks. The manually curated confident set contains 14 vertebrate
proteins and forms a supported OGN tree split. Opossum and zebrafish were
retained; spotted gar is tentative because SignalP did not support its
N-terminus; amphioxus was removed from the confident OGN set because alignment,
tree, reciprocal-hit, and synteny evidence did not jointly support OGN
orthology. OGN is therefore best presented as an exploratory vertebrate
candidate: its conservation analysis is informative, but its specific
growth-plate role should not be overstated.

## OMD — osteomodulin/osteoadherin (context gene)

### What the molecule is and what it does

OMD is a class-II SLRP also called osteoadherin. It is an acidic, secreted,
mineral-associated proteoglycan produced by osteoblast-lineage cells. Purified
osteoadherin binds hydroxyapatite and supports osteoblast attachment through an
alpha-v-beta-3-integrin-associated mechanism. In fetal bovine growth plates it
was localized to the primary spongiosa rather than to all cartilage zones
([Wendel et al., 1998](https://doi.org/10.1083/jcb.141.3.839)). This places OMD
at the mineralizing cartilage--bone interface.

In osteogenic cell and mandibular-defect models, OMD enhanced BMP2/SMAD
signalling and bone formation, and terminal LRR regions contributed to BMP2
interaction
([Lin et al., 2021](https://doi.org/10.1038/s41419-021-03404-5)). These models
support an osteogenic function but do not prove a direct role in longitudinal
growth-plate cartilage.

### What is not fully known and relation to this thesis

OMD's mineral binding, cell adhesion, and BMP2-associated effects are
biologically coherent, but their relative contributions in normal development
are not fully separated. In this project OMD has expression and representative
gene-structure data only; it did not undergo the canonical protein/domain,
SignalP, MSA, phylogeny, ProtSpace, or SynVoy workflow. It should remain a
context gene unless the complete comparative pipeline is deliberately added.

## The defensible biological conclusion

The genes do not all have the same level of functional knowledge:

- **BGN and FMOD** have the strongest combined skeletal genetics, collagen-
  matrix, and remodelling evidence.
- **FMOD and EPYC** have the clearest direct growth-plate localization evidence,
  although EPYC's molecular mechanism is less well defined.
- **PRELP** has a distinctive and experimentally supported N-terminal
  heparan-sulfate/osteoclast biology, but its normal in-vivo growth-plate
  mechanism remains incomplete.
- **LUM** has strong collagen-fibril and skeletal-matrix evidence, with important
  context dependence and overlap with FMOD.
- **OGN** has convincing collagen-fibrillogenesis biology and an informative
  vertebrate evolutionary history, but the weakest direct growth-plate-specific
  functional evidence among the main six.
- **OMD** is a strong cartilage--bone-interface context gene, not a completed
  member of the comparative panel.

The thesis therefore should not conclude merely that “everything is conserved.”
Its stronger conclusion is that an expression-guided set of SLRPs retains the
molecular features expected for secreted ECM proteins across vertebrates, while
integrated sequence, tree, gene-structure, and synteny evidence exposes specific
annotation, paralogy, and deep-lineage exceptions. Literature makes these genes
biologically plausible candidates; the computational work establishes ortholog
plausibility and conservation, not conserved function by itself.
