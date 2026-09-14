    # Gene evidence profiles

Last updated: 2026-09-03

> **Current panel:** DCN is a full comparative member and the decontaminated
> 13-protein BGN set has been approved and rebuilt. The historical filename is
> retained so existing links do not break.

The expanded comparative panel is **BGN, DCN, FMOD, PRELP, EPYC, LUM, and
OGN**. OMD is included below as a biologically relevant context gene, not as a
member of the full protein/tree/synteny panel. The profiles link literature-supported
function to the project's expression, protein, MSA, tree, synteny, and
gene-structure results. They do not infer that a known mammalian function is
experimentally demonstrated in every sampled species.

Use this file as the concise gene-by-gene bridge to the thesis results. For the
fuller biological explanation, experimental evidence, and explicit unknowns,
read [`GENE_FUNCTION_REFERENCE.md`](GENE_FUNCTION_REFERENCE.md) first or keep it
open alongside these profiles.

The current supplementary mature-protein motif analysis assigned all 103
seven-gene proteins most strongly to the model learned from their recorded gene.
Because the same small curated sets were used for model learning and evaluation,
this is supportive classification evidence rather than an independent accuracy
test. Curated MGI/HPO phenotypes are mentioned below as an independent
biological cross-check; their
counts reflect annotation depth and are not gene-effect sizes.

## BGN — class-I anchor

**Known function and relevance.** BGN is a secreted class-I SLRP involved in
collagen-rich skeletal ECM and signalling. In mouse genetics and biochemical
experiments, combined Bgn/Fmod deficiency caused low bone mass and increased
osteoclastogenesis; BGN and FMOD interacted with TNF-alpha/RANKL and attenuated
osteoclast formation ([Kram et al., 2017](https://doi.org/10.1038/s41598-017-12651-6)).

**Computational evidence.** BGN ranked first in the local P0 screen
(3012.762 TPM) and had strong mouse/rat growth-plate mean ranks (2.0/1.5).
All 13 canonical proteins had SLRP-like Pfam support; 11 have retained
SignalP-positive calls and the two replacement shark models are pending. Mean
forward identity was 80.69% and mean pairwise MSA identity was 72.14%. The
combined tree did not place all 13 in one pure split, but every BGN tip had a
closer BGN than DCN neighbour in the focused class-I diagnostic. Representative
coding structure was seven CDS exons in all five reference species.
MGI contained 14 skeletal/cartilage/joint terms for BGN, and the pinned HPO
release contained human disease-phenotype annotations for BGN.

The updated BGN SynVoy run passed all 14 genome-QC checks and retained 10 HIGH
and one MEDIUM annotation, with nine consistency flags versus 36 historically.
Best calls were HIGH in nine genomes. No post-filter candidate was retained in
amphioxus, catshark, chicken, elephant shark, or opossum. Detailed review
retained catshark `XP_038638277.1` tentatively using broad-block support,
rejected the chicken ASPN-like record, and left amphioxus, elephant shark, and
opossum ambiguous. The opossum result additionally requires rerun because its
FASTA/GFF sequence versions do not match.

**Curation and conclusion.** Chicken `XP_414298.2` was excluded from the protein
panel: despite a BGN locus label, the annotated product, reciprocal hit, and
focused tree supported ASPN. BGN remains the strongest biological and technical
anchor; the chicken record is an annotation/synteny conflict, not evidence
against vertebrate BGN conservation. Historical opossum, elephant-shark,
catshark, and whale-shark DCN/DCN-like records were also removed from BGN;
catshark `XP_038638277.1` and whale-shark `XP_048475905.1` are the tentative
replacement BGN models. The catshark locus has partial broad-block support, but
both partial proteins still require focused terminal/model review.

## DCN — class-I cartilage-matrix candidate and BGN paralog control

**Known function and relevance.** Decorin is a secreted class-I SLRP that binds
collagen-rich matrices and contributes to fibril organization. Developing-human
tissue localization showed DCN in type-I- and type-II-collagen-rich connective
tissues, including spatially restricted skeletal sites, and also showed that DCN
and BGN are not interchangeable expression markers
([Bianco et al., 1990](https://pubmed.ncbi.nlm.nih.gov/2212616/)). In decorin-null
mouse cartilage, the pericellular matrix had less aggrecan, mildly increased
collagen-II fibril diameter, reduced micromodulus, and impaired chondrocyte
calcium responses to mechanical/osmotic stimulation
([Chery et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7902451/)). These
experiments support a cartilage-matrix and mechanotransduction role, but they do
not establish that DCN alone controls longitudinal growth.

**Computational evidence.** DCN ranked fourth in the local P0 screen (581.246
TPM), with direct mouse/rat growth-plate mean ranks of 5.5/5.75 and detection in
all four summarized conditions in both species. All 15 curated proteins have
significant SLRP/LRR-family Pfam support. The MSA has 71.64% mean pairwise
identity, 85.9% mean occupancy, and 73 completely conserved ungapped columns.
The compact DCN tree completed with 1,000 ultrafast-bootstrap and 1,000 SH-aLRT
replicates. Across human, mouse, cow, chicken, and zebrafish, DCN retains seven
CDS blocks and all six homologous junctions retain intron phase. The four finite
human-target codon comparisons all have omega below one (median 0.109756),
compatible with purifying constraint. SignalP is not yet available for this
panel and must remain marked pending. ProtSpace assigns 15/15 DCN proteins to
the correct gene centroid, while opossum and catshark remain outlier priorities.
The combined compact tree recovers a pure 15/15 DCN split. The updated SynVoy
run passed all 14 target-genome QC checks and retained 10 HIGH and one MEDIUM
post-filter candidate. Amphioxus, spotted gar, and zebrafish had no retained
post-filter DCN candidate. Manual review nevertheless recovered canonical,
ordered DCN loci in coelacanth, spotted gar, and zebrafish; all remain tentative
until the pending SignalP layer is completed. Amphioxus is unresolved, and the
current opossum synteny result is excluded because of the input mismatch.

**Curation and conclusion.** Opossum `XP_001363160.3` is tentative because the
current protein is a 212-aa partial model. Catshark `XP_038636759.1` is a
tentative long isoform whose canonical SynVoy locus is confirmed but whose
isoform choice still needs direct sequence comparison. DCN is biologically justified in the expanded panel and is also
methodologically essential: distinguishing DCN from its class-I paralog BGN
exposed four contaminated or DCN-like historical BGN assignments. Current MGI
contains 24 DCN single-gene terms, including four skeletal/cartilage/joint
terms, and HPO contains 10 DCN phenotype terms; counts reflect database curation
depth rather than gene importance.

## FMOD — class-II growth-plate anchor

**Known function and relevance.** Developmental mouse data showed highest Fmod
expression in growing knee epiphyses and localized FMOD around late-hypertrophic
chondrocytes, the secondary ossification centre, and growth plate
([Säämänen et al., 2001](https://pmc.ncbi.nlm.nih.gov/articles/PMC1221771/)).
FMOD also contributes to collagen organization and skeletal remodelling together
with BGN.

**Computational evidence.** FMOD ranked third locally (605.864 TPM) and had
strong mouse/rat growth-plate mean ranks (2.0/1.5). All 14 nonredundant proteins
had SLRP-like domains and SignalP-positive N termini. Mean forward identity was
73.15%, mean MSA identity was 64.96%, and the 14/14 unrooted FMOD split had a
94.9/99 label. The representative coding structure was two CDS exons across all
five species.
MGI contributed four skeletal/cartilage/joint terms for FMOD.

**Curation and conclusion.** A historically duplicated lamprey sequence was
removed from FMOD and retained once under LUM. Zebrafish `fmodb` is the better
full-length human-FMOD match, whereas `fmoda` retains the ancestral FMOD--PRELP
locus; both are retained as co-ortholog evidence. Spotted-gar `fmoda` was
accepted at an alternative retained HIGH locus with eight shared human anchors.
Catshark remains tentative. The 684-aa elephant-shark FMOD-like record is a
likely compound/fused annotation: its N-terminal half is DCN/BGN-like and its
C-terminal half is FMOD-best, with a second LRR N-terminal cap and hydrophobic
segment near the join. Its locus is retained as tentative FMOD evidence, but
the protein is excluded from the confident sequence panel. The amphioxus
fragment was rejected as a BRD3-like locus. FMOD
is a strong class-II candidate with unusually direct growth-plate localization
evidence, but its deep-lineage models illustrate annotation and duplication
limits.

## PRELP — class-II growth-plate and bone-interface candidate

**Known function and relevance.** PRELP was characterized as a cartilage matrix
gene, and mouse studies conserved its LRR, cysteine, and glycosylation-related
features ([Grover et al., 1996](https://doi.org/10.1006/geno.1996.0605);
[Grover et al., 2002](https://doi.org/10.1016/S0945-053X(01)00165-2)). Its
N-terminal region can inhibit osteoclastogenesis
([Rucci et al., 2009](https://doi.org/10.1083/jcb.200906014)), while cell work
also links PRELP to osteoblast differentiation.

**Computational evidence.** PRELP ranked sixth locally (189.840 TPM), but direct
mouse/rat growth-plate mean ranks were strong (2.0/3.0). All 16 proteins had
SLRP-like domain and SignalP support. Mean forward identity was 75.36%, mean MSA
identity was 66.20%, and all 16 formed a pure unrooted split (96/100 label).
Two CDS exons were conserved across the five reference species.
MGI contained PRELP phenotype annotations, but none fell into the broad
skeletal/cartilage/joint category in the downloaded release; this does not
override the direct expression and experimental evidence.

**Curation and conclusion.** Zebrafish PRELP was accepted at a locus adjacent
to `fmoda`; catshark remains tentative because its annotation is sparse, and a
gene-specific amphioxus locus remains unresolved. The current opossum synteny
call is excluded pending a matched-input rerun. PRELP shows why direct
anatomical support can outweigh moderate abundance in the broader juvenile-
lineage sample.

## EPYC — class-III cartilage-specific candidate

**Known function and relevance.** EPYC expression accompanies cartilage
development, and protein was localized throughout embryonic growth-plate ECM
around resting, proliferating, and hypertrophic chondrocytes
([Johnson et al., 1999](https://pubmed.ncbi.nlm.nih.gov/10633869/)). This is the
panel's clearest direct anatomical literature rationale for a class-III member.

**Computational evidence.** EPYC ranked eighth locally (99.059 TPM) and had
intermediate mouse/rat growth-plate ranks (5.0/4.75). All 15 proteins had
SLRP-like domains and SignalP-positive N termini. Mean forward identity was
72.05%, mean MSA identity was 62.17%, and the 15/15 unrooted split had a
92.8/82 label. Six CDS exons were conserved in the five-species structure
screen.
MGI independently recorded short-femur and osteoarthritis phenotypes for EPYC.

**Curation and conclusion.** Whale-shark `XP_048463897.1` was retained as
tentative. Manual inspection found a similar early N-terminal beginning and a
strongly conserved middle-to-C-terminal region, but numerous earlier mismatches,
reduced query coverage, an alternative reciprocal hit, and an EPYC-to-OGN
ProtSpace tendency prevent an unqualified call. EPYC remains a justified main
candidate with this sequence-level caveat recorded explicitly. Locus review
accepted coelacanth, elephant-shark, and spotted-gar EPYC, retained zebrafish
tentatively, and rejected the amphioxus SynVoy interval because it overlaps an
unrelated RUN-domain protein-like gene.

## LUM — strong class-II candidate with completed automated synteny

**Known function and relevance.** LUM is associated with skeletal matrix
development and cartilage collagen deposition/fibrillogenesis
([Raouf et al., 2002](https://doi.org/10.1016/S0945-053X(02)00027-6);
[Kafienah et al., 2008](https://doi.org/10.1016/j.matbio.2008.04.002)). These
functions fit a secreted LRR-rich ECM protein.

**Computational evidence.** LUM ranked fifth locally (272.457 TPM), with mouse
and rat growth-plate mean ranks of 4.5. All 16 proteins had SLRP-like domains and
SignalP-positive N termini. Mean forward identity was 68.00% and mean MSA
identity was 60.21%. Coding structure was especially stable: two CDS exons and
339--345-aa estimated proteins across all five reference species.
MGI was dominated by corneal and skin/connective-tissue phenotypes for LUM and
included one tendon-related skeletal term.

**Curation and conclusion.** The main unrooted split contained 15/16 LUM tips;
lamprey `XP_075930353.1` lay among FMOD-like class-II sequences. It is retained
once as a tentative LUM-like deep-lineage assignment. The completed 14-target
SynVoy run passed genome QC for all targets and recovered a best HIGH call in
10 genomes and a best MEDIUM call in four. Manual review accepted the canonical
catshark, coelacanth, and spotted-gar loci using ordered flanking blocks plus
protein/tree support even though the run summary selected different fragments.
The amphioxus interval was intergenic and remains ambiguous. LUM is therefore a
strong main vertebrate candidate, not a claim of 14 proven one-to-one orthologs.

## OGN — exploratory vertebrate class-III candidate

**Known function and relevance.** BMP1/Tolloid processing enhances OGN's
regulation of type-I collagen fibrillogenesis
([Ge et al., 2004](https://doi.org/10.1074/jbc.M406630200)). Vertebrate analyses
show conservation alongside teleost duplication and subfunctionalization
([Costa et al., 2018](https://doi.org/10.1186/s12862-018-1310-2)), making OGN
useful for studying both conserved ECM architecture and lineage-specific history.

**Computational evidence.** OGN ranked second locally (657.412 TPM), but its
direct mouse/rat growth-plate mean ranks were weaker (7.5/7.0), and rat lacked
detection in all replicates. After manual curation, the confident vertebrate set
contains 14 proteins. Mean forward identity is 69.82%, mean MSA identity is
58.04%, and all 14 form a pure unrooted OGN split with a 97.7/100 label.
Thirteen of 14 have SLRP-like Pfam support and 13/14 are SignalP-positive. Six
CDS exons are conserved across the five structure-screen species.
MGI supported broader corneal, skin/collagen, cardiovascular, and metabolic
phenotypes, but no term fell into the broad skeletal category in this release.

**Manual decisions.** Opossum `XP_056660002.1` and zebrafish
`NP_001013588.1` were retained based on intact termini, conserved cysteines, and
continuous LRR-rich alignment. Spotted-gar `XP_006631165.2` was retained as
tentative because the aligned protein is structurally plausible but lacks a
SignalP-positive N terminus. Amphioxus `XP_066265713.1` was excluded from the
confident OGN set because its divergent termini, large sequence-specific
insertions, alternative reciprocal hit, tree position, and synteny did not
support gene-specific orthology; it remains only as uncertain SLRP-like
provenance. The sequence is from *Branchiostoma lanceolatum*, whereas the
SynVoy target was *B. floridae*, so it cannot be mapped onto the reviewed
SynVoy assembly. The opossum sequence decision remains valid at protein level,
but its current synteny score is excluded until the mismatched opossum target is
rerun.

**Conclusion.** OGN is the least direct growth-plate candidate but is worth
retaining as an explicitly exploratory vertebrate member. The exclusions and
tentative calls are an informative comparative result, not a failed analysis.

## OMD — skeletal context gene, not full-panel member

**Known function and relevance.** OMD/osteoadherin is secreted by osteoblasts,
binds hydroxyapatite, supports integrin-mediated osteoblast attachment, and was
localized to fetal growth-plate primary spongiosa
([Sommarin et al., 1998](https://pmc.ncbi.nlm.nih.gov/articles/PMC2132750/)).
OMD can also enhance BMP2/SMAD osteogenic signalling, with terminal LRRs
contributing to BMP2 interaction
([Lin et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7862363/)).

**Computational evidence.** OMD ranked ninth locally (96.317 TPM), had a mouse
growth-plate mean rank of 7.75, and had no all-replicate rat detection in this
dataset. Its representative coding structure was stable across human, mouse,
cow, chicken, and zebrafish: two CDS exons, estimated proteins of 402--424 aa,
and no flagged structural outlier.

**Scope and conclusion.** OMD did not undergo the canonical protein/domain,
SignalP, MSA, phylogeny, ProtSpace, or SynVoy workflow. It should therefore be
used as relevant cartilage--bone-interface context, not presented as an eighth
co-equal comparative candidate unless the full pipeline is deliberately added.

## Final selection argument

The original six genes were selected by integrating juvenile-lineage expression, direct
zone-resolved growth-plate evidence, known skeletal/cartilage biology,
representation across SLRP classes, and technical orthology/conservation quality.
The audit adds DCN because it meets the same biological and computational logic
and provides the critical class-I paralog control for BGN. The expanded
interpretation supports BGN, DCN, FMOD, PRELP, EPYC, and LUM as principal genes,
with OGN as a transparent exploratory member. RNA abundance alone is neither a
sufficient candidate-selection rule nor evidence of orthology.
