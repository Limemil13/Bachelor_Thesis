# Manual review of SynVoy loci

Last updated: 2026-09-21

**Current status:** all 98 SynVoy gene-by-species rows are reviewed: 26 detailed
P1 loci, 69 exact-locus P2 confirmations and four P3 controls. All live pending
queues contain zero rows. Use
`diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv` for the
complete curated decisions and `tables/synvoy_gene_species_evidence.tsv` for
the merged evidence view. The checklist below is retained for reproducibility
and for future rerun exceptions.

This completed review was not 98 separate visual GDV confirmations. All 98 rows
received exact-assembly GFF/accession checks; the 26 P1 exceptions received
detailed feature-overlap and oriented-neighbour comparisons, while the 68 P2
and four P3 rows received reproducible accession-to-locus confirmation.
Your GDV screenshots and notes informed selected difficult cases; the exact-GFF
checks are the reproducible basis of the 98-row audit. GDV remains the
recommended interface for a reader's own visual spot-checks. A row
without an archived GDV screenshot should not be described as visually
confirmed in GDV. New SynVoy reports or corrected opossum inputs reopen only
the rows whose underlying evidence changes.

**Opossum warning:** the current opossum SynVoy results are invalid for synteny
interpretation because the GFF chromosome versions do not match the paired
FASTA sequence IDs. Opossum proteins can still be used in protein/MSA/tree
analyses, but no current opossum SynVoy score should be cited until that target
is rerun with a synchronized input pair.

## Quick start: what you do for every species

Use this section while you are reviewing. The later sections explain difficult
cases in more detail.

### The nine-step checklist

| Step | Action | What to write down |
|---:|---|---|
| 1 | Choose the next P1 row from the priority table below. | Gene, species, SynVoy run, candidate count. |
| 2 | Open NCBI Genome Data Viewer (GDV) and select the same species and assembly used by SynVoy. | Assembly and chromosome/scaffold accession. |
| 3 | Search the exact SynVoy coordinate. If there is no coordinate, search the protein accession and then the gene symbol. | What gene/model overlaps the location: expected symbol, another gene, `LOC...`, fragment, or no model. |
| 4 | Search the canonical protein accession separately. | Whether the protein maps to the same locus, a different nearby locus, or nowhere. |
| 5 | Record the gene model. | Confirmed symbol, strand, transcript/protein accession, and whether the model is complete, partial, or fragmented. |
| 6 | Zoom out to roughly 100 kb on each side. | Two or three informative non-SLRP neighbours on the left and right, if the assembly permits it; also note their order/orientation. |
| 7 | Open the SynVoy SVG/HTML for the same species and gene. | Whether SynVoy and GDV show the same locus and flanking block. |
| 8 | Compare the protein/MSA/tree evidence already present in the project. | `phylogeny_consistent`: `yes`, `no`, or `not available`. No new tree run is required for this check. |
| 9 | Make one final call using all evidence together. | `accepted`, `tentative`, `rejected`, or `ambiguous`, followed by one or two evidence-based sentences. |

Do not call a gene absent merely because its name is missing from GDV or
SynVoy. Missing annotation, a fragmented assembly, or a scaffold boundary must
first be ruled out.

### What to record for each row

| Field | Allowed entry | Meaning |
|---|---|---|
| `confirmed_gene_symbol` | Symbol, `LOC...`, `uncharacterized`, or `not found` | The annotation actually shown at the reviewed locus. |
| `confirmed_protein_accession` | RefSeq accession or `not available` | The protein linked to that locus, not merely the expected protein. |
| `coordinate_to_accession_status` | `same locus`, `different locus`, `not mapped`, or `not assessable` | Whether the SynVoy coordinate and canonical protein resolve to the same place. |
| `neighbor_order_consistent` | `yes`, `partial`, `no`, or `not assessable` | Whether informative flanking anchors support the corresponding locus. A whole reversed block can still be consistent. |
| `phylogeny_consistent` | `yes`, `no`, or `not available` | Whether the candidate falls in the expected gene-specific clade. |
| `review_status` | `accepted`, `tentative`, `rejected`, or `ambiguous` | Overall decision using locus, protein, and tree evidence. |
| `manual_notes` | One or two sentences | State the positive evidence and the unresolved conflict. |
| `reviewer` / `review_date` | Your name or initials / `YYYY-MM-DD` | Audit trail for the thesis. |

Enter completed decisions in
`analyses/synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`.
The master evidence and P1/P2 queue files are generated views and must not be
edited as separate sources of truth.
The filled template and supporting screenshots can be used to update the
master row reproducibly.

### Copy-and-fill template

```text
Gene/species:
Assembly:
SynVoy coordinate:
SynVoy candidate count:
GDV symbol/model:
Protein accession:
Strand/model completeness:
Left flanking anchors:
Right flanking anchors:
Coordinate and accession: same locus / different locus / not mapped / not assessable
Neighbour order: yes / partial / no / not assessable
Tree: yes / no / not available
Decision: accepted / tentative / rejected / ambiguous
Manual note:
Reviewer/date:
```

### Exact assemblies used by SynVoy

These accessions come from the `#!genome-build-accession` lines in the exact
GFF files paired with the SynVoy FASTA inputs. In NCBI GDV, search/select the
assembly accession rather than relying only on a similar-looking build name.

| SynVoy label | Species | NCBI assembly accession | Build name | Full NCBI GDV URL | Role |
|---|---|---|---|---|---|
| human | *Homo sapiens* | `GCF_000001405.40` | GRCh38.p14 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000001405.40` | Home/reference locus |
| mouse | *Mus musculus* | `GCF_000001635.27` | GRCm39 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000001635.27` | Target |
| dog | *Canis lupus familiaris* | `GCF_011100685.1` | UU_Cfam_GSD_1.0 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_011100685.1` | Target; alternate assembly |
| cow | *Bos taurus* | `GCF_002263795.3` | ARS-UCD2.0 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_002263795.3` | Target |
| opossum | *Monodelphis domestica* | `GCF_027887165.2` | mMonDom1.pri | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_027887165.2` | Target; rerun required because the local FASTA/GFF sequence IDs mismatch |
| chicken | *Gallus gallus* | `GCF_016699485.2` | bGalGal1.mat.broiler.GRCg7b | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_016699485.2` | Target |
| anole | *Anolis carolinensis* | `GCF_035594765.1` | rAnoCar3.1.pri | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_035594765.1` | Target |
| frog | *Xenopus tropicalis* | `GCF_000004195.4` | UCB_Xtro_10.0 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000004195.4` | Target |
| zebrafish | *Danio rerio* | `GCF_049306965.1` | GRCz12tu | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_049306965.1` | Target |
| spotted gar | *Lepisosteus oculatus* | `GCF_040954835.1` | fLepOcu1.hap2 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_040954835.1` | Target |
| coelacanth | *Latimeria chalumnae* | `GCF_037176945.1` | fLatCha1.pri | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_037176945.1` | Target |
| elephant shark | *Callorhinchus milii* | `GCF_018977255.1` | IMCB_Cmil_1.0 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_018977255.1` | Target |
| catshark | *Scyliorhinus canicula* | `GCF_902713615.2` | sScyCan1.2 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_902713615.2` | Target |
| whale shark | *Rhincodon typus* | `GCF_021869965.1` | sRhiTyp1.1 | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_021869965.1` | Target |
| amphioxus | *Branchiostoma floridae* | `GCF_000003815.2` | Bfl_VNyyK | `https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000003815.2` | Target/deep-lineage comparison |

### Completed P1 review outcome

| Gene / species | Status | Canonical locus or conclusion |
|---|---|---|
| BGN / amphioxus | ambiguous | No gene-specific model; absence cannot be inferred. |
| BGN / catshark | tentative | `LOC119955784`, `NC_052166.1:4337032-4377582`, `XP_038638277.1`; partial broad-block support. |
| BGN / chicken | rejected | `XP_414298.2` is ASPN-like by product, focused tree and locus neighbours. |
| BGN / elephant shark | ambiguous | The only SynVoy hit was DCN; no annotated BGN model was resolved. |
| BGN / opossum | ambiguous | Conserved surrounding block but no BGN model; input mismatch requires rerun. |
| DCN / amphioxus | ambiguous | No gene-specific DCN locus resolved. |
| DCN / coelacanth | accepted | `NC_088145.1:114048728-114126931`, `XP_006004204.1`; eight ordered anchors and completed positive SignalP call. |
| DCN / spotted gar | accepted | `NC_090702.1:16396890-16423580`, `XP_015208366.1`; eight ordered anchors and completed positive SignalP call. |
| DCN / zebrafish | accepted | `NC_133179.1:16955805-16993841`, `NP_571772.1`; conserved EPYC-KERA-LUM-DCN block and completed positive SignalP call. |
| EPYC / amphioxus | rejected | SynVoy interval overlaps a RUN-domain protein-like gene, not EPYC. |
| EPYC / coelacanth | accepted | `NC_088145.1:113783413-113830371`, `XP_006009490.1`; ten ordered anchors. |
| EPYC / elephant shark | accepted | `NW_024704746.1:33285071-33298586`, `XP_007893502.1`; canonical rescue locus and seven ordered anchors. |
| EPYC / spotted gar | accepted | `NC_090702.1:16496431-16523089`, `XP_006633678.1`; 11 ordered anchors. |
| EPYC / zebrafish | tentative | `NC_133179.1:16906814-16920982`, `XP_021330689.1`; conserved block but reduced protein query coverage. |
| FMOD / amphioxus | rejected | SynVoy fragment lies in a BRD3-like gene. |
| FMOD / catshark | tentative | `NC_052160.1:142715179-142734639`, `XP_038675394.1`; adjacent PRELP but sparse annotation. |
| FMOD / zebrafish | tentative | `fmodb` is the fuller sequence match; `fmoda` retains the ancestral FMOD-PRELP block. |
| LUM / amphioxus | ambiguous | SynVoy fragment is intergenic; no gene-specific LUM model resolved. |
| LUM / catshark | accepted | `NC_052165.1:37258201-37272456`, `XP_038636708.1`; conserved KERA-LUM-DCN block. |
| LUM / coelacanth | accepted | `NC_088145.1:113890958-113908026`, `XP_006009488.1`; nine ordered anchors. |
| LUM / spotted gar | accepted | `NC_090702.1:16450380-16457870`, `XP_015208177.1`; nine ordered anchors. |
| OGN / amphioxus | rejected | All raw *B. floridae* calls map to other genes; reviewed sequence is from *B. lanceolatum*. |
| PRELP / amphioxus | ambiguous | No gene-specific PRELP locus resolved. |
| PRELP / catshark | tentative | `NC_052160.1:142404483-142432733`, `XP_038675391.1`; FMOD adjacency, sparse annotation. |
| PRELP / zebrafish | accepted | `NC_133186.1:23181387-23192171`, `XP_001923590.2`; adjacent `fmoda` and ordered anchors. |

Counts: 7 accepted, 8 tentative, 6 ambiguous and 4 rejected. See
`diagnostics/priority_locus_reviews/README.md` for the method, main findings and
file map.

<details>
<summary>Archived original P1 work order and questions</summary>

### Original P1 review order

This table is generated from the current 25-row P1 queue. A coordinate shown in
the third column is the **best candidate retained after SynVoy ownership and
paralog filtering**, unless the row explicitly says otherwise. `None retained`
means that the final report has no accepted HIGH/MEDIUM GOI coordinate; a GOI
label may still be visible in the provisional HTML. In that case, use the
listed NCBI protein or annotated locus rather than treating the HTML label as a
confirmed ortholog.

| Order | Gene / species | Coordinate or first search | Candidates | Main question |
|---:|---|---|---:|---|
| 1 | LUM / coelacanth | `NC_088145.1:113634160-113634471`; `XP_006009488.1` | 2 | Do the best fragment and canonical LUM protein map to the same locus? Inspect both candidate loci. |
| 2 | LUM / spotted gar | `NC_090702.1:17048413-17048754`; `XP_015208177.1` | 2 | Resolve the two fragment-level candidates and compare with the annotated LUM locus. |
| 3 | LUM / catshark | `NC_052156.1:105386109-105386660`; `XP_038636708.1` | 2 | Check both candidates, the partial model, and the many paralog/self-consistency flags. |
| 4 | LUM / amphioxus | `NC_049983.1:12602582-12603646`; then `LUM` | 1 | Is there a defensible LUM-like model or only a deep-lineage fragment? |
| 5 | BGN / chicken | Search `BGN`, `ASPN`, then the candidate locus | 0 | Distinguish a true BGN locus from the previously observed ASPN-like product conflict. |
| 6 | OGN / amphioxus | Search `OGN`, `osteoglycin`, and `mimecan` | 0 | Document whether the locus is only SLRP-like; do not force it into vertebrate OGN. |
| 7 | EPYC / spotted gar | **None retained.** NCBI EPYC: `NC_090702.1:16496431-16523089` (`XP_006633678.1`). Rejected provisional HTML hits: `NC_090721.1:6915943-6916575` (reassigned to KERA) and `NC_090702.1:17107093-17107848` (reassigned to LUM). | 0 | **Completed—accepted.** Eleven exact human flanking genes are shared and all 11 preserve oriented order; the canonical protein also has reciprocal, tree, domain, and SignalP support. |
| 8 | FMOD / zebrafish | `XP_073767056.1`; then `FMOD` | 0 | Check for an annotated/duplicated FMOD locus and conserved neighbours. |
| 9 | PRELP / zebrafish | `XP_001923590.2`; then `PRELP` or `prolargin` | 0 | Does the supported protein map to a defensible PRELP locus? |
| 10 | PRELP / amphioxus | Search `PRELP`, `prolargin`, then SLRP-like models | 0 | Is there any locus-level support, or is the case annotation-limited? |
| 11 | DCN / coelacanth | `NC_088145.1:114049722-114102359`; `XP_006004204.1` | 1 | Completed: accepted after ordered-flank, protein/tree, and positive SignalP review. |
| 12 | DCN / spotted gar | `XP_015208366.1`; then `DCN` | 0 | Locate the protein-supported gene and determine why SynVoy retained no locus. |
| 13 | DCN / zebrafish | `NP_571772.1`; then `DCN` | 0 | Locate the established protein and assess teleost locus rearrangement/annotation. |
| 14 | DCN / amphioxus | Search `DCN`, `decorin`, then SLRP-like models | 0 | Determine whether evidence is absent or simply not assessable in this deep lineage. |
| 15 | BGN / catshark | `XP_038638277.1`; then `BGN` | 0 | Completed: tentative partial BGN; signal peptide positive, C terminus incomplete, and locus support partial. |
| 16 | BGN / elephant shark | Search `BGN`, `biglycan`, then expected neighbours | 0 | Check whether an unlabelled/fragmented BGN locus exists. |
| 17 | BGN / opossum | Search `BGN`, `biglycan`, then expected neighbours | 0 | Check the apparent missing SynVoy call in a well-annotated mammal. |
| 18 | BGN / amphioxus | Search `BGN`, `biglycan`, then SLRP-like models | 0 | Treat absence cautiously; assess whether the locus is informative at all. |
| 19 | EPYC / elephant shark | `NW_024704760.1:2095489-2096193`; `XP_007893502.1` | 5 | Inspect every available candidate locus; do not decide from the best fragment alone. |
| 20 | EPYC / coelacanth | `NC_088145.1:113783923-113792675`; `XP_006009490.1` | 1 | Confirm the partial model and resolve the ProtSpace concern against locus/tree support. |
| 21 | EPYC / zebrafish | `NC_133200.1:20614002-20614649`; `XP_021330689.1` | 1 | Determine whether the fragment and protein-QC warning still support EPYC. |
| 22 | EPYC / amphioxus | `NC_049983.1:12339566-12340792`; then `EPYC` | 1 | Assess whether the partial deep-lineage locus is ortholog-specific or only SLRP-like. |
| 23 | FMOD / catshark | `NC_052156.1:108595627-108595680`; `XP_038675394.1` | 1 | Reconcile the tiny fragment with the full protein and ProtSpace concern. |
| 24 | FMOD / amphioxus | `NC_049986.1:10317928-10317972`; then `FMOD` | 1 | Decide whether the tiny fragment supplies meaningful locus support. |
| 25 | PRELP / catshark | `NC_052160.1:136025837-136038669`; `XP_038675391.1` | 1 | Confirm the partial locus and resolve the ProtSpace concern. |

Four current P1 rows contain more than one retained candidate: LUM in
coelacanth, spotted gar, and catshark, and EPYC in elephant shark. For these
rows, the single `best_*` coordinate is not enough. Search the canonical protein
and gene name as well and keep the call tentative until the alternative loci
have been compared.

</details>

#### Why the spotted-gar EPYC coordinates disagree

The supplied spotted-gar GFF contains the annotated `epyc` gene and transcript
at `NC_090702.1:16496431-16523089`; its CDS encodes `XP_006633678.1`. SynVoy's
HTML instead displays two provisional MEDIUM hits named `GOI_EPYC`. They are
sequence-and-neighbourhood candidates produced before the final family-ownership
check, not two annotated EPYC genes. The final report reassigns the hit on
`NC_090721.1:6915943-6916575` to `gene-KERA` and the hit on
`NC_090702.1:17107093-17107848` to `gene-LUM`. Consequently the authoritative
post-filter result is `best_confidence = NONE`, with no retained EPYC
coordinate.

This is a SynVoy detection miss rather than an assembly mismatch: both the NCBI
page and the input GFF use assembly `GCF_040954835.1`. The annotated EPYC locus
was present in the input annotation but was not recovered among the selected
EPYC candidate regions. The provisional HTML is therefore useful for inspecting
why a hit was considered, but `synvoy_report.json` and
`tables/synvoy_gene_species_evidence.tsv` are the sources of truth for the final
call.

The completed manual comparison supports the canonical locus strongly. Human
and spotted-gar EPYC share 11 exact flanking-gene symbols, all 11 remain in
consistent order after orienting the blocks by EPYC strand, and the local
`EPYC-KERA-LUM-DCN` SLRP cluster is preserved. `ATP2B1/atp2b1a` and
`KITLG/kitlga` provide two additional positionally corresponding anchors whose
symbols differ only because of teleost-style suffixes. With the existing
reciprocal EPYC hit, EPYC-clade placement, LRR-domain architecture, and strong
SignalP result, `XP_006633678.1` is **accepted as spotted-gar EPYC** despite the
automated SynVoy miss. Reproducible tables and the exact command are in
`diagnostics/spotted_gar_epyc_neighborhood/`.

The pasted chat classified BGN zebrafish as a multi-candidate P1 case. That is
out of date: the current master table classifies BGN zebrafish as a P3 positive
control with one complete HIGH call (`NP_001001825.2`).

### Decision rule in one table

| Decision | Use when |
|---|---|
| `accepted` | Locus/neighbours, protein identity, and tree all support the same ortholog. |
| `tentative` | The main evidence agrees, but one important weakness remains, such as a partial model, multiple loci, weak neighbours, or missing SignalP. |
| `rejected` | Another paralog/product is better supported, or locus, protein, and tree jointly contradict the assignment. |
| `ambiguous` | The assembly or annotation does not provide enough evidence to choose. This is not proof of gene absence. |

### Files and SynVoy result folders

| Gene | SynVoy result folder |
|---|---|
| BGN | `results/bgn_human_15species_dev_20260905` |
| DCN | `results/dcn_human_15species_dev_20260903` |
| EPYC | `results/epyc_human_15species` |
| FMOD | `results/fmod_human_15species` until the stopped refresh is resumed and completed later |
| LUM | `results/lum_human_15species` |
| OGN | `results/ogn_human_15species_standard_1408` |
| PRELP | `results/prelp_human_15species_fresh` |

From WSL, open a result directory in Windows Explorer with, for example:

```bash
cd ~/projects/bachelor_thesis/synvoy_tool/SynVoy
explorer.exe "$(wslpath -w results/lum_human_15species)"
```

For each reviewed case, keep or send three views: (1) the exact SynVoy
coordinate in GDV, (2) the wider approximately +/-100 kb GDV neighbourhood,
and (3) the matching SynVoy species plot. If the protein accession maps
elsewhere, also capture that locus.

## When to review

SynVoy reports are parsed for BGN, DCN, FMOD, PRELP, EPYC, LUM, and OGN. The
All 98 current rows are complete. New review is needed only if a target is
rerun, if a candidate/protein set changes, or if a new coordinate/protein
conflict appears. The opossum target specifically requires a clean rerun with
matching FASTA and GFF sequence IDs before it can be reviewed.

The authoritative curated-decision file is:

`analyses/synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`

The 98-row master table and the P1/P2 files are generated views. Rebuild them
from the curated-decision file rather than entering decisions only in a queue.

## What synteny can establish

SynVoy supplies locus-level support: it asks whether a candidate occurs in a
genomic neighbourhood compatible with the human reference locus. It does not
prove gene function and it does not independently distinguish every close SLRP
paralog. A final orthology call should combine:

1. target-locus coordinates and annotation;
2. order and orientation of flanking genes;
3. protein accession, length, domains, and SignalP;
4. reciprocal protein hit;
5. placement in the gene-specific tree.

## Optional visual review after the completed tables

The P2 and P3 coordinate/accession checks are complete. If time permits, open a
small representative subset in GDV for visual figure-quality confirmation:
one clean accepted locus, one tentative fragmented/duplicate locus, and one
rejected paralog. Escalate any new symbol/protein disagreement to the full P1
procedure.

## Detailed NCBI procedure for a P1 row

### 1. Start from the exact assembly and locus

Open [NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/gdv/) and select the
same species and assembly version used by SynVoy. Search for the chromosome or
scaffold accession and the `best_start`--`best_end` range from the master table.
If SynVoy has no retained coordinate, search for the canonical protein
accession, gene symbol, and known neighbouring genes instead.

Do not compare coordinates between different assembly versions. A coordinate is
meaningful only for the assembly from which it was generated.

### 2. Inspect annotation and genomic context

Turn on the RefSeq/NCBI gene, transcript, and protein tracks. Zoom out until the
nearest informative genes on both sides are visible. Record:

- the annotated symbol or `LOC...`/uncharacterized model;
- gene strand and approximate locus extent;
- transcript and protein accession, if present;
- at least two conserved flanking anchors when the assembly permits it;
- whether order and orientation agree exactly, agree after a small local
  rearrangement, or do not agree.

The target SLRP name alone is weak evidence because related SLRPs can be
clustered and annotations can inherit the wrong family label. Flanking genes
outside the SLRP family are more informative anchors.

### 3. Compare with SynVoy and protein evidence

Open the corresponding SynVoy SVG/HTML neighbourhood plot and compare the same
anchors. Then check the canonical protein accession in the final candidate
table and gene-specific tree. A locus with the right neighbours but an
ASPN-like protein, for example, is an annotation/product conflict rather than a
confident BGN ortholog.

### 4. If the gene is missing or poorly annotated

Do **not** record `gene absent` merely because the gene name is not displayed.
Use this escalation path:

1. inspect the raw GFF3 around the SynVoy region for `LOC`, predicted, partial,
   transcript, CDS, or pseudogene features;
2. search the expected protein accession and its linked Gene/assembly record;
3. run protein-to-genome `tblastn` against the exact target assembly, or use the
   GDV BLAST track, and check whether co-linear exon-like hits fall between the
   expected flanking genes;
4. inspect another, better assembly of the same species if available, but keep
   that evidence separate from the assembly actually analysed by SynVoy;
5. classify the result as `ambiguous/annotation-limited` if sequence or
   neighbourhood evidence is suggestive but no reliable model exists.

Annotation failure, assembly fragmentation, a gap, or a scaffold boundary can
all remove visible genes without representing biological gene loss.

## Decision vocabulary

- `accepted`: locus, protein, and tree are jointly compatible with the expected
  ortholog; small rearrangements or a generic `LOC` label are explained.
- `tentative`: most evidence supports the ortholog, but one material issue
  remains, such as a partial N terminus, missing SignalP, fragmented model, or
  weak neighbourhood.
- `rejected`: another paralog/product is better supported, or locus, protein,
  and tree evidence jointly contradict the assignment.
- `ambiguous`: available assembly/annotation is insufficient to decide.

Use `yes`, `partial`, `no`, or `not assessable` for
`neighbor_order_consistent`. Use `yes`, `no`, or `not available` for
`phylogeny_consistent`. Give one or two evidence-based sentences in
`manual_notes`; do not write only “looks good.”

## Highest-value existing cases

- **Chicken BGN:** the locus label and protein disagree; `XP_414298.2` is
  ASPN-like by reciprocal and tree evidence. Review flanking anchors, but keep
  the protein excluded from the BGN panel unless independent evidence resolves
  the conflict.
- **Amphioxus OGN:** already excluded from the confident OGN set; use the locus
  review to document why it is annotation-/orthology-limited, not to force a
  vertebrate OGN assignment.
- **Spotted-gar EPYC, zebrafish FMOD, and amphioxus/zebrafish PRELP:** check for
  fragmented or unnamed models and whether conserved flanks survive despite the
  missing retained SynVoy call.
- **Spotted-gar OGN:** structurally plausible and retained as tentative because
  SignalP is negative; synteny is supporting evidence, not a substitute for the
  N-terminal caveat.
- **Whale-shark EPYC:** record the independent alignment decision separately;
  its HIGH rescue-level synteny call does not settle the protein-level concern.

## Thesis wording

Report SynVoy as qualitative microsynteny support. Use “no retained SynVoy
candidate” or “annotation-limited locus” until gene loss has been independently
demonstrated. Do not convert HIGH/MEDIUM candidate counts into a count of proven
orthologs.
