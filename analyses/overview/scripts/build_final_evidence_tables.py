#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "analyses/overview/tables"
GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")

SPECIES_TO_SYNVOY = {
    "anolis_carolinensis": "anole",
    "bos_taurus": "cow",
    "callorhinchus_milii": "elephant_shark",
    "canis_lupus_familiaris": "dog",
    "chicken": "chicken",
    "frog": "frog",
    "latimeria_chalumnae": "coelacanth",
    "lepisosteus_oculatus": "spotted_gar",
    "monodelphis_domestica": "opossum",
    "mouse": "mouse",
    "scyliorhinus_canicula": "catshark",
    "whale_shark": "whale_shark",
    "zebrafish": "zebrafish",
    "amphioxus": "amphioxus",
}

LITERATURE = {
    "BGN": (
        "Mouse studies show that BGN and FMOD affect the skeletal matrix and "
        "osteoclast formation (Kram2017). BGN was used as the class-I reference "
        "in this study. MGI lists 14 skeletal, cartilage or joint phenotype "
        "terms for BGN, and HPO also contains human phenotype annotations."
    ),
    "DCN": (
        "DCN is a class-I SLRP involved in collagen fibril formation. It was "
        "included because it is the closest comparison to BGN within the selected "
        "genes. MGI and HPO contain skeletal and connective-tissue phenotypes for DCN."
    ),
    "FMOD": (
        "In developing mouse knees, FMOD was found around late-hypertrophic "
        "chondrocytes, the secondary ossification centre and the growth plate "
        "(Saamanen2001). The combined loss of BGN and FMOD also affects bone "
        "remodelling (Kram2017). MGI lists four skeletal, cartilage or joint terms."
    ),
    "PRELP": (
        "PRELP is a cartilage matrix protein with LRR and cysteine features. "
        "Previous studies connect it to osteoclast and osteoblast regulation "
        "(Grover1996; Grover2002; Rucci2009; Li2016PRELP). MGI contains phenotype "
        "annotations for PRELP, but no skeletal-category term in the used release."
    ),
    "EPYC": (
        "EPYC protein was found throughout the mouse growth-plate matrix around "
        "resting, proliferating and hypertrophic chondrocytes (Johnson1999). "
        "MGI also contains short-femur and osteoarthritis terms for EPYC."
    ),
    "LUM": (
        "LUM has been linked to skeletal matrix development and to collagen "
        "deposition and fibril formation in cartilage (Raouf2002; Kafienah2008). "
        "Most MGI terms concern the cornea and skin, but one tendon term is present."
    ),
    "OGN": (
        "OGN affects type-I collagen fibril formation after BMP1/Tolloid processing "
        "and has a known duplication history in vertebrates and teleosts "
        "(Ge2004; Costa2018). MGI contains ECM, corneal, skin and cardiometabolic "
        "phenotypes, but no skeletal-category term."
    ),
    "OMD": (
        "OMD is secreted by osteoblasts and can bind to mineral. It was found in "
        "the primary spongiosa of fetal growth plates and can increase BMP2/SMAD "
        "osteogenic signalling (Sommarin1998; Lin2021OMD)."
    ),
}

OVERALL = {
    "BGN": (
        "BGN was kept as a main class-I candidate. However, not every species could "
        "be resolved with the same confidence. Catshark remains tentative, the "
        "chicken record is ASPN-like, and elephant shark and opossum remain ambiguous."
    ),
    "DCN": (
        "DCN was kept as the full class-I comparison. The coelacanth, spotted-gar "
        "and zebrafish loci have conserved neighbouring genes. SignalP supports "
        "secretion for 14 of 15 proteins, while the opossum sequence lacks its N terminus."
    ),
    "FMOD": (
        "FMOD was kept as a main class-II candidate. The amphioxus assignment was "
        "rejected, catshark remains tentative, and zebrafish has two possible "
        "duplicated co-orthologs supported by different parts of the evidence."
    ),
    "PRELP": (
        "PRELP was kept as a main class-II candidate. Zebrafish was accepted, "
        "catshark remains tentative, and no clear gene-specific amphioxus ortholog was found."
    ),
    "EPYC": (
        "EPYC was kept as a class-III candidate with direct anatomical support from "
        "the literature. Most vertebrate loci were accepted, while zebrafish and "
        "whale shark remain tentative and the amphioxus SynVoy assignment was rejected."
    ),
    "LUM": (
        "LUM was kept as a vertebrate class-II candidate because the protein, tree "
        "and gene-structure results agree. Manual review supported the catshark, "
        "coelacanth and spotted-gar loci even when SynVoy selected another fragment. "
        "The amphioxus assignment remains unresolved."
    ),
    "OGN": (
        "OGN was kept as an exploratory vertebrate class-III candidate. Amphioxus "
        "was excluded from the confident set and spotted gar remains tentative. "
        "The Branchiostoma SynVoy result was not enough to identify an OGN ortholog."
    ),
    "OMD": (
        "OMD provides useful skeletal and growth-plate context. It was not included "
        "in the full protein, tree and synteny workflow and is therefore not treated "
        "as one of the main comparative genes."
    ),
}


def read_tsv(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def index(
    rows: list[dict[str, str]], *fields: str
) -> dict[tuple[str, ...], dict[str, str]]:
    return {tuple(row[field] for field in fields): row for row in rows}


def compact_number(value: str) -> str:
    if value == "":
        return ""
    try:
        return f"{float(value):.2f}".rstrip("0").rstrip(".")
    except ValueError:
        return value


def expression_text(row: dict[str, str]) -> str:
    return (
        f"In the GSE305415 WT P0 samples, the mean expression was "
        f"{compact_number(row['local_p0_mean_tpm'])} TPM "
        f"(rank {row['local_p0_rank_of_9']} of 9). In the growth-plate dataset, "
        f"the mean ranks were {compact_number(row['mouse_growth_plate_mean_rank'])} "
        f"in mouse and {compact_number(row['rat_growth_plate_mean_rank'])} in rat."
    )


def protein_text(row: dict[str, str]) -> str:
    return (
        f"The panel contained {row['Protein count']} proteins. "
        f"{row['SLRP-like domain count']} of {row['Protein count']} had SLRP-like "
        f"domains, and SignalP predicted a signal peptide for "
        f"{row['SignalP-positive count']} of {row['SignalP-tested count']} tested "
        f"proteins. The mean identity to the reference sequence was "
        f"{compact_number(row['Mean forward identity percent'])}%. Candidate-level "
        "motif results are listed in the separate candidate table."
    )


def msa_tree_text(
    gene: str,
    protein: dict[str, str],
    msa: dict[str, str],
    tree: dict[str, str],
    codon: dict[str, str],
) -> str:
    outside = tree["tips_outside_largest_pure_gene_split"].strip()
    tree_part = (
        f"The largest gene-specific split contained "
        f"{tree['largest_pure_gene_split_tip_count']} of {tree['tip_count']} sequences"
    )
    if tree["unrooted_monophyletic"] == "yes":
        tree_part += ", which formed one group in the unrooted tree"
    elif outside:
        tree_part += f". The sequences outside this split were: {outside}"
    split_label = tree.get("largest_pure_gene_split_label", "").strip()
    if split_label:
        tree_part += f". Its SH-aLRT/UFBoot support was {split_label}"
    return (
        f"The mean pairwise identity in the alignment was "
        f"{compact_number(msa['mean_pairwise_identity_percent'])}%. {tree_part}. "
        f"All {codon['interpretable_omega_count']} interpretable human-target NG86 "
        f"ratios were below 1, with a median of "
        f"{compact_number(codon['median_dN_dS_omega'])}."
    )


def synvoy_gene_text(
    gene: str,
    summary: dict[str, dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> str:
    row = summary.get(gene)
    if row is None:
        return (
            f"Not available: no completed 14-target SynVoy report was found for {gene}."
        )
    reviewed = [
        evidence
        for evidence in evidence_rows
        if evidence["gene"] == gene and evidence["review_status"]
    ]
    counts = Counter(evidence["review_status"] for evidence in reviewed)
    manual_part = ""
    if reviewed:
        ordered_statuses = ("accepted", "tentative", "ambiguous", "rejected")
        status_text = ", ".join(
            f"{status}={counts[status]}"
            for status in ordered_statuses
            if counts[status]
        )
        manual_part = f" Manual review gave the following results: {status_text}."
    return (
        f"SynVoy processed 14 target genomes. After filtering, "
        f"{row['final_high_goi']} HIGH and {row['final_medium_goi']} MEDIUM "
        f"candidate calls remained.{manual_part} The matched-input opossum row "
        "was included after the same locus-level review used for the other species."
    )


def structure_text(row: dict[str, str], extended: dict[str, str]) -> str:
    text = (
        f"The gene-structure comparison included five reference species. The CDS "
        f"contained {row['cds_exon_count_range']} exons and the estimated proteins "
        f"were {row['protein_length_min_aa']}-{row['protein_length_max_aa']} aa long. "
        f"There were {row['structural_outlier_count']} representative structural "
        f"outliers. Intron phase was conserved at "
        f"{extended['phase_conserved_junction_count']} of "
        f"{extended['splice_junction_count']} comparable splice junctions, and all "
        f"{extended['complete_start_count']} of {extended['species_count']} CDS "
        f"models passed the completeness check."
    )
    if row["gene"] == "BGN":
        text += (
            " The annotated chicken locus does not independently support a BGN "
            "ortholog because its protein is ASPN-like. The zebrafish structure "
            "uses bgna, while the protein panel uses bgnb."
        )
    return text


def build_gene_table() -> None:
    candidate = {
        r["gene"]: r for r in read_tsv("analyses/overview/candidate_gene_evidence.tsv")
    }
    protein = {
        r["Gene"]: r
        for r in read_tsv(
            "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_gene_summary.tsv"
        )
    }
    msa = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/protein_analysis/tables/msa_conservation_gene_summary.tsv"
        )
    }
    tree = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/phylogenetics/combined_trees/seven_gene_tree/compact_tree_gene_clade_review.tsv"
        )
    }
    structure = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/gene_structure/tables/gene_structure_extended_summary.tsv"
        )
    }
    extended_structure = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/gene_structure/extended/tables/extended_gene_structure_gene_summary.tsv"
        )
    }
    codon = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/evolutionary_rates/tables/codon_constraint_gene_summary.tsv"
        )
    }
    synvoy = {
        r["gene"]: r for r in read_tsv("analyses/synteny/tables/synvoy_run_summary.tsv")
    }
    synvoy_evidence = read_tsv(
        "analyses/synteny/tables/synvoy_gene_species_evidence.tsv"
    )
    manual_issues = {
        "BGN": (
            "Chicken XP_414298.2 was excluded because the protein is ASPN-like. "
            "Earlier opossum, elephant-shark and catshark candidates were removed "
            "because they were DCN sequences. Catshark XP_038638277.1 lacks part of "
            "its C terminus, and whale-shark XP_048475905.1 is a discontinuous fragment "
            "without a detected signal peptide. Both were therefore kept as tentative records."
        ),
        "DCN": (
            "Opossum XP_001363160.3 only contains the conserved C-terminal part of "
            "DCN and lacks the N-terminal signal peptide. Its matched-assembly locus "
            "is supported by the conserved EPYC--KERA--LUM--DCN neighbourhood. "
            "Catshark XP_038636759.1 was kept, "
            "but its N-terminal extension and later SignalP cleavage site were noted."
        ),
        "FMOD": (
            "No sequence-level decision remains unresolved. Predicted proteins and "
            "sequences from early-diverging lineages still have the usual annotation limitations."
        ),
        "PRELP": (
            "No sequence-level decision remains unresolved. The position of the lamprey "
            "sequence should be interpreted carefully because the tree is unrooted."
        ),
        "EPYC": (
            "Whale-shark XP_048463897.1 was kept as tentative because its earlier "
            "sequence region is less consistent and its query coverage is reduced."
        ),
        "LUM": (
            "Lamprey XP_075930353.1 lies outside the main LUM split and was therefore "
            "kept as a tentative assignment from an early-diverging lineage."
        ),
        "OGN": (
            "The opossum and zebrafish sequences were kept. Spotted gar remains "
            "tentative because SignalP did not predict a cleavage site, although its "
            "N terminus is strongly hydrophobic. Amphioxus was excluded from the "
            "confident OGN set and kept only as an uncertain SLRP-like sequence."
        ),
    }

    rows: list[dict[str, str]] = []
    for gene in GENES:
        rows.append(
            {
                "Gene": gene,
                "expression support": expression_text(candidate[gene]),
                "literature/function support": LITERATURE[gene],
                "protein/domain conservation": protein_text(protein[gene]),
                "MSA/tree support": msa_tree_text(
                    gene, protein[gene], msa[gene], tree[gene], codon[gene]
                ),
                "SynVoy synteny support": synvoy_gene_text(
                    gene, synvoy, synvoy_evidence
                ),
                "gene-structure support": structure_text(
                    structure[gene], extended_structure[gene]
                ),
                "manual-review issues": manual_issues[gene],
                "overall interpretation": OVERALL[gene],
            }
        )

    omd = candidate["OMD"]
    rows.append(
        {
            "Gene": "OMD",
            "expression support": expression_text(omd),
            "literature/function support": LITERATURE["OMD"],
            "protein/domain conservation": (
                "OMD was not included in the seven-gene protein, domain and SignalP "
                "analysis. However, the translated sequence from each of the five "
                "reference species had LRR-family Pfam support."
            ),
            "MSA/tree support": (
                "OMD was not included in the seven-gene alignment or tree. All "
                f"{codon['OMD']['interpretable_omega_count']} interpretable human-target "
                f"NG86 ratios were below 1, with a median of "
                f"{compact_number(codon['OMD']['median_dN_dS_omega'])}."
            ),
            "SynVoy synteny support": (
                "SynVoy was not run for OMD because it was included only as a context gene."
            ),
            "gene-structure support": structure_text(
                structure["OMD"], extended_structure["OMD"]
            ),
            "manual-review issues": (
                "No candidate-level manual review table was made because OMD was not "
                "included in the full protein workflow."
            ),
            "overall interpretation": OVERALL["OMD"],
        }
    )

    fields = [
        "Gene",
        "expression support",
        "literature/function support",
        "protein/domain conservation",
        "MSA/tree support",
        "SynVoy synteny support",
        "gene-structure support",
        "manual-review issues",
        "overall interpretation",
    ]
    write_tsv(OUT_DIR / "final_gene_level_evidence.tsv", rows, fields)


def synteny_text(
    gene: str, species: str, synteny: dict[tuple[str, str], dict[str, str]]
) -> str:
    if species == "human":
        return "Human reference locus."
    alias = SPECIES_TO_SYNVOY.get(species)
    if alias is None:
        return "This species was not included in the 14-target SynVoy analysis."
    row = synteny.get((gene, alias))
    if row is None:
        return "No completed SynVoy result was available for this gene."
    updated_confidence = {
        ("FMOD", "catshark"): "NONE",
        ("FMOD", "zebrafish"): "MEDIUM",
        ("FMOD", "whale_shark"): "MEDIUM",
    }.get((gene, alias))
    if row["review_status"]:
        confirmed = row["confirmed_gene_symbol"]
        accession = row["confirmed_protein_accession"]
        confirmed_part = ""
        if confirmed and confirmed != "not found":
            confirmed_part = f"; confirmed {confirmed}"
            if accession and accession != "not available":
                confirmed_part += f" ({accession})"
        synvoy_result = (
            f"{updated_confidence} in the updated SynVoy report"
            if updated_confidence
            else f"{row['best_confidence']} (grade {row['evidence_grade']})"
        )
        return (
            f"Manual review: {row['review_status']}; position: "
            f"{row['coordinate_to_accession_status']}; expected neighbour order: "
            f"{row['neighbor_order_consistent']}{confirmed_part}; SynVoy result: "
            f"{synvoy_result}."
        )
    return (
        f"SynVoy result: {row['best_confidence']} (grade "
        f"{row['evidence_grade']}); review priority: {row['review_priority']}."
    )


def tree_text(gene: str, accession: str, tree: dict[str, str]) -> str:
    outside = tree["tips_outside_largest_pure_gene_split"]
    if accession and accession in outside:
        return (
            f"Outside the largest {gene}-only split, which contained "
            f"{tree['largest_pure_gene_split_tip_count']} of {tree['tip_count']} sequences."
        )
    label = tree.get("largest_pure_gene_split_label", "").strip()
    if not label:
        label = tree["supporting_split_label"].strip()
    support = f" The SH-aLRT/UFBoot support was {label}." if label else ""
    return (
        f"Inside the largest {gene}-only split, which contained "
        f"{tree['largest_pure_gene_split_tip_count']} of {tree['tip_count']} "
        f"sequences.{support}"
    )


def reviewed_synteny_note(row: dict[str, str]) -> str:
    """Use short wording for routine matches and keep details for difficult loci."""
    if (
        row["coordinate_to_accession_status"] == "same locus"
        and row["neighbor_order_consistent"] == "yes"
    ):
        symbol = row["confirmed_gene_symbol"]
        locus = f"the annotated {symbol} locus" if symbol else "the same annotated locus"
        if row["review_status"] == "accepted":
            return (
                f"The selected protein and the SynVoy candidate match {locus}. "
                "The neighbouring genes also remain in the expected order."
            )
        if row["review_status"] == "tentative":
            return (
                f"The selected protein and the SynVoy candidate match {locus}, and "
                "the neighbouring genes remain in the expected order. Other protein "
                "or annotation concerns keep the assignment tentative."
            )
    return row["manual_notes"]


def decision_text(value: str) -> str:
    replacements = {
        "retain": "retained",
        "retain as tentative": "retained as tentative",
        "exclude from confident OGN ortholog set": (
            "excluded from the confident OGN ortholog set"
        ),
        "retain as confident DCN ortholog with N-terminal model caveat": (
            "retained as a confident DCN ortholog, with an N-terminal model caveat"
        ),
        "retain as partial DCN sequence; locus-level orthology unresolved": (
            "retained as a partial DCN sequence; the locus assignment remains unresolved"
        ),
        "retain as partial DCN sequence; locus-level orthology supported": (
            "retained as a partial DCN sequence; the locus assignment is supported"
        ),
        "retain as tentative BGN locus/fragment; exclude from complete-protein claims": (
            "retained as a tentative BGN locus and fragment; not used as a complete protein"
        ),
        "retain as tentative partial BGN ortholog": (
            "retained as a tentative partial BGN ortholog"
        ),
    }
    return replacements.get(value, value)


def domain_text(quality: str, hit_count: str) -> str:
    if quality == "Good SLRP-like domain architecture":
        return f"Expected SLRP-like domain pattern with {hit_count} LRR-related hits."
    return f"Needs inspection; {hit_count} LRR-related hits were found."


def signalp_text(prediction: str, quality: str) -> str:
    descriptions = {
        "strong_signal_peptide": "strong signal-peptide prediction",
        "weak_signal_peptide": "weak signal-peptide prediction",
        "low_cleavage_probability": "low cleavage-site probability",
        "watch_cleavage_probability": "cleavage-site probability should be checked",
        "no_signal_peptide": "no signal peptide predicted",
    }
    if prediction != "SP" and quality == "no_signal_peptide":
        return "No signal peptide was predicted."
    result = "Signal peptide predicted" if prediction == "SP" else "No signal peptide predicted"
    return f"{result}; {descriptions.get(quality, quality.replace('_', ' '))}."


def qc_text(value: str) -> str:
    if not value:
        return ""
    descriptions = {
        "high_alignment_gap": "large gap in the alignment",
        "no_signalp_prediction": "no SignalP prediction",
        "low_query_coverage": "low query coverage",
        "reciprocal_hit_differs_from_reference": (
            "reciprocal hit differs from the expected reference"
        ),
        "reduced_query_coverage": "reduced query coverage",
        "signalp_watch": "SignalP result should be checked",
        "weak_domain_support": "weak domain support",
    }
    items = [descriptions.get(item, item.replace("_", " ")) for item in value.split(";")]
    return "QC points: " + "; ".join(items) + "."


def evidence_text(value: str) -> str:
    return (
        value.replace("support BGN provenance", "support a BGN assignment")
        .replace("Retain only as", "It was kept only as")
        .replace("provenance candidate", "record")
    )


def build_candidate_table() -> None:
    proteins = read_tsv(
        "analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv"
    )
    manual = {
        r["accession"]: r
        for r in read_tsv(
            "analyses/protein_analysis/manual_review/manual_sequence_review_decisions.tsv"
        )
    }
    old_manual = {
        r["Protein accession"]: r
        for r in read_tsv(
            "analyses/protein_analysis/manual_review/manual_sequence_review_5.tsv"
        )
    }
    trees = {
        r["gene"]: r
        for r in read_tsv(
            "analyses/phylogenetics/combined_trees/seven_gene_tree/compact_tree_gene_clade_review.tsv"
        )
    }
    motif_assignments = index(
        read_tsv(
            "analyses/protein_analysis/motifs/tables/"
            "protein_candidate_motif_model_assignments.tsv"
        ),
        "gene",
        "species",
        "accession",
    )
    synteny = index(
        read_tsv("analyses/synteny/tables/synvoy_gene_species_evidence.tsv"),
        "gene",
        "species",
    )
    locus_decisions = index(
        read_tsv(
            "analyses/synteny/diagnostics/priority_locus_reviews/"
            "all_synvoy_review_decisions.tsv"
        ),
        "gene",
        "species",
    )
    for key, decision in locus_decisions.items():
        if key not in synteny:
            continue
        for field in (
            "coordinate_to_accession_status",
            "review_status",
            "confirmed_gene_symbol",
            "confirmed_protein_accession",
            "neighbor_order_consistent",
            "phylogeny_consistent",
            "reviewer",
            "review_date",
            "manual_notes",
        ):
            synteny[key][field] = decision[field]

    rows: list[dict[str, str]] = []
    for p in proteins:
        gene, species, accession = p["Gene"], p["Species"], p["Protein accession"]
        manual_row = manual.get(accession)
        motif = motif_assignments[(gene, species, accession)]
        motif_status = (
            "the model from the same gene scored highest"
            if motif["motif_model_status"] == "own gene model best"
            else motif["motif_model_status"]
        )
        manual_status = f"MSA: {p['MSA quality']}; motif comparison: {motif_status}"
        final_decision = "retained as a confident ortholog"
        decision_note = ""
        if manual_row:
            plain_decision = decision_text(manual_row["manual_decision"])
            manual_status += f"; manual decision: {plain_decision}"
            final_decision = plain_decision
            decision_note = evidence_text(manual_row["manual_evidence"])
        if gene == "LUM" and species == "lamprey":
            final_decision = "retained as a tentative early-lineage LUM assignment"
            decision_note = (
                "The sequence lies outside the main LUM split, but its annotation "
                "and reciprocal comparison still support a lumican-like assignment."
            )

        alias = SPECIES_TO_SYNVOY.get(species)
        synteny_row = synteny.get((gene, alias)) if alias else None
        if synteny_row and synteny_row["review_status"]:
            locus_status = synteny_row["review_status"]
            if (
                locus_status == "tentative"
                and final_decision == "retained as a confident ortholog"
            ):
                final_decision = "retained as a tentative ortholog or co-ortholog"
            elif (
                locus_status == "ambiguous"
                and final_decision == "retained as a confident ortholog"
            ):
                final_decision = (
                    "protein retained provisionally; locus assignment remains ambiguous"
                )
            elif (
                locus_status == "rejected"
                and final_decision == "retained as a confident ortholog"
            ):
                final_decision = "excluded because manual locus review rejected the assignment"
            decision_note = " ".join(
                x
                for x in [
                    decision_note,
                    reviewed_synteny_note(synteny_row),
                ]
                if x
            )

        concerns = " ".join(x for x in [qc_text(p["QC concerns"]), decision_note] if x)
        rows.append(
            {
                "Gene": gene,
                "species": species,
                "accession": accession,
                "protein length": p["Analyzed protein length"],
                "domain status": domain_text(
                    p["Domain quality"], p["Number of LRR-related domains"]
                ),
                "SignalP status": signalp_text(
                    p["SignalP prediction"], p["SignalP quality"]
                ),
                "MSA/manual status": manual_status,
                "tree support": tree_text(gene, accession, trees[gene]),
                "synteny support": synteny_text(gene, species, synteny),
                "final decision": final_decision,
                "notes": concerns,
            }
        )

    p = old_manual["XP_066265713.1"]
    m = manual["XP_066265713.1"]
    rows.append(
        {
            "Gene": "OGN",
            "species": "amphioxus",
            "accession": "XP_066265713.1",
            "protein length": p["Analyzed protein length"],
            "domain status": (
                f"{domain_text(p['Domain quality'], p['Number of LRR-related domains'])} "
                "The pattern is SLRP-like but does not specifically identify OGN."
            ),
            "SignalP status": signalp_text(
                p["SignalP prediction"], p["SignalP quality"]
            ),
            "MSA/manual status": (
                f"Historical MSA: {p['MSA quality']}; manual decision: "
                f"{decision_text(m['manual_decision'])}"
            ),
            "tree support": (
                "The sequence was outside the historical OGN-only split containing "
                "14 of 15 sequences and was not included in the final 103-sequence tree."
            ),
            "synteny support": synteny_text("OGN", "amphioxus", synteny),
            "final decision": (
                "excluded from the confident OGN set; kept as an uncertain SLRP-like record"
            ),
            "notes": (
                "The sequence is SLRP-like, but its terminal structure and several "
                "large insertions differ from the vertebrate OGN sequences. The tree "
                "and synteny results did not support a clear OGN assignment, so it "
                "was kept only as an uncertain SLRP-like record."
            ),
        }
    )

    rows.append(
        {
            "Gene": "FMOD",
            "species": "zebrafish",
            "accession": "NP_001025243.1",
            "protein length": "342",
            "domain status": "Not included in the canonical domain analysis.",
            "SignalP status": "Not included in the canonical SignalP analysis.",
            "MSA/manual status": (
                "Additional duplicated co-ortholog; not yet integrated into the canonical panel."
            ),
            "tree support": "Not included in the final FMOD tree.",
            "synteny support": (
                "Manual review: tentative. fmoda lies next to prelp in the conserved "
                "ancestral FMOD-PRELP region."
            ),
            "final decision": (
                "retained as an additional teleost FMOD co-ortholog; protein QC is "
                "needed before inclusion in the canonical panel"
            ),
            "notes": (
                "The BLASTP comparison to human FMOD gave 56.6% identity across 80% "
                "of the query. The selected fmodb protein has broader sequence "
                "coverage, while fmoda has stronger support from the ancestral locus."
            ),
        }
    )

    rows.append(
        {
            "Gene": "FMOD",
            "species": "callorhinchus_milii",
            "accession": "XP_007897806.2",
            "protein length": "684 (exact GFF/FASTA reconstruction; compound model)",
            "domain status": (
                "Both halves contain several Pfam LRR hits and significant LRRNT "
                "regions at aa 54-76 and 383-412."
            ),
            "SignalP status": (
                "SignalP was not run. The model has a hydrophobic N terminus and a "
                "second signal-peptide-like region at aa 330-347."
            ),
            "MSA/manual status": (
                "Manual and BLASTP checks found a DCN/BGN-like N-terminal half and "
                "an FMOD-like C-terminal half."
            ),
            "tree support": "Not included in the final FMOD tree.",
            "synteny support": (
                "Manual review: tentative. The SynVoy candidate overlaps LOC103182549. "
                "Five human neighbouring genes are shared, and four keep their order "
                "after accounting for the reverse orientation."
            ),
            "final decision": (
                "locus retained as tentative FMOD evidence; compound protein model "
                "excluded from the confident FMOD panel"
            ),
            "notes": (
                "The model contains nine CDS exons and a complete 2055-bp CDS, but "
                "its two SLRP-like halves suggest that two models may have been joined. "
                "Residues 1-330 matched DCN best at 47.0% identity, while residues "
                "376-684 matched FMOD best at 57.0% identity and covered 82.2% of "
                "human FMOD. Transcript evidence or a revised model is needed before "
                "this protein can be included as an ortholog."
            ),
        }
    )

    rows.append(
        {
            "Gene": "BGN",
            "species": "chicken",
            "accession": "XP_414298.2",
            "protein length": "371",
            "domain status": (
                "Good SLRP-like architecture (7 LRR-related hits), but the domains "
                "do not distinguish BGN from related SLRPs."
            ),
            "SignalP status": (
                "A strong signal peptide was predicted in the earlier excluded-panel run."
            ),
            "MSA/manual status": (
                "Historical MSA: Good; the annotation and reciprocal result disagree."
            ),
            "tree support": (
                "In the focused reference tree, the sequence was sister to human ASPN "
                "with 99.3 SH-aLRT and 100 UFBoot support."
            ),
            "synteny support": synteny_text("BGN", "chicken", synteny),
            "final decision": "excluded from the canonical BGN protein panel",
            "notes": (
                "The sequence was kept only to document the annotation conflict. "
                "It should not be described as chicken BGN unless separate locus "
                "evidence supports that assignment."
            ),
        }
    )

    rows.sort(key=lambda r: (GENES.index(r["Gene"]), r["species"], r["accession"]))
    fields = [
        "Gene",
        "species",
        "accession",
        "protein length",
        "domain status",
        "SignalP status",
        "MSA/manual status",
        "tree support",
        "synteny support",
        "final decision",
        "notes",
    ]
    write_tsv(OUT_DIR / "final_candidate_level_confidence.tsv", rows, fields)


def main() -> None:
    build_gene_table()
    build_candidate_table()
    print("Built final gene-level and candidate-level evidence tables.")


if __name__ == "__main__":
    main()
