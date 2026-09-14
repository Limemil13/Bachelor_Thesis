#!/usr/bin/env python3
"""Build thesis-facing gene- and candidate-level evidence tables.

The canonical protein table is the source of retained candidates.  Explicitly
excluded provenance records are appended so that curation decisions are not
lost when the confident panel is regenerated.
"""

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
        "BGN and FMOD regulate skeletal matrix and osteoclastogenesis in mouse "
        "genetics/biochemistry (Kram2017); BGN is used as the class-I anchor. "
        "Current MGI adds 14 skeletal/cartilage/joint MP terms, and pinned HPO "
        "contains BGN human disease-phenotype annotations."
    ),
    "DCN": (
        "DCN is a canonical class-I SLRP and collagen-fibrillogenesis regulator "
        "with broad cartilage and connective-tissue relevance. It provides the "
        "closest functional and evolutionary comparator for BGN; current MGI/HPO "
        "evidence includes skeletal and connective-tissue phenotypes."
    ),
    "FMOD": (
        "Developmental mouse data localize FMOD around late-hypertrophic "
        "chondrocytes, the secondary ossification centre and growth plate "
        "(Saamanen2001); BGN/FMOD double deficiency alters bone remodelling "
        "(Kram2017). Current MGI adds four skeletal/cartilage/joint MP terms."
    ),
    "PRELP": (
        "PRELP is a cartilage matrix protein with conserved LRR/cysteine "
        "features and reported roles in osteoclast and osteoblast regulation "
        "(Grover1996; Grover2002; Rucci2009; Li2016PRELP). Current MGI contains "
        "PRELP phenotype annotations but no skeletal-category MP term in this release."
    ),
    "EPYC": (
        "Developmental mouse data place EPYC protein throughout growth-plate "
        "ECM around resting, proliferating and hypertrophic chondrocytes "
        "(Johnson1999). Current MGI adds short-femur and osteoarthritis terms."
    ),
    "LUM": (
        "LUM is linked to skeletal matrix development and cartilage collagen "
        "deposition/fibrillogenesis (Raouf2002; Kafienah2008). Current MGI is "
        "dominated by corneal/skin phenotypes and includes a tendon term."
    ),
    "OGN": (
        "OGN regulates type-I collagen fibrillogenesis after BMP1/Tolloid "
        "processing and has a documented vertebrate/teleost duplication history "
        "(Ge2004; Costa2018). Current MGI supports broader ECM/corneal/skin and "
        "cardiometabolic phenotypes but has no skeletal-category term."
    ),
    "OMD": (
        "OMD is an osteoblast-secreted mineral-binding protein localized to "
        "fetal growth-plate primary spongiosa and can enhance BMP2/SMAD "
        "osteogenic signalling (Sommarin1998; Lin2021OMD)."
    ),
}

OVERALL = {
    "BGN": (
        "Strong class-I reference candidate, but the cross-species locus set is not uniformly resolved: "
        "catshark is tentative, the chicken record is ASPN-like, and elephant shark/opossum remain ambiguous."
    ),
    "DCN": (
        "Retain as a full class-I comparator. Canonical coelacanth, spotted-gar and zebrafish "
        "loci have conserved neighbourhoods, and the completed SignalP screen supports secretion "
        "for 14 of 15 proteins; the negative opossum record is an N-terminally incomplete fragment."
    ),
    "FMOD": (
        "Strong class-II candidate. Amphioxus is rejected, catshark remains tentative, and zebrafish "
        "contains two plausible duplicated co-orthologs with complementary sequence and synteny support."
    ),
    "PRELP": (
        "Strong class-II growth-plate candidate; zebrafish is accepted, catshark is tentative, "
        "and a gene-specific amphioxus ortholog is unresolved."
    ),
    "EPYC": (
        "Biologically specific class-III candidate with strong anatomical literature support. "
        "Several vertebrate loci are accepted, zebrafish and whale shark remain tentative, and the amphioxus SynVoy call is rejected."
    ),
    "LUM": (
        "Strong vertebrate class-II candidate with clean protein, phylogenetic and gene-structure support. "
        "Manual review accepted the canonical catshark, coelacanth and spotted-gar loci despite incorrect SynVoy fragment selection; amphioxus is unresolved."
    ),
    "OGN": (
        "Retain as an exploratory vertebrate class-III candidate; the confident set excludes amphioxus "
        "and keeps spotted gar tentative. The Branchiostoma SynVoy evidence does not establish an OGN ortholog."
    ),
    "OMD": "Useful skeletal/growth-plate context gene, but not part of the full canonical protein/tree/synteny panel and therefore not a co-equal main candidate.",
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
        f"Local P0 {compact_number(row['local_p0_mean_tpm'])} TPM "
        f"(rank {row['local_p0_rank_of_9']}/9); mouse/rat growth-plate mean "
        f"ranks {compact_number(row['mouse_growth_plate_mean_rank'])}/"
        f"{compact_number(row['rat_growth_plate_mean_rank'])}."
    )


def protein_text(row: dict[str, str]) -> str:
    return (
        f"{row['Protein count']} proteins; {row['SLRP-like domain count']}/"
        f"{row['Protein count']} SLRP-like; {row['SignalP-positive count']}/"
        f"{row['SignalP-tested count']} tested SignalP-positive; "
        f"{row['SignalP-pending count']} SignalP pending; mean forward identity "
        f"{compact_number(row['Mean forward identity percent'])}%. "
        "Protein-motif support is reported separately at candidate level."
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
        f"{tree['largest_pure_gene_split_tip_count']}/{tree['tip_count']} tips "
        f"in largest pure split"
    )
    if tree["unrooted_monophyletic"] == "yes":
        tree_part += "; unrooted monophyletic"
    elif outside:
        tree_part += f"; outside: {outside}"
    return (
        f"Mean pairwise MSA identity "
        f"{compact_number(msa['mean_pairwise_identity_percent'])}%; {tree_part}; "
        f"all {codon['interpretable_omega_count']} finite human-target NG86 "
        f"ratios below 1 (median {compact_number(codon['median_dN_dS_omega'])})."
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
        manual_part = f" All target reviews recorded ({status_text})."
    return (
        f"14 targets processed; {row['final_high_goi']} HIGH and "
        f"{row['final_medium_goi']} MEDIUM retained gene-of-interest calls; "
        f"{manual_part} Opossum output is excluded pending a rerun with matching FASTA/GFF sequence IDs."
    )


def structure_text(row: dict[str, str], extended: dict[str, str]) -> str:
    return (
        f"Five reference species; CDS-exon count {row['cds_exon_count_range']}; "
        f"estimated protein {row['protein_length_min_aa']}-"
        f"{row['protein_length_max_aa']} aa; "
        f"{row['structural_outlier_count']} representative outliers; "
        f"{extended['phase_conserved_junction_count']}/"
        f"{extended['splice_junction_count']} comparable splice junctions "
        f"preserve intron phase; complete CDS QC in "
        f"{extended['complete_start_count']}/{extended['species_count']} species."
    )


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
            "analyses/phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv"
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
    manual = {
        r["accession"]: r
        for r in read_tsv(
            "analyses/protein_analysis/manual_review/manual_sequence_review_decisions.tsv"
        )
    }

    manual_issues = {
        "BGN": (
            "Chicken XP_414298.2 is excluded as ASPN-like; historical opossum, "
            "elephant-shark and catshark entries were DCN contaminants. Catshark "
            "XP_038638277.1 has a strong signal peptide but is C-terminally incomplete; "
            "whale-shark XP_048475905.1 is a discontinuous fragment without a detectable "
            "signal peptide. Both remain tentative BGN locus/protein records."
        ),
        "DCN": (
            "Opossum XP_001363160.3 is a conserved C-terminal fragment without the "
            "N-terminal signal peptide and its locus cannot be assessed with the mismatched "
            "SynVoy input pair. Catshark XP_038636759.1 is retained confidently, with its "
            "N-terminal extension and later SignalP cleavage recorded as a model caveat."
        ),
        "FMOD": "No unresolved sequence-level decision; deep-lineage and predicted records retain normal annotation caveats.",
        "PRELP": "No unresolved sequence-level decision; lamprey/deep-lineage placement should be described without overinterpreting root direction.",
        "EPYC": f"Whale shark XP_048463897.1: {manual['XP_048463897.1']['manual_decision']}.",
        "LUM": "Lamprey XP_075930353.1 is outside the main LUM split and is retained only as a tentative deep-lineage assignment.",
        "OGN": (
            "Opossum and zebrafish retained; spotted gar retained as tentative because SignalP is negative; "
            "amphioxus excluded from the confident OGN set and retained only as uncertain SLRP-like provenance."
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
                "Not run in the canonical seven-gene protein/domain/SignalP panel; "
                "all five representative translations have LRR-family Pfam support "
                "in the extended structure/domain map."
            ),
            "MSA/tree support": (
                "Not run in the canonical seven-gene MSA/phylogeny panel; all "
                f"{codon['OMD']['interpretable_omega_count']} finite human-target NG86 "
                f"ratios are below 1 (median {compact_number(codon['OMD']['median_dN_dS_omega'])})."
            ),
            "SynVoy synteny support": "Not run; OMD is a context gene rather than a full comparative-panel member.",
            "gene-structure support": structure_text(
                structure["OMD"], extended_structure["OMD"]
            ),
            "manual-review issues": "No candidate-level manual packet because OMD was not advanced to the full protein pipeline.",
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
        return "human reference locus"
    alias = SPECIES_TO_SYNVOY.get(species)
    if alias is None:
        return "not analyzed in 14-target SynVoy panel"
    row = synteny.get((gene, alias))
    if row is None:
        return "pending (gene not yet processed by SynVoy)"
    if alias == "opossum":
        return (
            "not interpretable: the opossum target FASTA/GFF sequence IDs do not "
            "match; exclude this SynVoy result until the target is rerun"
        )
    if row["review_status"]:
        confirmed = row["confirmed_gene_symbol"]
        accession = row["confirmed_protein_accession"]
        confirmed_part = ""
        if confirmed and confirmed != "not found":
            confirmed_part = f"; confirmed {confirmed}"
            if accession and accession != "not available":
                confirmed_part += f" ({accession})"
        return (
            f"manual {row['review_status']}; coordinate relation "
            f"{row['coordinate_to_accession_status']}; neighbour order "
            f"{row['neighbor_order_consistent']}{confirmed_part}; automated "
            f"SynVoy {row['best_confidence']} (grade {row['evidence_grade']})"
        )
    return (
        f"{row['best_confidence']}; evidence grade {row['evidence_grade']}; "
        f"review {row['review_priority']}"
    )


def tree_text(gene: str, accession: str, tree: dict[str, str]) -> str:
    outside = tree["tips_outside_largest_pure_gene_split"]
    if accession and accession in outside:
        return (
            f"outside largest pure {gene} split "
            f"({tree['largest_pure_gene_split_tip_count']}/{tree['tip_count']})"
        )
    label = tree["supporting_split_label"].strip()
    support = f"; split label {label}" if label else ""
    return (
        f"inside largest pure {gene} split "
        f"({tree['largest_pure_gene_split_tip_count']}/{tree['tip_count']}){support}"
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
            "analyses/phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv"
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

    rows: list[dict[str, str]] = []
    for p in proteins:
        gene, species, accession = p["Gene"], p["Species"], p["Protein accession"]
        manual_row = manual.get(accession)
        motif = motif_assignments[(gene, species, accession)]
        manual_status = f"MSA {p['MSA quality']}; motif: {motif['motif_model_status']}"
        final_decision = "retain in confident ortholog panel"
        decision_note = ""
        if manual_row:
            manual_status += f"; manual: {manual_row['manual_decision']}"
            final_decision = manual_row["manual_decision"]
            decision_note = manual_row["manual_evidence"]
        if gene == "LUM" and species == "lamprey":
            final_decision = "retain as tentative deep-lineage LUM assignment"
            decision_note = "Outside the main LUM split; annotation and reciprocal evidence still support a lumican-like assignment."

        alias = SPECIES_TO_SYNVOY.get(species)
        synteny_row = synteny.get((gene, alias)) if alias else None
        if synteny_row and synteny_row["review_status"]:
            locus_status = synteny_row["review_status"]
            if locus_status == "tentative" and final_decision == "retain in confident ortholog panel":
                final_decision = "retain as tentative ortholog/co-ortholog assignment"
            elif locus_status == "ambiguous" and final_decision == "retain in confident ortholog panel":
                final_decision = "retain protein provisionally; locus-level orthology remains ambiguous"
            elif locus_status == "rejected" and final_decision == "retain in confident ortholog panel":
                final_decision = "exclude: manual locus review rejected the gene assignment"
            decision_note = "; ".join(
                x
                for x in [
                    decision_note,
                    f"Manual synteny review: {synteny_row['manual_notes']}",
                ]
                if x
            )

        concerns = "; ".join(
            x for x in [p["QC concerns"], p["Notes"], decision_note] if x
        )
        rows.append(
            {
                "Gene": gene,
                "species": species,
                "accession": accession,
                "protein length": p["Analyzed protein length"],
                "domain status": f"{p['Domain quality']} ({p['Number of LRR-related domains']} LRR-related hits)",
                "SignalP status": f"{p['SignalP prediction']}; {p['SignalP quality']}",
                "MSA/manual status": manual_status,
                "tree support": tree_text(gene, accession, trees[gene]),
                "synteny support": synteny_text(gene, species, synteny),
                "final decision": final_decision,
                "notes": concerns,
            }
        )

    # Preserve the explicitly excluded amphioxus OGN as an auditable provenance row.
    p = old_manual["XP_066265713.1"]
    m = manual["XP_066265713.1"]
    rows.append(
        {
            "Gene": "OGN",
            "species": "amphioxus",
            "accession": "XP_066265713.1",
            "protein length": p["Analyzed protein length"],
            "domain status": f"{p['Domain quality']} ({p['Number of LRR-related domains']} LRR-related hits); SLRP-like but not OGN-specific",
            "SignalP status": f"{p['SignalP prediction']}; {p['SignalP quality']}",
            "MSA/manual status": f"historical MSA {p['MSA quality']}; manual: {m['manual_decision']}",
            "tree support": "outside the historical 14/15 pure OGN split; omitted from the final 103-tip tree",
            "synteny support": synteny_text("OGN", "amphioxus", synteny),
            "final decision": "exclude from confident OGN ortholog set; retain as uncertain SLRP-like provenance",
            "notes": m["manual_evidence"],
        }
    )

    # Record the second zebrafish FMOD co-ortholog revealed by locus review.
    rows.append(
        {
            "Gene": "FMOD",
            "species": "zebrafish",
            "accession": "NP_001025243.1",
            "protein length": "342",
            "domain status": "not yet included in the canonical domain run",
            "SignalP status": "not yet included in the canonical SignalP run",
            "MSA/manual status": "supplementary duplicated co-ortholog; canonical-panel integration pending",
            "tree support": "not yet included in the final FMOD tree",
            "synteny support": (
                "manual tentative; fmoda lies beside prelp in the conserved ancestral FMOD-PRELP block"
            ),
            "final decision": "retain as supplementary teleost FMOD co-ortholog; complete protein QC before canonical-panel inclusion",
            "notes": (
                "Human FMOD BLASTP: 56.6% identity over 80% query coverage. "
                "The selected fmodb protein has broader coverage, while fmoda has stronger ancestral synteny."
            ),
        }
    )

    # Preserve the annotation-limited elephant-shark FMOD locus identified by
    # the completed all-target synteny review.
    rows.append(
        {
            "Gene": "FMOD",
            "species": "callorhinchus_milii",
            "accession": "XP_007897806.2",
            "protein length": "684 (exact GFF/FASTA reconstruction; compound model)",
            "domain status": "Pfam LRR-rich across both halves; two significant LRRNT caps (aa 54-76 and 383-412)",
            "SignalP status": "SignalP not run; hydrophobic N terminus plus a second signal-peptide-like segment at aa 330-347",
            "MSA/manual status": "manual/protein-BLAST QC supports a DCN/BGN-like N half and FMOD-best C half",
            "tree support": "not yet included in the final FMOD tree",
            "synteny support": (
                "manual tentative; SynVoy rescue overlaps LOC103182549 and the broad block retains "
                "five human anchors, four in order after reversal"
            ),
            "final decision": "retain locus as tentative FMOD evidence; exclude the compound protein model from the confident FMOD panel",
            "notes": (
                "The nine-CDS model is complete (2055-bp CDS, terminal stop, no internal stop), but its two SLRP-like "
                "halves indicate a probable fused/compound annotation. Full-protein BLASTP aligns residues 1-330 best "
                "to DCN (47.0% identity) and residues 376-684 best to FMOD (57.0% identity; 82.2% human-reference "
                "coverage). Transcript evidence or a revised gene model is required before ortholog inclusion."
            ),
        }
    )

    # Preserve the chicken BGN-locus/asporin-product conflict as a second provenance row.
    bgn_chicken = index(
        read_tsv(
            "analyses/protein_analysis/candidates/BGN_chicken_annotation_conflict_curation.tsv"
        ),
        "gene",
        "species",
    )[("BGN", "chicken")]
    rows.append(
        {
            "Gene": "BGN",
            "species": "chicken",
            "accession": "XP_414298.2",
            "protein length": "371",
            "domain status": "Good SLRP-like architecture (7 LRR-related hits); not gene-specific",
            "SignalP status": "SP; strong_signal_peptide (historical excluded-panel run)",
            "MSA/manual status": "historical MSA Good; product/reciprocal-hit conflict",
            "tree support": "focused reference tree: sister to human ASPN (99.3 SH-aLRT/100 UFBoot)",
            "synteny support": synteny_text("BGN", "chicken", synteny),
            "final decision": "exclude from canonical BGN protein panel",
            "notes": bgn_chicken["interpretation"],
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
