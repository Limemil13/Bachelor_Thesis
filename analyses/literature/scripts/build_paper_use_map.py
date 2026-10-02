#!/usr/bin/env python3
"""Build the literature-to-thesis map used during thesis writing."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LITERATURE_DIR = ROOT / "analyses" / "literature"
THESIS_DIR = ROOT / "thesis"


EXTRA_PAPERS = {
    "Gil2025": {
        "scope": "general/evolution",
        "evidence_type": "primary comparative genomics",
        "main_relevant_finding": (
            "Comparative genomic analyses support expansion of clustered vertebrate SLRP genes "
            "through whole-genome and local duplication, followed by lineage-specific retention, "
            "loss, rearrangement, and transcriptional divergence."
        ),
        "how_it_is_used": (
            "Supports the Introduction's family-evolution model and the interpretation that current "
            "clusters reflect duplication history rather than defining SLRP classes by chromosome."
        ),
        "caveat": (
            "It provides a family-level model; the present thesis does not independently reconstruct "
            "or date every duplication event."
        ),
        "source_location": "Abstract; Results and Discussion on clustered genes and duplication history (open-access article).",
    },
    "John2026": {
        "scope": "expression dataset",
        "evidence_type": "primary mouse genetics and transcriptomics",
        "main_relevant_finding": (
            "PPP1R15B-dependent control of ER-stress signalling and lipid metabolism is required for "
            "normal chondrocyte development, growth-plate organization, and skeletal growth."
        ),
        "how_it_is_used": (
            "Cited primarily as the source study for GSE305415/SRA samples of P0 Prx1-lineage "
            "chondroprogenitors used in the initial expression screen."
        ),
        "caveat": (
            "The selected samples are sorted neonatal knee-joint progenitors, not purified anatomical "
            "growth-plate zones; GEO sample wording about culture versus immediate extraction remains ambiguous."
        ),
        "source_location": "Abstract and Methods; GEO GSE305415 sample metadata.",
    },
    "Bi1999Sox9": {
        "scope": "growth-plate background",
        "evidence_type": "primary mouse developmental genetics",
        "main_relevant_finding": (
            "Sox9 is required for chondrogenic lineage formation and cartilage development and controls "
            "key cartilage-matrix genes."
        ),
        "how_it_is_used": "Supports the core non-SLRP network table describing SOX9 and cartilage ECM production.",
        "caveat": "Establishes a central developmental regulator, not a direct mechanism for the selected SLRPs.",
        "source_location": "Abstract and developmental phenotype/results sections.",
    },
    "Vortkamp1996IhhPthrp": {
        "scope": "growth-plate background",
        "evidence_type": "primary developmental signalling study",
        "main_relevant_finding": (
            "The IHH-PTHrP feedback circuit regulates the pace and spatial progression of chondrocyte hypertrophy."
        ),
        "how_it_is_used": "Supports the non-SLRP growth-plate network table and zonation background.",
        "caveat": "Defines a core signalling circuit but is not direct evidence for SLRP function.",
        "source_location": "Abstract and experimental results on cartilage differentiation.",
    },
    "Chen2014Runx2": {
        "scope": "growth-plate background",
        "evidence_type": "primary mouse developmental genetics",
        "main_relevant_finding": (
            "Runx2 regulates chondrocyte proliferation and maturation during endochondral ossification."
        ),
        "how_it_is_used": "Supports the RUNX2/COL10A1 entry in the core growth-plate gene-network table.",
        "caveat": "RUNX2 regulation provides developmental context rather than evidence about SLRP orthology.",
        "source_location": "Abstract; Results and Discussion (open-access full text, PMCID PMC4535340).",
    },
    "Deng1996Fgfr3": {
        "scope": "growth-plate background",
        "evidence_type": "primary mouse genetics",
        "main_relevant_finding": (
            "Fgfr3 disruption increases endochondral bone growth and expands proliferating and hypertrophic "
            "growth-plate compartments, establishing FGFR3 as a negative regulator of bone growth."
        ),
        "how_it_is_used": "Supports the FGFR3 entry in the core growth-plate gene-network table.",
        "caveat": "The cited mouse loss-of-function study does not by itself cover every human activating phenotype.",
        "source_location": "Abstract and skeletal phenotype/results sections.",
    },
    "Tamura2004Npr2": {
        "scope": "growth-plate background",
        "evidence_type": "primary mouse genetics",
        "main_relevant_finding": (
            "Loss of the guanylyl cyclase-B receptor NPR2 impairs endochondral ossification and longitudinal growth."
        ),
        "how_it_is_used": "Supports the NPPC/NPR2 entry in the core growth-plate gene-network table.",
        "caveat": "The paper directly tests NPR2; the table's ligand-receptor framing also relies on established pathway knowledge.",
        "source_location": "Abstract; Results and Discussion (open-access full text, PMCID PMC534612).",
    },
    "Inada2004Mmp13": {
        "scope": "growth-plate background",
        "evidence_type": "primary mouse genetics",
        "main_relevant_finding": (
            "Mmp13 loss causes abnormal hypertrophic cartilage remodelling, collagen accumulation, and delayed "
            "ossification and vascular invasion."
        ),
        "how_it_is_used": "Supports the MMP13 entry in the core growth-plate gene-network table.",
        "caveat": "Describes a matrix-remodelling enzyme and does not show that SLRPs are its only relevant substrates.",
        "source_location": "Abstract; Results and Discussion (open-access full text, PMCID PMC535367).",
    },
    "Gerber1999Vegf": {
        "scope": "growth-plate background",
        "evidence_type": "primary developmental/angiogenesis study",
        "main_relevant_finding": (
            "VEGF from hypertrophic cartilage couples vascular invasion to cartilage remodelling and trabecular bone formation."
        ),
        "how_it_is_used": "Supports the VEGFA entry in the core growth-plate gene-network table.",
        "caveat": "Provides transition-zone context rather than direct evidence for any selected SLRP.",
        "source_location": "Abstract and experimental results on vascular invasion and ossification.",
    },
    "Bailey2015MEME": {
        "scope": "method",
        "evidence_type": "software/method reference",
        "main_relevant_finding": "Describes the MEME Suite for motif discovery and motif-based sequence analysis.",
        "how_it_is_used": "Cites the software suite used for protein motif discovery and MAST scanning.",
        "caveat": "A discovered statistical motif is not automatically a biochemical functional motif.",
        "source_location": "Abstract and MEME Suite tool descriptions.",
    },
    "Bailey2021STREME": {
        "scope": "method",
        "evidence_type": "algorithm reference",
        "main_relevant_finding": "Introduces STREME for discriminative discovery of enriched sequence motifs.",
        "how_it_is_used": "Cites the method used to compare each gene's protein sequences against the other genes.",
        "caveat": "No holdout set was used in this exploratory thesis analysis, so motif enrichment is descriptive.",
        "source_location": "Abstract, algorithm description, and benchmarking sections.",
    },
    "Ovek2026JASPAR": {
        "scope": "method/database",
        "evidence_type": "database release reference",
        "main_relevant_finding": (
            "Documents the 2026 JASPAR release and its expanded curated transcription-factor binding-profile collections."
        ),
        "how_it_is_used": "Cites the JASPAR 2026 vertebrate CORE profiles used for promoter motif comparison.",
        "caveat": "Motif matches predict binding compatibility, not transcription-factor occupancy or regulation in growth plates.",
        "source_location": "Abstract and JASPAR CORE/database-update sections.",
    },
    "Baldarelli2024MGI": {
        "scope": "method/database",
        "evidence_type": "database reference",
        "main_relevant_finding": "Describes Mouse Genome Informatics as an integrated mouse genetics and phenotype knowledgebase.",
        "how_it_is_used": "Cites the MGI genotype-phenotype data used for the mouse phenotype summary.",
        "caveat": "Annotation counts reflect available curated evidence and study intensity, not biological effect size.",
        "source_location": "Abstract and MGI data/resource description.",
    },
    "Gargano2024HPO": {
        "scope": "method/database",
        "evidence_type": "ontology/database reference",
        "main_relevant_finding": "Describes the Human Phenotype Ontology resource and its gene-disease-phenotype annotations.",
        "how_it_is_used": "Cites the HPO release used for human phenotype and disease annotation.",
        "caveat": "Absence of an HPO record is absence of curated evidence, not proof of no human phenotype.",
        "source_location": "Abstract and resource/data-content sections.",
    },
    "NeiGojobori1986": {
        "scope": "method",
        "evidence_type": "evolutionary-rate method",
        "main_relevant_finding": (
            "Introduces a counting-based method for estimating synonymous and nonsynonymous substitutions."
        ),
        "how_it_is_used": "Cites the pairwise dN/dS calculation used for human-reference comparisons.",
        "caveat": "Short or divergent pairwise sequences and low or undefined dS make individual ratios unstable.",
        "source_location": "Method derivation and simulation/accuracy sections.",
    },
    "Bianco1990DecorinBiglycan": {
        "scope": "BGN/DCN",
        "evidence_type": "primary developmental localization",
        "main_relevant_finding": (
            "Biglycan and decorin show distinct but overlapping localization patterns in developing human skeletal "
            "and non-skeletal connective tissues."
        ),
        "how_it_is_used": "Not currently cited; available for a concise DCN/BGN developmental-localization statement.",
        "caveat": "Localization supports tissue relevance but does not establish a gene-specific molecular mechanism.",
        "source_location": "Abstract and tissue-localization results.",
    },
    "Chery2021DecorinCartilage": {
        "scope": "DCN",
        "evidence_type": "primary cartilage mechanobiology",
        "main_relevant_finding": (
            "Decorin contributes to cartilage pericellular-matrix organization and micromechanical properties."
        ),
        "how_it_is_used": "Not currently cited; suitable support for the DCN cartilage-matrix profile.",
        "caveat": "Cartilage micromechanics evidence does not by itself prove a growth-plate-zone-specific role.",
        "source_location": "Abstract; Results and Discussion (open-access full text, PMCID PMC7902451).",
    },
    "Danielson1997DecorinKnockout": {
        "scope": "DCN",
        "evidence_type": "primary mouse genetics",
        "main_relevant_finding": (
            "Decorin loss produces abnormal collagen fibrils and tissue fragility, demonstrating an in-vivo role in fibrillogenesis."
        ),
        "how_it_is_used": "Not currently cited; suitable support for a conserved DCN collagen-organization function.",
        "caveat": "The strongest phenotype is in skin; transfer to growth-plate cartilage requires caution.",
        "source_location": "Abstract; Results and Discussion (open-access full text, PMCID PMC2134287).",
    },
    "LumicanReview2024": {
        "scope": "LUM/background",
        "evidence_type": "review",
        "main_relevant_finding": (
            "Reviews lumican as a multifunctional extracellular proteoglycan with roles in collagen fibrillogenesis, "
            "cell behaviour, inflammation, and disease."
        ),
        "how_it_is_used": "Background PDF only; not currently present in the bibliography or cited thesis text.",
        "caveat": "Broad review evidence should not replace primary cartilage or bone studies for thesis-specific claims.",
        "source_location": "Local PDF pp. 1, 7, and 14 (overview and collagen-fibrillogenesis sections).",
        "title": "Lumican, a multifunctional cell instructive biomarker proteoglycan",
        "year": "2024",
        "doi": "10.3390/ijms25052825",
        "source_access": "Local PDF: modern_biology_behindLUM.pdf",
    },
}


SOURCE_LOCATIONS = {
    "IozzoSchaefer2015": "Local PDF pp. 61-68 (SLRP nomenclature, classes, and LRR-rich architecture).",
    "SchaeferIozzo2008": "Local PDF pp. 1 and 4-5 (overview and signalling/growth-factor functions).",
    "Kram2017": "Local PDF pp. 1-2 and 8-13 (abstract, functional assays, and discussion).",
    "Saamanen2001": "Local PDF pp. 1 and 4-8 (developmental expression/localization results).",
    "Johnson1999": "Local PDF pp. 1 and 5-10 (growth-plate expression and protein-localization results).",
    "Grover1996": "Local PDF pp. 1 and 3-7 (gene organization and cartilage-expression results).",
    "Grover2002": "Local PDF pp. 1 and 4-9 (three-exon gene structure and cartilage expression).",
    "Rucci2009": "Local PDF pp. 1-3 and 5-14 (osteoclastogenesis mechanism and discussion).",
    "Li2016PRELP": "Local PDF pp. 1-4 (osteoblastic differentiation and mineralization experiments).",
    "Costa2018": "Local PDF pp. 1-3 and 6-15 (family evolution, teleost duplication, and expression analyses).",
    "Ge2004": "Local PDF pp. 1-2 and 4-7 (processing and collagen-fibrillogenesis assays).",
    "Raouf2002": "Local PDF pp. 1-4 and 6 (osteoblast differentiation and bone-matrix localization).",
    "Kafienah2008": "Local PDF pp. 1 and 4-8 (type-II collagen accumulation and fibril measurements).",
    "Nikitovic2012": "Local PDF pp. 1-8 (bone-focused SLRP review).",
    "Wilda2000": "Local PDF pp. 1 and 4-10 (comparative developmental expression patterns).",
    "Kram2020": "Local PDF pp. 1-4 and 14-22 (OPG-Fc experiment and age-dependent outcomes).",
    "Rucci2013": "Local PDF pp. 1-2 and 4-12 (PRELP-derived peptide and bone-loss models).",
    "Shu2019": "Local PDF pp. 1-4 and 10-13 (growth-plate localization and proteolysis).",
    "Lee2021Lum": "Local PDF pp. 1-6 and 8-10 (Akt, osteoclastogenesis, and resorption assays).",
    "Park2008": "Abstract, Results, and Discussion (open-access full text, PMCID PMC2637281).",
    "Sommarin1998": "Abstract and Results/Discussion (open-access full text, PMCID PMC2132750).",
    "Lin2021OMD": "Abstract and Results/Discussion (open-access full text, PMCID PMC7862363).",
    "Corsi2002": "Abstract and phenotype results (PubMed record PMID 12102052).",
    "Svensson1999": "Abstract and collagen-fibril results (PubMed record PMID 10092650).",
    "Nuka2010": "Abstract and Results/Discussion (open-access full text, PMCID PMC3013283).",
    "Bengtsson2000": "Abstract and biochemical binding results (PubMed record PMID 11007795).",
    "Tasheva2002": "Abstract and collagen-fibril phenotype results (PubMed record PMID 12432342).",
}


METHOD_KEYS = {
    "Bailey2015MEME",
    "Bailey2021STREME",
    "Ovek2026JASPAR",
    "Baldarelli2024MGI",
    "Gargano2024HPO",
    "NeiGojobori1986",
}


def parse_bibliography(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, dict[str, str]] = {}
    for match in re.finditer(r"(?ms)^@\w+\{([^,]+),\s*(.*?)^\}", text):
        key, body = match.groups()
        fields: dict[str, str] = {}
        for line in body.splitlines():
            field_match = re.match(r"\s*(\w+)\s*=\s*\{(.*)\},?\s*$", line)
            if field_match:
                fields[field_match.group(1).lower()] = field_match.group(2)
        entries[key] = fields
    return entries


def clean_tex(text: str) -> str:
    text = re.sub(r"\\cite\w*\{[^}]+\}", "", text)
    for command in ("textit", "textbf", "texttt", "emph"):
        text = re.sub(rf"\\{command}\{{([^{{}}]*)\}}", r"\1", text)
    text = text.replace("\\&", "&").replace("\\%", "%")
    text = text.replace("\\", " ").replace("~", " ")
    text = re.sub(r"\\[A-Za-z]+\*?(?:\[[^\]]*\])?", "", text)
    text = text.replace("{", "").replace("}", "")
    return " ".join(text.split()).strip(" &")


def citation_uses(
    chapter_dir: Path,
) -> tuple[set[str], dict[str, list[str]], dict[str, list[str]]]:
    cited: set[str] = set()
    locations: dict[str, list[str]] = defaultdict(list)
    statements: dict[str, list[str]] = defaultdict(list)
    heading_pattern = re.compile(r"\\(section|subsection|subsubsection)\{([^}]*)\}")
    citation_pattern = re.compile(r"\\cite\w*\{([^}]+)\}")

    for path in sorted(chapter_dir.glob("*.tex")):
        lines = path.read_text(encoding="utf-8").splitlines()
        heading = path.stem.capitalize()
        for index, line in enumerate(lines):
            heading_match = heading_pattern.search(line)
            if heading_match:
                heading = clean_tex(heading_match.group(2))
            citation_matches = list(citation_pattern.finditer(line))
            if not citation_matches:
                continue
            start = max(0, index - 3)
            end = min(len(lines), index + 2)
            context = clean_tex(" ".join(lines[start:end]))
            if len(context) > 520:
                context = context[:517].rstrip() + "..."
            for citation_match in citation_matches:
                for key in (
                    item.strip() for item in citation_match.group(1).split(",")
                ):
                    cited.add(key)
                    location = f"{path.name} > {heading} (line {index + 1})"
                    if location not in locations[key]:
                        locations[key].append(location)
                    if context and context not in statements[key]:
                        statements[key].append(context)
    return cited, locations, statements


def load_evidence(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return {
            row["bibtex_key"]: row for row in csv.DictReader(handle, delimiter="\t")
        }


def source_identifier(fields: dict[str, str], evidence: dict[str, str] | None) -> str:
    if evidence and evidence.get("doi_or_pmid"):
        return evidence["doi_or_pmid"]
    identifiers = []
    if fields.get("doi"):
        identifiers.append(f"DOI:{fields['doi']}")
    if fields.get("pmid"):
        identifiers.append(f"PMID:{fields['pmid']}")
    if fields.get("pmcid"):
        identifiers.append(f"PMCID:{fields['pmcid']}")
    return "; ".join(identifiers)


def source_access(fields: dict[str, str], evidence: dict[str, str] | None) -> str:
    if evidence and evidence.get("private_pdf_or_open_url"):
        value = evidence["private_pdf_or_open_url"]
        return f"Local PDF: {value}" if not value.startswith("http") else value
    if fields.get("url"):
        return fields["url"]
    if fields.get("pmcid"):
        return f"https://pmc.ncbi.nlm.nih.gov/articles/{fields['pmcid']}/"
    if fields.get("pmid"):
        return f"https://pubmed.ncbi.nlm.nih.gov/{fields['pmid']}/"
    if fields.get("doi"):
        return f"https://doi.org/{fields['doi']}"
    return "Bibliography metadata only"


def markdown_escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def main() -> None:
    bibliography = parse_bibliography(THESIS_DIR / "bibliography.bib")
    evidence = load_evidence(LITERATURE_DIR / "slrp_literature_evidence.tsv")
    cited, locations, statements = citation_uses(THESIS_DIR / "chapters")

    keys = list(bibliography)
    keys.append("LumicanReview2024")
    rows = []
    for key in keys:
        fields = bibliography.get(key, {})
        evidence_row = evidence.get(key)
        extra = EXTRA_PAPERS.get(key, {})

        scope = (evidence_row or {}).get("gene") or extra.get("scope", "general")
        evidence_type = (evidence_row or {}).get("evidence_type") or extra.get(
            "evidence_type", "paper"
        )
        finding = (evidence_row or {}).get("paraphrased_finding") or extra.get(
            "main_relevant_finding", ""
        )
        thesis_use = (evidence_row or {}).get(
            "link_to_computational_result"
        ) or extra.get("how_it_is_used", "")
        caveat = (evidence_row or {}).get("caveat") or extra.get("caveat", "")
        title = extra.get("title") or clean_tex(fields.get("title", key))
        year = extra.get("year") or fields.get("year", "")
        identifier = source_identifier(fields, evidence_row)
        access = extra.get("source_access") or source_access(fields, evidence_row)
        source_location = SOURCE_LOCATIONS.get(key) or extra.get(
            "source_location",
            "Abstract and relevant Results/Discussion sections; exact page not recorded.",
        )

        if key in cited:
            status = "cited"
            current_location = "; ".join(locations[key])
            current_statement = " || ".join(statements[key])
        else:
            status = "not currently cited"
            current_location = "Not currently cited in thesis chapter TeX"
            current_statement = ""

        category = (
            "method_or_data_source"
            if key in METHOD_KEYS or key == "John2026"
            else "biology_and_interpretation"
        )
        if status != "cited":
            category = "relevant_not_currently_cited"

        rows.append(
            {
                "bibtex_key": key,
                "year": year,
                "title": title,
                "scope": scope,
                "category": category,
                "evidence_type": evidence_type,
                "main_relevant_finding": finding,
                "how_it_is_used": thesis_use,
                "current_thesis_location": current_location,
                "current_thesis_statement": current_statement,
                "source_identifier": identifier or extra.get("doi", ""),
                "source_access": access,
                "source_location_checked": source_location,
                "caveat": caveat,
                "citation_status": status,
            }
        )

    category_order = {
        "biology_and_interpretation": 0,
        "method_or_data_source": 1,
        "relevant_not_currently_cited": 2,
    }
    rows.sort(
        key=lambda row: (
            category_order[row["category"]],
            row["current_thesis_location"],
            row["bibtex_key"],
        )
    )

    tsv_path = LITERATURE_DIR / "paper_evidence_and_thesis_use.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    md_path = LITERATURE_DIR / "PAPER_EVIDENCE_AND_THESIS_USE.md"
    lines = [
        "# Paper evidence and use in the thesis",
        "",
        "This is the writing-oriented literature map. It distinguishes what each paper supports, where it is currently cited, and what the evidence does not prove. The line numbers describe the current LaTeX draft and can change during rewriting; rerun `scripts/build_paper_use_map.py` to refresh them.",
        "",
        f"Coverage: {len(rows)} unique papers ({sum(row['citation_status'] == 'cited' for row in rows)} currently cited; {sum(row['citation_status'] != 'cited' for row in rows)} relevant but not currently cited). The second local copy of the Costa et al. teleost article is recorded as a duplicate, not as a separate paper.",
        "",
    ]

    section_titles = {
        "biology_and_interpretation": "Biological and evolutionary evidence currently cited",
        "method_or_data_source": "Method and data-source references currently cited",
        "relevant_not_currently_cited": "Relevant papers not currently cited",
    }
    for category in section_titles:
        lines.extend(
            [
                f"## {section_titles[category]}",
                "",
                "| Paper | Main thesis-relevant message | Current thesis use and location | Source check | Main caution |",
                "|---|---|---|---|---|",
            ]
        )
        for row in (item for item in rows if item["category"] == category):
            paper = f"`{row['bibtex_key']}` ({row['year']}), {row['title']}"
            thesis_cell = f"{row['how_it_is_used']} **Location:** {row['current_thesis_location']}"
            source_cell = f"{row['source_location_checked']} {row['source_access']}"
            lines.append(
                "| "
                + " | ".join(
                    markdown_escape(value)
                    for value in (
                        paper,
                        row["main_relevant_finding"],
                        thesis_cell,
                        source_cell,
                        row["caveat"],
                    )
                )
                + " |"
            )
        lines.append("")

    lines.extend(
        [
            "## How to use this during rewriting",
            "",
            "For each paragraph, first identify the biological claim. Then use `current_thesis_location` to find the present wording, check `main_relevant_finding`, and retain the caveat. A citation supports only the claim described here; it should not be stretched into proof of conserved function across all species.",
            "",
            "The detailed TSV contains the current citation context as well as the compact summary shown above.",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote {tsv_path.relative_to(ROOT)} ({len(rows)} rows)")
    print(f"Wrote {md_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
