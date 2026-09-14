#!/usr/bin/env python3
"""Fetch current NCBI records for the six predicted transcript representatives.

This produces a dated annotation audit. It does not replace sequence/synteny
judgment where NCBI's gene symbol and translated product are discordant.
"""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "analyses"
    / "overview"
    / "gene_structure_predicted_transcript_ncbi_review.tsv"
)

TARGETS = [
    ("BGN", "chicken", "XM_414298.8"),
    ("FMOD", "zebrafish", "XM_005167241.6"),
    ("OMD", "chicken", "XM_004944504.5"),
    ("OMD", "zebrafish", "XM_017353401.4"),
    ("PRELP", "chicken", "XM_418054.8"),
    ("PRELP", "zebrafish", "XM_001923555.9"),
]

CONCLUSIONS = {
    "XM_414298.8": (
        "unresolved annotation conflict",
        "Current NCBI MODEL record assigns the locus/gene symbol BGN but names "
        "the translated product asporin; the compact protein also reverse-hits "
        "human ASPN. Retain only as a tentative chicken BGN-locus representative "
        "pending neighborhood/synteny confirmation.",
    ),
    "XM_005167241.6": (
        "database review complete; retain with teleost-paralog wording",
        "Current NCBI MODEL transcript is fmodb (fibromodulin b), variant X1. "
        "Use the explicit fmodb/teleost FMOD-family label rather than implying "
        "that zebrafish has a single unsuffixed FMOD gene.",
    ),
    "XM_004944504.5": (
        "database review complete; retain predicted model",
        "Current NCBI MODEL transcript is chicken OMD/osteomodulin; NCBI also "
        "lists an alternative predicted splice record, so this is a selected "
        "representative model rather than a curated NM transcript.",
    ),
    "XM_017353401.4": (
        "database review complete; retain predicted model",
        "Current NCBI MODEL transcript is zebrafish omd/osteomodulin and matches "
        "the expected gene assignment.",
    ),
    "XM_418054.8": (
        "database review complete; retain predicted model",
        "Current NCBI MODEL transcript encodes chicken PRELP/prolargin isoform X2; "
        "it is a predicted representative, not a curated NM transcript.",
    ),
    "XM_001923555.9": (
        "database review complete; retain predicted model",
        "Current NCBI MODEL transcript encodes zebrafish prelp/prolargin and "
        "matches the expected gene assignment.",
    ),
}


def fetch_records(accessions: list[str]):
    params = urlencode(
        {
            "db": "nuccore",
            "id": ",".join(accessions),
            "rettype": "gb",
            "retmode": "text",
            "tool": "bachelor_thesis_local_audit",
            "email": "local-audit@example.invalid",
        }
    )
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    with urlopen(url, timeout=120) as response:
        text = response.read().decode("utf-8")
    return list(SeqIO.parse(StringIO(text), "genbank"))


def first_qualifier(feature, key: str) -> str:
    values = feature.qualifiers.get(key, [])
    return values[0] if values else ""


def main() -> None:
    records = fetch_records([accession for _, _, accession in TARGETS])
    by_accession = {record.id: record for record in records}
    expected = {accession for _, _, accession in TARGETS}
    if set(by_accession) != expected:
        raise SystemExit(
            f"NCBI record mismatch: missing={sorted(expected - set(by_accession))}"
        )

    rows = []
    for expected_gene, species, accession in TARGETS:
        record = by_accession[accession]
        gene_feature = next(
            feature for feature in record.features if feature.type == "gene"
        )
        cds = next(feature for feature in record.features if feature.type == "CDS")
        ncbi_gene = first_qualifier(gene_feature, "gene")
        product = first_qualifier(cds, "product")
        protein_id = first_qualifier(cds, "protein_id")
        translation = first_qualifier(cds, "translation")
        status, conclusion = CONCLUSIONS[accession]
        rows.append(
            {
                "gene": expected_gene,
                "species": species,
                "transcript_accession": accession,
                "record_date": record.annotations.get("date", ""),
                "definition": record.description,
                "refseq_model_status": "MODEL / predicted (XM_)",
                "ncbi_gene_symbol": ncbi_gene,
                "ncbi_product": product,
                "protein_accession": protein_id,
                "translated_length_aa": len(translation),
                "expected_gene_symbol_compatible": (
                    "yes"
                    if ncbi_gene.upper() == expected_gene
                    or (expected_gene == "FMOD" and ncbi_gene.lower() == "fmodb")
                    else "no"
                ),
                "review_status": status,
                "review_conclusion": conclusion,
                "ncbi_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
            }
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} current NCBI transcript reviews to {OUTPUT}")


if __name__ == "__main__":
    main()
