#!/usr/bin/env python3
"""Remove DCN contamination from the compact BGN protein panel.

Four historical BGN records are DCN proteins. Annotated biglycan-like models
replace the catshark and whale-shark loci; the workflow removes
the opossum/elephant-shark records because their matched NCBI annotations do not
provide a separate supported BGN model. The two replacements remain tentative
partial/deep-lineage candidates.
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq
from repair_bgn_whale_shark import read_fasta, run_blast, write_fasta

ROOT = Path(__file__).resolve().parents[3]
if not (synvoy_root := os.environ.get("SYNVOY_ROOT")):
    raise RuntimeError("Set SYNVOY_ROOT to the SynVoy checkout")
SYNVOY = Path(synvoy_root)
REFERENCE = (
    ROOT
    / "analyses/phylogenetics/combined_trees/slrp_family_reference_tree"
    / "SLRP_reference_context.faa"
)
TARGET_FILES = (
    ROOT / "analyses/protein_analysis/candidates/BGN_domain_input.faa",
    ROOT
    / "analyses/protein_analysis/candidates/BGN_domain_input_chicken_asporin_excluded.faa",
)
OUTPUT_FILES = (
    ROOT
    / "analyses/protein_analysis/candidates/BGN_domain_input_decontaminated_review.faa",
    ROOT
    / "analyses/protein_analysis/candidates/BGN_domain_input_canonical_decontaminated.faa",
)
AUDIT = ROOT / "analyses/protein_analysis/candidates/BGN_deep_lineage_correction.tsv"


@dataclass(frozen=True)
class Replacement:
    species: str
    gff_name: str
    fna_name: str
    accession: str
    removed_accession: str
    status: str


REPLACEMENTS = (
    Replacement(
        "scyliorhinus_canicula",
        "catshark.gff",
        "catshark.fna",
        "XP_038638277.1",
        "XP_038636760.1",
        "tentative_biglycan_like_model",
    ),
    Replacement(
        "whale_shark",
        "whale_shark.gff",
        "whale_shark.fna",
        "XP_048475905.1",
        "XP_048463588.1",
        "tentative_partial_biglycan_like_model",
    ),
)
EXCLUSIONS = {
    "monodelphis_domestica": (
        "XP_001363160.2",
        "excluded_DCN_contamination_no_separate_supported_BGN_model_in_target_GFF",
    ),
    "callorhinchus_milii": (
        "NP_001279248.1",
        "excluded_DCN_contamination_no_separate_supported_BGN_model_in_target_GFF",
    ),
}


def parse_attrs(text: str) -> dict[str, str]:
    return dict(
        item.split("=", 1) for item in text.rstrip(";").split(";") if "=" in item
    )


def reconstruct(replacement: Replacement) -> tuple[str, str, str, int, int]:
    gff = SYNVOY / "pro_panel/targets_15/gff" / replacement.gff_name
    genome = SYNVOY / "pro_panel/targets_15/fna" / replacement.fna_name
    cds: list[tuple[str, int, int, str]] = []
    with gff.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") or replacement.accession not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("protein_id") == replacement.accession:
                cds.append((fields[0], int(fields[3]), int(fields[4]), fields[6]))
    if not cds:
        raise ValueError(f"No CDS for {replacement.accession} in {gff}")
    strands, seqids = {row[3] for row in cds}, {row[0] for row in cds}
    if len(strands) != 1 or len(seqids) != 1:
        raise ValueError(f"Inconsistent locus for {replacement.accession}")
    strand = next(iter(strands))
    ordered = sorted(cds, key=lambda row: row[1], reverse=strand == "-")
    pieces: list[str] = []
    for seqid, start, end, _strand in ordered:
        output = subprocess.check_output(
            ["samtools", "faidx", str(genome), f"{seqid}:{start}-{end}"], text=True
        )
        sequence = "".join(output.splitlines()[1:]).upper()
        pieces.append(
            str(Seq(sequence).reverse_complement()) if strand == "-" else sequence
        )
    coding = "".join(pieces)
    if len(coding) % 3:
        raise ValueError(f"Non-triplet CDS for {replacement.accession}: {len(coding)}")
    protein = str(Seq(coding).translate())
    protein = protein.removesuffix("*")
    if "*" in protein:
        raise ValueError(f"Internal stop in {replacement.accession}")
    return protein, ordered[0][0], strand, len(ordered), len(coding)


def validate(
    protein: str, human: str, family_records: list[tuple[str, str]], accession: str
) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="bgn_deep_repair_") as temp_name:
        temp = Path(temp_name)
        human_fasta, target_fasta, family_fasta = (
            temp / "human.faa",
            temp / "target.faa",
            temp / "family.faa",
        )
        write_fasta(human_fasta, [("human_BGN", human)])
        write_fasta(target_fasta, [(accession, protein)])
        write_fasta(family_fasta, family_records)
        forward = run_blast(
            human_fasta, target_fasta, "pident qcovs length evalue bitscore"
        )
        reverse = run_blast(
            target_fasta, family_fasta, "sseqid pident qcovs evalue bitscore"
        )
    reverse_gene, reverse_accession = reverse[0].split("__", 1)
    if reverse_gene != "BGN":
        raise ValueError(f"{accession} reverse best gene is {reverse_gene}, not BGN")
    return {
        "forward_identity_percent": f"{float(forward[0]):.2f}",
        "forward_query_coverage_percent": f"{float(forward[1]):.2f}",
        "forward_alignment_length": forward[2],
        "forward_evalue": forward[3],
        "forward_bitscore": forward[4],
        "reverse_best_gene": reverse_gene,
        "reverse_best_accession": reverse_accession,
        "reverse_identity_percent": f"{float(reverse[1]):.2f}",
        "reverse_query_coverage_percent": f"{float(reverse[2]):.2f}",
        "reverse_evalue": reverse[3],
        "reverse_bitscore": reverse[4],
    }


def main() -> None:
    current = read_fasta(TARGET_FILES[0])
    human = next(
        sequence for header, sequence in current if header.startswith("human|BGN|")
    )
    family_records: list[tuple[str, str]] = []
    for header, sequence in read_fasta(REFERENCE):
        gene, _class_name, accession = header.split()[0].split("|")
        family_records.append((f"{gene}__{accession}", sequence))

    replacement_data: dict[str, tuple[str, str, dict[str, str]]] = {}
    audit_rows: list[dict[str, object]] = []
    for replacement in REPLACEMENTS:
        protein, seqid, strand, block_count, cds_length = reconstruct(replacement)
        metrics = validate(protein, human, family_records, replacement.accession)
        header = (
            f"{replacement.species}|BGN|{replacement.accession} "
            f"forward_pident={metrics['forward_identity_percent']} "
            f"forward_qcovs={metrics['forward_query_coverage_percent']} "
            f"reverse_hit={metrics['reverse_best_accession']}"
        )
        replacement_data[replacement.species] = (header, protein, metrics)
        audit_rows.append(
            {
                "gene": "BGN",
                "species": replacement.species,
                "action": "replace",
                "removed_accession": replacement.removed_accession,
                "replacement_accession": replacement.accession,
                "final_status": replacement.status,
                "seqid": seqid,
                "strand": strand,
                "cds_block_count": block_count,
                "cds_length_bp": cds_length,
                "protein_length_aa": len(protein),
                **metrics,
                "evidence": "replacement reconstructed from the exact NCBI target GFF3/genome; reciprocal best human SLRP gene is BGN",
            }
        )
    for species, (accession, status) in EXCLUSIONS.items():
        audit_rows.append(
            {
                "gene": "BGN",
                "species": species,
                "action": "exclude",
                "removed_accession": accession,
                "replacement_accession": "",
                "final_status": status,
                "seqid": "",
                "strand": "",
                "cds_block_count": "",
                "cds_length_bp": "",
                "protein_length_aa": "",
                "forward_identity_percent": "",
                "forward_query_coverage_percent": "",
                "forward_alignment_length": "",
                "forward_evalue": "",
                "forward_bitscore": "",
                "reverse_best_gene": "DCN",
                "reverse_best_accession": "",
                "reverse_identity_percent": "",
                "reverse_query_coverage_percent": "",
                "reverse_evalue": "",
                "reverse_bitscore": "",
                "evidence": "historical protein is a DCN candidate and no separate supported BGN model was found in the matched target GFF",
            }
        )

    for path, output_path in zip(TARGET_FILES, OUTPUT_FILES, strict=False):
        output: list[tuple[str, str]] = []
        seen_replacements: set[str] = set()
        for header, sequence in read_fasta(path):
            species, gene, _accession = header.split()[0].split("|")
            if gene == "BGN" and species in EXCLUSIONS:
                continue
            if gene == "BGN" and species in replacement_data:
                if species not in seen_replacements:
                    new_header, new_sequence, _metrics = replacement_data[species]
                    output.append((new_header, new_sequence))
                    seen_replacements.add(species)
                continue
            output.append((header, sequence))
        if seen_replacements != set(replacement_data):
            raise ValueError(f"Did not replace every deep-lineage BGN record in {path}")
        expected_count = 13 if "chicken_asporin_excluded" in path.name else 14
        if len(output) != expected_count:
            raise ValueError(
                f"Review BGN panel from {path.name} must contain {expected_count} "
                f"proteins, found {len(output)}"
            )
        write_fasta(output_path, output)

    fields = list(audit_rows[0])
    with AUDIT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)
    print("Proposed decontaminated BGN review panel: 13 proteins")
    for row in audit_rows:
        print(
            f"{row['species']}: {row['action']} -> {row['replacement_accession'] or 'no confident BGN protein'}"
        )
    print(f"Wrote {AUDIT}")
    for output_path in OUTPUT_FILES:
        print(f"Wrote review-only panel {output_path}")


if __name__ == "__main__":
    main()
