#!/usr/bin/env python3
"""Replace the historical whale-shark DCN sequence mislabeled as BGN.

NCBI GFF3 assigns XP_048463588.1 to dcn (GeneID 109920536). The actual
biglycan-like locus on the thesis whale-shark assembly encodes the partial
protein XP_048475905.1 (GeneID 125487751). The replacement is reconstructed
from the matching GFF3/genome pair and validated against the human SLRP
reference panel with BLASTP.
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[3]
if not (synvoy_root := os.environ.get("SYNVOY_ROOT")):
    raise RuntimeError("Set SYNVOY_ROOT to the SynVoy checkout")
SYNVOY = Path(synvoy_root)
GFF = SYNVOY / "pro_panel/targets_15/gff/whale_shark.gff"
GENOME = SYNVOY / "pro_panel/targets_15/fna/whale_shark.fna"
TARGET_ACCESSION = "XP_048475905.1"
WRONG_ACCESSION = "XP_048463588.1"
TARGET_FILES = (
    ROOT / "analyses/protein_analysis/candidates/BGN_domain_input.faa",
    ROOT
    / "analyses/protein_analysis/candidates/BGN_domain_input_chicken_asporin_excluded.faa",
)
REFERENCE = (
    ROOT
    / "analyses/phylogenetics/combined_trees/slrp_family_reference_tree"
    / "SLRP_reference_context.faa"
)
AUDIT = ROOT / "analyses/protein_analysis/candidates/BGN_whale_shark_correction.tsv"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    # Preserve source headers so the exact replaced accession can be audited.
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    # Write the repaired panel to a new FASTA; never overwrite the source panel.
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def parse_attrs(text: str) -> dict[str, str]:
    # Parse the GFF3 attribute column for CDS parent and protein identifiers.
    return dict(
        item.split("=", 1) for item in text.rstrip(";").split(";") if "=" in item
    )


def reconstruct_protein() -> tuple[str, str, str, int, int]:
    # Find CDS records for the target accession, extract them from the matching
    # genome assembly in transcript order and translate the reconstructed CDS.
    cds: list[tuple[str, int, int, str, int]] = []
    with GFF.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") or TARGET_ACCESSION not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("protein_id") != TARGET_ACCESSION:
                continue
            cds.append(
                (fields[0], int(fields[3]), int(fields[4]), fields[6], int(fields[7]))
            )
    if not cds:
        raise ValueError(f"No {TARGET_ACCESSION} CDS records in {GFF}")
    strands = {item[3] for item in cds}
    seqids = {item[0] for item in cds}
    if len(strands) != 1 or len(seqids) != 1:
        raise ValueError("Target CDS is not on one sequence/strand")
    strand = strands.pop()
    ordered = sorted(cds, key=lambda item: item[1], reverse=strand == "-")
    pieces: list[str] = []
    for seqid, start, end, _strand, _phase in ordered:
        result = subprocess.run(
            ["samtools", "faidx", str(GENOME), f"{seqid}:{start}-{end}"],
            check=True,
            capture_output=True,
            text=True,
        )
        genomic = "".join(result.stdout.splitlines()[1:]).upper()
        pieces.append(
            str(Seq(genomic).reverse_complement()) if strand == "-" else genomic
        )
    coding = "".join(pieces)
    if len(coding) % 3:
        raise ValueError(
            f"Reconstructed CDS length is not divisible by three: {len(coding)}"
        )
    protein = str(Seq(coding).translate())
    protein = protein.removesuffix("*")
    if "*" in protein:
        raise ValueError("Internal stop in reconstructed whale-shark BGN-like model")
    return protein, ordered[0][0], strand, len(cds), len(coding)


def run_blast(query: Path, subject: Path, fields: str) -> list[str]:
    # Validate the replacement against a fixed reference FASTA with BLASTP.
    result = subprocess.run(
        [
            "blastp",
            "-query",
            str(query),
            "-subject",
            str(subject),
            "-max_target_seqs",
            "1",
            "-max_hsps",
            "1",
            "-outfmt",
            f"6 {fields}",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not result.stdout.strip():
        raise ValueError("BLASTP returned no hit")
    return result.stdout.splitlines()[0].split("\t")


def main() -> None:
    # Reconstruct, validate and substitute the whale-shark BGN-like protein while
    # writing all supporting coordinates and BLAST fields to an audit table.
    protein, seqid, strand, cds_count, coding_length = reconstruct_protein()
    original = read_fasta(TARGET_FILES[0])
    human = next(seq for header, seq in original if header.startswith("human|BGN|"))
    references = read_fasta(REFERENCE)

    with tempfile.TemporaryDirectory(prefix="bgn_whale_repair_") as temp_name:
        temp = Path(temp_name)
        human_fasta = temp / "human_bgn.faa"
        target_fasta = temp / "whale_bgn.faa"
        family_fasta = temp / "family.faa"
        write_fasta(human_fasta, [("human_BGN", human)])
        write_fasta(target_fasta, [(TARGET_ACCESSION, protein)])
        family_records = []
        for header, sequence in references:
            gene, _class_name, accession = header.split()[0].split("|")
            family_records.append((f"{gene}__{accession}", sequence))
        write_fasta(family_fasta, family_records)
        pident, qcovs, alignment_length, evalue, bitscore = run_blast(
            human_fasta, target_fasta, "pident qcovs length evalue bitscore"
        )
        reverse_id, reverse_pident, reverse_qcovs, reverse_evalue, reverse_bitscore = (
            run_blast(target_fasta, family_fasta, "sseqid pident qcovs evalue bitscore")
        )
    reverse_gene, reverse_accession = reverse_id.split("__", 1)
    if reverse_gene != "BGN":
        raise ValueError(f"Correct GFF model reverses to {reverse_gene}, not BGN")

    new_header = (
        f"whale_shark|BGN|{TARGET_ACCESSION} "
        f"forward_pident={float(pident):.2f} "
        f"forward_qcovs={float(qcovs):.2f} reverse_hit={reverse_accession}"
    )
    for path in TARGET_FILES:
        records = read_fasta(path)
        matches = [
            i
            for i, (header, _sequence) in enumerate(records)
            if header.startswith("whale_shark|BGN|")
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one whale-shark BGN record in {path}")
        records[matches[0]] = (new_header, protein)
        write_fasta(path, records)

    row = {
        "gene": "BGN",
        "species": "whale_shark",
        "removed_mislabeled_accession": WRONG_ACCESSION,
        "replacement_accession": TARGET_ACCESSION,
        "replacement_status": "tentative_partial_biglycan_like_model",
        "seqid": seqid,
        "strand": strand,
        "cds_block_count": cds_count,
        "cds_length_bp": coding_length,
        "protein_length_aa": len(protein),
        "forward_identity_percent": f"{float(pident):.2f}",
        "forward_query_coverage_percent": f"{float(qcovs):.2f}",
        "forward_alignment_length": alignment_length,
        "forward_evalue": evalue,
        "forward_bitscore": bitscore,
        "reverse_best_gene": reverse_gene,
        "reverse_best_accession": reverse_accession,
        "reverse_identity_percent": f"{float(reverse_pident):.2f}",
        "reverse_query_coverage_percent": f"{float(reverse_qcovs):.2f}",
        "reverse_evalue": reverse_evalue,
        "reverse_bitscore": reverse_bitscore,
        "evidence": (
            "NCBI GFF3 assigns XP_048463588.1 to dcn; XP_048475905.1 is the "
            "partial biglycan-like model at the SynVoy BGN locus"
        ),
    }
    with AUDIT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(
        f"Replaced {WRONG_ACCESSION} with {TARGET_ACCESSION} in {len(TARGET_FILES)} BGN FASTAs"
    )
    print(
        f"Reconstructed {len(protein)} aa from {cds_count} CDS blocks; reverse best hit {reverse_gene}"
    )
    print(f"Wrote {AUDIT}")


if __name__ == "__main__":
    main()
