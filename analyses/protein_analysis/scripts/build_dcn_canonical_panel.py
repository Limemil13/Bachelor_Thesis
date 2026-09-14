"""Build the compact DCN ortholog panel from the pinned NCBI ortholog download.

The panel mirrors the vertebrate sampling used for the other thesis SLRPs. An
amphioxus sequence is intentionally not forced into the set: the old DCN hit
XP_066272062.1 is also the historical LUM hit and is not a confident DCN
ortholog. Each retained protein is checked against the human DCN query and a
human SLRP reference panel with BLASTP.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
PROJECT_ROOT = BASE.parents[1]
SOURCE_FASTA = (
    PROJECT_ROOT
    / "analyses/phylogenetics/per_gene_trees/DCN/ncbi_download/DCN_orthologs"
    / "ncbi_dataset/data/protein.faa"
)
REFERENCE_FASTA = (
    PROJECT_ROOT
    / "analyses/phylogenetics/combined_trees/slrp_family_reference_tree"
    / "SLRP_reference_context.faa"
)
OUT_FASTA = BASE / "candidates/DCN_domain_input.faa"
OUT_SELECTION = BASE / "candidates/DCN_candidate_selection.tsv"


@dataclass(frozen=True)
class Choice:
    species: str
    organism: str
    accession: str
    status: str = "retained"
    note: str = "NCBI DCN ortholog-group representative"


CHOICES = (
    Choice(
        "human", "Homo sapiens", "NP_001911.1", note="RefSeq isoform a; human reference"
    ),
    Choice("mouse", "Mus musculus", "NP_031859.1"),
    Choice("rattus_norvegicus", "Rattus norvegicus", "NP_077043.1"),
    Choice("bos_taurus", "Bos taurus", "NP_776331.2"),
    Choice("canis_lupus_familiaris", "Canis lupus familiaris", "NP_001003228.1"),
    Choice("chicken", "Gallus gallus", "NP_001025918.2"),
    Choice("anolis_carolinensis", "Anolis carolinensis", "XP_003221092.1"),
    Choice("frog", "Xenopus tropicalis", "NP_001093704.1"),
    Choice("zebrafish", "Danio rerio", "NP_571772.1"),
    Choice("lepisosteus_oculatus", "Lepisosteus oculatus", "XP_015208366.1"),
    Choice("latimeria_chalumnae", "Latimeria chalumnae", "XP_006004204.1"),
    Choice(
        "monodelphis_domestica",
        "Monodelphis domestica",
        "XP_001363160.3",
        status="tentative_partial",
        note="Current NCBI DCN ortholog is only 212 aa; inspect against SynVoy locus model",
    ),
    Choice("callorhinchus_milii", "Callorhinchus milii", "NP_001279248.1"),
    Choice(
        "scyliorhinus_canicula",
        "Scyliorhinus canicula",
        "XP_038636759.1",
        status="tentative_isoform",
        note="Longest annotated X1 isoform; compare with 350-aa X2 during MSA review",
    ),
    Choice(
        "whale_shark",
        "Rhincodon typus",
        "XP_048463588.1",
        note="Same protein sequence as alternate accession XP_048463587.1",
    ),
)


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(chunks))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def run_blast(query: Path, subject: Path, fields: str) -> list[list[str]]:
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
    return [line.split("\t") for line in result.stdout.splitlines() if line.strip()]


def main() -> None:
    if not SOURCE_FASTA.exists():
        raise FileNotFoundError(SOURCE_FASTA)
    if not REFERENCE_FASTA.exists():
        raise FileNotFoundError(REFERENCE_FASTA)

    source = read_fasta(SOURCE_FASTA)
    references = read_fasta(REFERENCE_FASTA)
    missing = [choice.accession for choice in CHOICES if choice.accession not in source]
    if missing:
        raise ValueError(f"Selected DCN proteins absent from NCBI download: {missing}")

    human_sequence = source["NP_001911.1"][1]
    family_records: list[tuple[str, str]] = []
    for header, sequence in references.values():
        parts = header.split()[0].split("|")
        if len(parts) != 3:
            raise ValueError(f"Unexpected SLRP reference header: {header}")
        gene, _, accession = parts
        family_records.append((f"{gene}__{accession}", sequence))

    output_records: list[tuple[str, str]] = []
    selection_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="dcn_panel_") as temp_name:
        temp = Path(temp_name)
        human_query = temp / "human_dcn.faa"
        family_fasta = temp / "human_slrp_references.faa"
        write_fasta(human_query, [("human_DCN_NP_001911.1", human_sequence)])
        write_fasta(family_fasta, family_records)

        for choice in CHOICES:
            original_header, sequence = source[choice.accession]
            expected_organism = f"[organism={choice.organism}]"
            if expected_organism not in original_header:
                raise ValueError(
                    f"{choice.accession} has wrong organism header: {original_header}"
                )

            candidate = temp / f"{choice.accession}.faa"
            write_fasta(candidate, [(choice.accession, sequence)])
            forward = run_blast(
                human_query,
                candidate,
                "pident qcovs length evalue bitscore",
            )
            reverse = run_blast(
                candidate,
                family_fasta,
                "sseqid pident qcovs evalue bitscore",
            )
            if not forward or not reverse:
                raise ValueError(f"BLASTP returned no result for {choice.accession}")

            pident, qcovs, aln_len, evalue, bitscore = forward[0]
            (
                reverse_id,
                reverse_pident,
                reverse_qcovs,
                reverse_evalue,
                reverse_bitscore,
            ) = reverse[0]
            reverse_gene, reverse_accession = reverse_id.split("__", 1)
            if reverse_gene != "DCN":
                raise ValueError(
                    f"{choice.accession} reverse family hit is {reverse_gene}, not DCN"
                )

            header = (
                f"{choice.species}|DCN|{choice.accession} "
                f"forward_pident={float(pident):.2f} "
                f"forward_qcovs={float(qcovs):.2f} "
                f"reverse_hit={reverse_accession}"
            )
            output_records.append((header, sequence))
            selection_rows.append(
                {
                    "gene": "DCN",
                    "species": choice.species,
                    "organism": choice.organism,
                    "accession": choice.accession,
                    "protein_length": len(sequence),
                    "selection_status": choice.status,
                    "forward_identity_percent": f"{float(pident):.2f}",
                    "forward_query_coverage_percent": f"{float(qcovs):.2f}",
                    "forward_alignment_length": aln_len,
                    "forward_evalue": evalue,
                    "forward_bitscore": bitscore,
                    "reverse_best_gene": reverse_gene,
                    "reverse_best_accession": reverse_accession,
                    "reverse_identity_percent": f"{float(reverse_pident):.2f}",
                    "reverse_query_coverage_percent": f"{float(reverse_qcovs):.2f}",
                    "reverse_evalue": reverse_evalue,
                    "reverse_bitscore": reverse_bitscore,
                    "notes": choice.note,
                    "original_header": original_header,
                }
            )

    if len(output_records) != 15 or len({row[0] for row in output_records}) != 15:
        raise ValueError("DCN panel must contain 15 unique proteins")

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(OUT_FASTA, output_records)
    fields = list(selection_rows[0])
    with OUT_SELECTION.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selection_rows)

    print(f"Wrote {OUT_FASTA} ({len(output_records)} proteins)")
    print(f"Wrote {OUT_SELECTION}")


if __name__ == "__main__":
    main()
