#!/usr/bin/env python3
"""Extract strand-aware promoter windows for representative SLRP transcripts.

The representative-transcript table is the same source used by the local gene-
structure analysis.  Promoters are therefore transcript-choice dependent and
must not be interpreted as experimentally validated regulatory elements.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
WINDOWS = {
    "core": (500, 100),
    "proximal": (2000, 200),
}
COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def ensure_fai(samtools: str, fasta: Path) -> dict[str, int]:
    fai = Path(f"{fasta}.fai")
    if not fai.is_file():
        subprocess.run([samtools, "faidx", str(fasta)], check=True)
    lengths: dict[str, int] = {}
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def fetch(samtools: str, fasta: Path, seqid: str, start: int, end: int) -> str:
    output = subprocess.check_output(
        [samtools, "faidx", "-n", "1000000000", str(fasta), f"{seqid}:{start}-{end}"],
        text=True,
    )
    lines = output.splitlines()
    return "".join(lines[1:]).upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument("--genome-dir", required=True, type=Path)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    rows = [
        row
        for row in read_tsv(args.representatives)
        if row["query_gene"] in GENES and row["species_or_file"] in SPECIES
    ]
    if len(rows) != len(GENES) * len(SPECIES):
        raise SystemExit(
            f"Expected {len(GENES) * len(SPECIES)} representative rows, found {len(rows)}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequence_sets: dict[str, dict[str, str]] = defaultdict(dict)
    metadata: list[dict[str, str | int | float]] = []
    fai_by_species: dict[str, dict[str, int]] = {}

    for species in SPECIES:
        fasta = args.genome_dir / f"{species}.fna"
        fai_by_species[species] = ensure_fai(args.samtools, fasta)

    for row in sorted(
        rows, key=lambda item: (item["query_gene"], item["species_or_file"])
    ):
        gene = row["query_gene"]
        species = row["species_or_file"]
        seqid = row["seqid"]
        strand = row["strand"]
        transcript = row["transcript_id"]
        transcript_start = int(row["gene_start"])
        transcript_end = int(row["gene_end"])
        tss = transcript_start if strand == "+" else transcript_end
        chromosome_length = fai_by_species[species][seqid]
        fasta = args.genome_dir / f"{species}.fna"

        for window_name, (upstream, downstream) in WINDOWS.items():
            if strand == "+":
                raw_start, raw_end = tss - upstream, tss + downstream
            else:
                raw_start, raw_end = tss - downstream, tss + upstream
            start = max(1, raw_start)
            end = min(chromosome_length, raw_end)
            sequence = fetch(args.samtools, fasta, seqid, start, end)
            if strand == "-":
                sequence = reverse_complement(sequence)
            expected_length = upstream + downstream + 1
            if len(sequence) != end - start + 1:
                raise RuntimeError(
                    f"Length mismatch for {gene}/{species}/{window_name}"
                )
            header = (
                f"{gene}|{species}|{transcript}|{window_name}|"
                f"{seqid}:{start}-{end}|strand={strand}|tss={tss}"
            )
            sequence_sets[window_name][header] = sequence
            ambiguous = sum(base not in {"A", "C", "G", "T"} for base in sequence)
            metadata.append(
                {
                    "gene": gene,
                    "species": species,
                    "transcript_id": transcript,
                    "seqid": seqid,
                    "strand": strand,
                    "tss_1based": tss,
                    "window": window_name,
                    "upstream_bp": upstream,
                    "downstream_bp": downstream,
                    "genomic_start_1based": start,
                    "genomic_end_1based": end,
                    "sequence_length_bp": len(sequence),
                    "expected_length_bp": expected_length,
                    "clipped_at_contig_edge": "yes"
                    if len(sequence) != expected_length
                    else "no",
                    "gc_percent": round(
                        100
                        * (sequence.count("G") + sequence.count("C"))
                        / len(sequence),
                        3,
                    ),
                    "ambiguous_percent": round(100 * ambiguous / len(sequence), 3),
                    "header": header,
                    "source_fasta": str(fasta),
                    "representative_source": str(args.representatives),
                }
            )

    for window_name, sequences in sequence_sets.items():
        combined = args.output_dir / f"slrp_{window_name}_promoters_5species.fna"
        with combined.open("w", encoding="utf-8") as handle:
            for header, sequence in sequences.items():
                handle.write(f">{header}\n{wrap(sequence)}\n")

        per_gene = args.output_dir / "per_gene" / window_name
        per_gene.mkdir(parents=True, exist_ok=True)
        for gene in GENES:
            positives = {h: s for h, s in sequences.items() if h.startswith(f"{gene}|")}
            controls = {
                h: s for h, s in sequences.items() if not h.startswith(f"{gene}|")
            }
            for label, selected in (
                ("positive", positives),
                ("other_genes_control", controls),
            ):
                path = per_gene / f"{gene}_{label}.fna"
                with path.open("w", encoding="utf-8") as handle:
                    for header, sequence in selected.items():
                        handle.write(f">{header}\n{wrap(sequence)}\n")

    fields = list(metadata[0])
    with (args.output_dir / "promoter_metadata.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata)

    print(
        f"Extracted {len(metadata)} promoter windows from {len(rows)} representative transcripts"
    )


if __name__ == "__main__":
    main()
