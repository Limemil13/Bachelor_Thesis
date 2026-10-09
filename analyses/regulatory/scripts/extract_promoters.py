#!/usr/bin/env python3
"""Takes the already selected transcript for each gene and species, identifies its annotated
transcription start site, and extracts two surrounding DNA regions from the genome for motif analysis
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

GENES = ("BGN", "DCN", "FMOD", "PRELP", "EPYC", "LUM", "OGN")
SPECIES = ("human", "mouse", "cow", "chicken", "zebrafish")
EXCLUDED_LOCI = {
    ("BGN", "chicken"): {
        "expected_transcript_id": "rna-XM_414298.8",
        "reason": (
            "excluded from BGN promoter analysis because the translated product "
            "XP_414298.2 is ASPN-like and does not provide independent BGN "
            "orthology support"
        ),
    }
}
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


def portable_source(path: Path) -> str:
    """Keep an analysis-relative path instead of a local user directory."""
    parts = path.parts
    if "analyses" in parts:
        return Path(*parts[parts.index("analyses") :]).as_posix()
    return path.name


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

    selected_rows = [
        row
        for row in read_tsv(args.representatives)
        if row["query_gene"] in GENES and row["species_or_file"] in SPECIES
    ]
    expected_all = {(gene, species) for gene in GENES for species in SPECIES}
    selected_keys = [
        (row["query_gene"], row["species_or_file"]) for row in selected_rows
    ]
    if len(selected_keys) != len(set(selected_keys)):
        raise SystemExit("Representative table contains duplicate gene/species rows")
    if set(selected_keys) != expected_all:
        missing = sorted(expected_all - set(selected_keys))
        unexpected = sorted(set(selected_keys) - expected_all)
        raise SystemExit(
            "Representative table does not contain the expected 35-row panel; "
            f"missing={missing}, unexpected={unexpected}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    excluded_rows: list[dict[str, str]] = []
    rows: list[dict[str, str]] = []
    for row in selected_rows:
        key = (row["query_gene"], row["species_or_file"])
        exclusion = EXCLUDED_LOCI.get(key)
        if exclusion is None:
            rows.append(row)
            continue
        if row["transcript_id"] != exclusion["expected_transcript_id"]:
            raise SystemExit(
                f"Refusing to exclude unexpected transcript for {key}: "
                f"found {row['transcript_id']}, expected "
                f"{exclusion['expected_transcript_id']}"
            )
        excluded_rows.append(
            {
                "gene": row["query_gene"],
                "species": row["species_or_file"],
                "transcript_id": row["transcript_id"],
                "seqid": row["seqid"],
                "representative_note": row.get("note", ""),
                "exclusion_reason": exclusion["reason"],
                "representative_source": portable_source(args.representatives),
            }
        )

    expected_retained = expected_all - set(EXCLUDED_LOCI)
    retained_keys = {(row["query_gene"], row["species_or_file"]) for row in rows}
    if retained_keys != expected_retained:
        raise SystemExit(
            "Promoter exclusion produced an unexpected retained panel; "
            f"missing={sorted(expected_retained - retained_keys)}, "
            f"unexpected={sorted(retained_keys - expected_retained)}"
        )

    exclusion_fields = [
        "gene",
        "species",
        "transcript_id",
        "seqid",
        "representative_note",
        "exclusion_reason",
        "representative_source",
    ]
    with (args.output_dir / "promoter_excluded_loci.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=exclusion_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded_rows)

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
                    "source_fasta": f"<genome-dir>/{fasta.name}",
                    "representative_source": portable_source(args.representatives),
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
        f"Extracted {len(metadata)} promoter windows from {len(rows)} retained "
        f"representative transcripts; documented {len(excluded_rows)} exclusion"
    )


if __name__ == "__main__":
    main()
