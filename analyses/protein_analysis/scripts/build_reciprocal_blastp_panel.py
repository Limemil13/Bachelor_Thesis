#!/usr/bin/env python3
"""homo sapiens as query -> against all proteomes to find best sequence, and reverse blasting
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import textwrap
from pathlib import Path

DEFAULT_SPECIES = (
    "human",
    "mouse",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "bos_taurus",
    "monodelphis_domestica",
    "chicken",
    "anolis_carolinensis",
    "frog",
    "zebrafish",
    "lepisosteus_oculatus",
    "latimeria_chalumnae",
    "callorhinchus_milii",
    "scyliorhinus_canicula",
    "lamprey",
    "amphioxus",
    "whale_shark",
)


def normalize_accession(value: str) -> str:
    token = value.split()[0]
    if "|" in token:
        parts = [part for part in token.split("|") if part]
        token = parts[-1]
    return token


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def find_record(path: Path, accession: str) -> tuple[str, str] | None:
    wanted = normalize_accession(accession)
    wanted_base = wanted.split(".")[0]
    for header, sequence in read_fasta(path):
        observed = normalize_accession(header)
        if observed == wanted or observed.split(".")[0] == wanted_base:
            return header, sequence
    return None


def ensure_blast_db(proteome: Path, prefix: Path) -> None:
    #keep search local for faster reruns
    if any(prefix.parent.glob(prefix.name + ".*")):
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["makeblastdb", "-in", str(proteome), "-dbtype", "prot", "-out", str(prefix)],
        check=True,
    )


def blastp_top1(query: Path, db_prefix: Path, evalue: float) -> dict[str, str] | None:
    command = [
        "blastp",
        "-query",
        str(query),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length evalue bitscore qcovs",
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    if not result.stdout.strip():
        return None
    columns = result.stdout.splitlines()[0].split("\t")
    return dict(
        zip(
            ("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs"),
            columns,
            strict=True,
        )
    )


def write_fasta_record(path: Path, header: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f">{header}\n" + "\n".join(textwrap.wrap(sequence, 80)) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--species", default=",".join(DEFAULT_SPECIES))
    parser.add_argument("--evalue", type=float, default=1e-5)
    args = parser.parse_args()

    root = args.project_root.resolve()
    query = args.query.resolve()
    output_dir = args.output_dir.resolve()
    species_names = tuple(
        item.strip() for item in args.species.split(",") if item.strip()
    )
    query_records = read_fasta(query)
    if len(query_records) != 1:
        raise ValueError(
            f"Expected exactly one query sequence in {query}, found {len(query_records)}"
        )
    query_accession = normalize_accession(query_records[0][0])

    human_proteome = root / "data" / "ncbi" / "proteomes" / "human" / "protein.faa"
    human_db = human_proteome.parent / "blastdb" / "proteome_prot"
    ensure_blast_db(human_proteome, human_db)

    best_hit_dir = output_dir / "target_best_hits"
    rows: list[dict[str, object]] = []
    candidate_records: list[tuple[str, str]] = []
    for species in species_names:
        proteome = root / "data" / "ncbi" / "proteomes" / species / "protein.faa"
        if not proteome.exists():
            rows.append(
                {
                    "gene": args.gene,
                    "species": species,
                    "status": "TARGET_PROTEOME_MISSING",
                }
            )
            continue
        database = proteome.parent / "blastdb" / "proteome_prot"
        ensure_blast_db(proteome, database)
        forward = blastp_top1(query, database, args.evalue)
        if forward is None:
            rows.append(
                {"gene": args.gene, "species": species, "status": "NO_FORWARD_HIT"}
            )
            continue

        accession = normalize_accession(forward["sseqid"])
        target_record = find_record(proteome, accession)
        if target_record is None:
            rows.append(
                {
                    "gene": args.gene,
                    "species": species,
                    "forward_best_hit": accession,
                    "status": "TARGET_SEQUENCE_NOT_FOUND",
                }
            )
            continue
        target_header, target_sequence = target_record
        target_query = best_hit_dir / f"{species}_{accession}.faa"
        write_fasta_record(target_query, target_header, target_sequence)
        reverse = blastp_top1(target_query, human_db, args.evalue)
        reverse_accession = normalize_accession(reverse["sseqid"]) if reverse else ""
        reverse_record = (
            find_record(human_proteome, reverse_accession)
            if reverse_accession
            else None
        )
        reverse_description = reverse_record[0] if reverse_record else ""
        reciprocal_pass = (
            bool(reverse_accession)
            and reverse_accession.split(".")[0] == query_accession.split(".")[0]
        )

        rows.append(
            {
                "gene": args.gene,
                "query_accession": query_accession,
                "species": species,
                "forward_best_hit": accession,
                "forward_pident": forward["pident"],
                "forward_qcovs": forward["qcovs"],
                "forward_length": forward["length"],
                "forward_evalue": forward["evalue"],
                "forward_bitscore": forward["bitscore"],
                "target_protein_length": len(target_sequence),
                "reverse_best_hit": reverse_accession,
                "reverse_description": reverse_description,
                "reverse_pident": reverse["pident"] if reverse else "",
                "reverse_qcovs": reverse["qcovs"] if reverse else "",
                "reverse_evalue": reverse["evalue"] if reverse else "",
                "reciprocal_pass": "yes" if reciprocal_pass else "no",
                "status": "PASS" if reciprocal_pass else "REVERSE_BEST_DIFFERS",
                "target_original_header": target_header,
            }
        )
        candidate_header = (
            f"{species}|{args.gene}|{accession} "
            f"forward_pident={float(forward['pident']):.2f} "
            f"forward_qcovs={float(forward['qcovs']):.1f} "
            f"reverse_hit={reverse_accession}"
        )
        #keeping failed candidates incase of paralog
        candidate_records.append((candidate_header, target_sequence))

    output_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "gene",
        "query_accession",
        "species",
        "forward_best_hit",
        "forward_pident",
        "forward_qcovs",
        "forward_length",
        "forward_evalue",
        "forward_bitscore",
        "target_protein_length",
        "reverse_best_hit",
        "reverse_description",
        "reverse_pident",
        "reverse_qcovs",
        "reverse_evalue",
        "reciprocal_pass",
        "status",
        "target_original_header",
    ]
    with (output_dir / "reciprocal_blastp.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    with (output_dir / f"{args.gene}_candidates_all.faa").open(
        "w", encoding="utf-8"
    ) as handle:
        for header, sequence in candidate_records:
            handle.write(f">{header}\n")
            handle.write("\n".join(textwrap.wrap(sequence, 80)) + "\n")

    if len(rows) != len(species_names) or len(candidate_records) != len(species_names):
        raise RuntimeError(
            f"Incomplete screen: {len(rows)} rows and {len(candidate_records)} sequences "
            f"for {len(species_names)} species"
        )
    print(f"Wrote {output_dir / 'reciprocal_blastp.tsv'}")
    print(f"Wrote {output_dir / f'{args.gene}_candidates_all.faa'}")
    print(f"Species screened: {len(rows)}")
    print(
        f"Reciprocal passes: {sum(row.get('reciprocal_pass') == 'yes' for row in rows)}"
    )


if __name__ == "__main__":
    main()
