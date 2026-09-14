"""Build and validate the canonical compact SLRP protein dataset.

The compact thesis dataset contains one selected protein per confident
gene/species key for BGN, DCN, EPYC, FMOD, OGN, PRELP, and LUM. Dog OGN is read
from the corrected FASTA; amphioxus OGN is retained only as manual-review
provenance. DCN is the promoted class-I comparator; its questionable amphioxus
hit is deliberately excluded from the compact ortholog set.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
PROJECT_ROOT = BASE.parents[1]
CANDIDATE_DIR = BASE / "candidates"

GENE_ORDER = ("BGN", "DCN", "EPYC", "FMOD", "OGN", "PRELP", "LUM")
EXPECTED_COUNTS = {
    "BGN": 13,
    "DCN": 15,
    "EPYC": 15,
    "FMOD": 14,
    "OGN": 14,
    "PRELP": 16,
    "LUM": 16,
}

CANONICAL_FASTAS = {
    "BGN": CANDIDATE_DIR / "BGN_domain_input_canonical_decontaminated.faa",
    "DCN": CANDIDATE_DIR / "DCN_domain_input.faa",
    "EPYC": CANDIDATE_DIR / "EPYC_domain_input.faa",
    "FMOD": CANDIDATE_DIR / "FMOD_domain_input_lamprey_deduplicated.faa",
    "OGN": CANDIDATE_DIR / "OGN_domain_input_dog_corrected.faa",
    "PRELP": CANDIDATE_DIR / "PRELP_domain_input.faa",
    "LUM": CANDIDATE_DIR / "LUM_domain_input.faa",
}

ORIGINAL_OGN_FASTA = CANDIDATE_DIR / "OGN_domain_input.faa"
COMBINED_FASTA = CANDIDATE_DIR / "slrp_candidates_canonical.faa"
MANIFEST = CANDIDATE_DIR / "canonical_candidate_manifest.tsv"
FILE_SELECTION = CANDIDATE_DIR / "canonical_file_selection.tsv"
BY_GENE_DIR = CANDIDATE_DIR / "canonical_by_gene"

MANUAL_EXCLUDED_KEYS = {
    ("OGN", "amphioxus", "XP_066265713.1"),
}

METRIC_RE = re.compile(r"\b(forward_pident|forward_qcovs|reverse_hit)=([^\s]+)")


@dataclass(frozen=True)
class Record:
    gene: str
    species: str
    accession: str
    header: str
    sequence: str
    source_file: Path
    original_length: int
    correction: str
    forward_pident: str
    forward_qcovs: str
    reverse_hit: str

    @property
    def key(self) -> tuple[str, str, str]:
        return self.gene, self.species, self.accession

    @property
    def canonical_id(self) -> str:
        return "|".join(self.key)


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


def parse_header(header: str) -> tuple[str, str, str, dict[str, str]]:
    first = header.split()[0]
    parts = first.split("|")
    if len(parts) < 3:
        raise ValueError(f"Expected species|gene|accession header, got: {header}")
    species, gene, accession = parts[:3]
    metrics = dict(METRIC_RE.findall(header))
    return gene, species, accession, metrics


def project_relative(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()


def load_records() -> list[Record]:
    missing = [path for path in CANONICAL_FASTAS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing canonical FASTA(s): " + ", ".join(map(str, missing))
        )

    original_ogn_lengths: dict[tuple[str, str, str], int] = {}
    for header, sequence in read_fasta(ORIGINAL_OGN_FASTA):
        gene, species, accession, _ = parse_header(header)
        original_ogn_lengths[(gene, species, accession)] = len(sequence)

    records: list[Record] = []
    for expected_gene in GENE_ORDER:
        source = CANONICAL_FASTAS[expected_gene]
        for header, sequence in read_fasta(source):
            gene, species, accession, metrics = parse_header(header)
            if gene != expected_gene:
                raise ValueError(f"{source} contains {gene}, expected {expected_gene}")

            key = (gene, species, accession)
            if key in MANUAL_EXCLUDED_KEYS:
                continue
            original_length = original_ogn_lengths.get(key, len(sequence))
            correction = "none"
            if key == ("OGN", "canis_lupus_familiaris", "XP_038383338.1"):
                correction = "trimmed_first_64_aa_from_MKTLQST"
                if original_length - len(sequence) != 64:
                    raise ValueError("Corrected dog OGN must be exactly 64 aa shorter")

            records.append(
                Record(
                    gene=gene,
                    species=species,
                    accession=accession,
                    header=header,
                    sequence=sequence,
                    source_file=source,
                    original_length=original_length,
                    correction=correction,
                    forward_pident=metrics.get("forward_pident", ""),
                    forward_qcovs=metrics.get("forward_qcovs", ""),
                    reverse_hit=metrics.get("reverse_hit", ""),
                )
            )

    keys = [record.key for record in records]
    duplicate_keys = [key for key, count in Counter(keys).items() if count > 1]
    if duplicate_keys:
        raise ValueError(f"Duplicate canonical protein keys: {duplicate_keys}")

    counts = Counter(record.gene for record in records)
    if dict(counts) != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected per-gene counts: {dict(counts)}")
    expected_total = sum(EXPECTED_COUNTS.values())
    if len(records) != expected_total:
        raise ValueError(
            f"Expected {expected_total} canonical proteins, found {len(records)}"
        )
    return records


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    records = load_records()

    with COMBINED_FASTA.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.canonical_id}\n")
            for start in range(0, len(record.sequence), 60):
                handle.write(record.sequence[start : start + 60] + "\n")

    BY_GENE_DIR.mkdir(parents=True, exist_ok=True)
    for gene in GENE_ORDER:
        with (BY_GENE_DIR / f"{gene}_canonical.faa").open(
            "w", encoding="utf-8"
        ) as handle:
            for record in records:
                if record.gene != gene:
                    continue
                handle.write(f">{record.header}\n")
                for start in range(0, len(record.sequence), 60):
                    handle.write(record.sequence[start : start + 60] + "\n")

    manifest_rows: list[dict[str, object]] = []
    for record in records:
        manifest_rows.append(
            {
                "gene": record.gene,
                "species": record.species,
                "accession": record.accession,
                "canonical_id": record.canonical_id,
                "original_length": record.original_length,
                "analyzed_length": len(record.sequence),
                "length_delta": len(record.sequence) - record.original_length,
                "sequence_correction": record.correction,
                "forward_pident": record.forward_pident,
                "forward_qcovs": record.forward_qcovs,
                "reverse_hit": record.reverse_hit,
                "source_file": project_relative(record.source_file),
                "original_header": record.header,
            }
        )

    manifest_fields = [
        "gene",
        "species",
        "accession",
        "canonical_id",
        "original_length",
        "analyzed_length",
        "length_delta",
        "sequence_correction",
        "forward_pident",
        "forward_qcovs",
        "reverse_hit",
        "source_file",
        "original_header",
    ]
    write_tsv(MANIFEST, manifest_rows, manifest_fields)

    selection_rows = [
        {
            "gene": gene,
            "canonical_source_file": project_relative(CANONICAL_FASTAS[gene]),
            "sequence_count": EXPECTED_COUNTS[gene],
            "status": "canonical",
        }
        for gene in GENE_ORDER
    ]
    write_tsv(
        FILE_SELECTION,
        selection_rows,
        ["gene", "canonical_source_file", "sequence_count", "status"],
    )

    print(f"Wrote {COMBINED_FASTA}")
    print(f"Wrote {MANIFEST}")
    print(f"Wrote {FILE_SELECTION}")
    print(f"Canonical proteins: {len(records)}")


if __name__ == "__main__":
    main()
