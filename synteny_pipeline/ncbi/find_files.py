from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GenomeFiles:
    genome_fasta: Path | None
    gff3: Path | None
    protein_faa: Path | None
    assembly_report: Path | None
    dataset_catalog: Path | None


def find_genome_files(unzipped_dir: Path) -> GenomeFiles:
    """Locate genome, annotation, protein, and metadata files in one NCBI datasets folder."""
    unzipped_dir = unzipped_dir.resolve()
    if not unzipped_dir.exists():
        raise FileNotFoundError(f"Unzipped directory does not exist: {unzipped_dir}")

    root = prefer_ncbi_dataset_root(unzipped_dir)

    fasta_candidates = sorted(
        find_by_suffix(root, {".fna", ".fa", ".fasta"}),
        key=lambda path: (not looks_like_genome_fasta(path.name), len(str(path))),
    )
    if not fasta_candidates:
        raise FileNotFoundError(f"No genome FASTA found under {root}")

    gff_candidates = sorted(
        find_by_suffix(root, {".gff", ".gff3"}),
        key=lambda path: (not looks_like_annotation_gff(path.name), len(str(path))),
    )

    protein_candidates = sorted(
        find_by_suffix(root, {".faa", ".fa", ".fasta"}),
        key=lambda path: (not looks_like_protein_fasta(path.name), len(str(path))),
    )

    return GenomeFiles(
        genome_fasta=fasta_candidates[0],
        gff3=gff_candidates[0] if gff_candidates else None,
        protein_faa=protein_candidates[0] if protein_candidates else None,
        assembly_report=first_match(root, ["assembly_report.txt", "assembly_report_"]),
        dataset_catalog=first_match(root, ["dataset_catalog.json"]),
    )


def prefer_ncbi_dataset_root(unzipped_dir: Path) -> Path:
    candidate = unzipped_dir / "ncbi_dataset"
    return candidate if candidate.exists() else unzipped_dir


def find_by_suffix(root: Path, suffixes: set[str]) -> list[Path]:
    return [
        path
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    ]


def looks_like_genome_fasta(name: str) -> bool:
    name = name.lower()
    if any(
        x in name
        for x in ["cds", "rna", "transcript", "protein", "proteins", "prot", "pep"]
    ):
        return False
    return any(
        x in name for x in ["genomic", "genome", "chromosome", "assembly"]
    ) or name.endswith(".fna")


def looks_like_annotation_gff(name: str) -> bool:
    name = name.lower()
    if any(x in name for x in ["cds", "rna", "protein", "pep"]):
        return False
    return any(x in name for x in ["genomic", "annotation", "genes", ".gff", ".gff3"])


def looks_like_protein_fasta(name: str) -> bool:
    name = name.lower()
    if any(x in name for x in ["genomic", "genome", "cds", "rna", "transcript"]):
        return False
    return any(x in name for x in ["protein", "proteins", "prot", "pep", ".faa"])


def first_match(root: Path, patterns: list[str]) -> Path | None:
    patterns = [pattern.lower() for pattern in patterns]
    for path in root.rglob("*"):
        if path.is_file() and any(pattern in path.name.lower() for pattern in patterns):
            return path
    return None
