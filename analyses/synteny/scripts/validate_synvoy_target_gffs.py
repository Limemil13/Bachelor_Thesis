"""Validate that SynVoy target FASTAs have matching, nonempty GFF annotations.

The check compares every FASTA sequence accession used by an annotated GFF
record. Checking only the first chromosome can miss a mixture of accession
versions within an otherwise correctly named assembly.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def fasta_accessions(path: Path) -> set[str]:
    # Collect all assembly sequence accessions present in the genome FASTA.
    accessions: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                accessions.add(line[1:].split(maxsplit=1)[0])
    if not accessions:
        raise ValueError(f"No FASTA header found in {path}")
    return accessions


def annotated_gff_seqids(path: Path) -> set[str]:
    # Collect sequence IDs from actual GFF features, not only header directives.
    seqids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            if len(fields) >= 3 and fields[2] in {"gene", "mRNA", "CDS"}:
                seqids.add(fields[0])
    return seqids


def main() -> None:
    # Validate each named FASTA/GFF pair before launching an expensive SynVoy run.
    parser = argparse.ArgumentParser()
    parser.add_argument("--fna-dir", required=True, type=Path)
    parser.add_argument("--gff-dir", required=True, type=Path)
    parser.add_argument("--expected", type=int, default=14)
    args = parser.parse_args()

    fastas = sorted(args.fna_dir.glob("*.fna"))
    failures: list[str] = []
    for fasta in fastas:
        fasta_seqids = fasta_accessions(fasta)
        gff = args.gff_dir / f"{fasta.stem}.gff"
        if not gff.is_file() or gff.stat().st_size == 0:
            failures.append(f"{fasta.stem}: missing or empty {gff}")
            print(f"FAIL\t{fasta.stem}\tmissing_gff")
            continue
        gff_seqids = annotated_gff_seqids(gff)
        missing = sorted(gff_seqids - fasta_seqids)
        if missing:
            sample = ",".join(missing[:10])
            failures.append(
                f"{fasta.stem}: {len(missing)} annotated GFF seqid(s) absent "
                f"from FASTA ({sample})"
            )
            print(
                f"FAIL\t{fasta.stem}\t{len(fasta_seqids)}\t{len(gff_seqids)}\t"
                f"{len(missing)}\taccession_mismatch"
            )
            continue
        print(
            f"PASS\t{fasta.stem}\t{len(fasta_seqids)}\t{len(gff_seqids)}\t"
            f"0\t{gff.stat().st_size}"
        )

    if len(fastas) != args.expected:
        failures.append(f"expected {args.expected} FASTAs, found {len(fastas)}")

    if failures:
        raise SystemExit("Validation failed:\n- " + "\n- ".join(failures))
    print(f"VALIDATED\t{len(fastas)}/{args.expected}")


if __name__ == "__main__":
    main()
