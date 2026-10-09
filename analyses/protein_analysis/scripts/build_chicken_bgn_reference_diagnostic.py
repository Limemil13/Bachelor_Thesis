#!/usr/bin/env python3
"""Add the flagged chicken BGN-locus protein to the human SLRP reference panel."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
REFERENCE = (
    ROOT
    / "analyses/phylogenetics/combined_trees/slrp_family_reference_tree/SLRP_reference_context.faa"
)
CANONICAL = ROOT / "analyses/protein_analysis/candidates/slrp_candidates_canonical.faa"
OUT = ROOT / "analyses/protein_analysis/diagnostics/BGN/chicken_annotation_conflict"
TARGET = "XP_414298.2"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    # Read reference and candidate FASTA records with complete headers.
    records = []
    header = None
    chunks = []
    for raw in path.read_text(encoding="utf-8").splitlines():
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
    # Write diagnostic subsets without changing canonical source files.
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    # Add the disputed chicken protein to the human SLRP reference panel so its
    # placement can be compared with BGN, ASPN and related paralogs.
    references = read_fasta(REFERENCE)
    matches = [record for record in read_fasta(CANONICAL) if TARGET in record[0]]
    if len(references) != 19 or len(matches) != 1:
        raise SystemExit(
            f"Expected 19 references and one {TARGET}; got {len(references)}, {len(matches)}"
        )
    chicken = ("CHICKEN_BGN_LOCUS|conflict|XP_414298.2", matches[0][1])
    OUT.mkdir(parents=True, exist_ok=True)
    write_fasta(
        OUT / "human_slrp_reference_plus_chicken_bgn_locus.faa", references + [chicken]
    )
    print(f"Wrote 20-sequence diagnostic input to {OUT}")


if __name__ == "__main__":
    main()
