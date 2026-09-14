#!/usr/bin/env python3
"""Reconstruct one annotated protein from an indexed genome FASTA and GFF3."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.parse import unquote

from Bio.Seq import Seq


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = unquote(value)
    return attributes


class IndexedFasta:
    def __init__(self, fasta: Path, needed_seqids: set[str]):
        self.handle = fasta.open("rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        self.sequences: dict[str, str] = {}
        index_path = fasta.with_suffix(fasta.suffix + ".fai")
        if index_path.is_file():
            with index_path.open(encoding="utf-8") as index_handle:
                for line in index_handle:
                    name, length, offset, line_bases, line_width, *_ = line.split("\t")
                    self.index[name] = (
                        int(length),
                        int(offset),
                        int(line_bases),
                        int(line_width),
                    )
        else:
            current = ""
            parts: list[str] = []
            with fasta.open(encoding="ascii") as fasta_handle:
                for line in fasta_handle:
                    if line.startswith(">"):
                        if current in needed_seqids:
                            self.sequences[current] = "".join(parts).upper()
                        current = line[1:].split(maxsplit=1)[0]
                        parts = []
                    elif current in needed_seqids:
                        parts.append(line.strip())
                if current in needed_seqids:
                    self.sequences[current] = "".join(parts).upper()
            missing = needed_seqids - set(self.sequences)
            if missing:
                raise ValueError(f"FASTA sequences not found: {sorted(missing)}")

    def fetch(self, name: str, start: int, end: int) -> str:
        if name in self.sequences:
            sequence = self.sequences[name]
            if start < 1 or end > len(sequence) or start > end:
                raise ValueError(f"Invalid interval {name}:{start}-{end}")
            return sequence[start - 1 : end]
        length, offset, line_bases, line_width = self.index[name]
        if start < 1 or end > length or start > end:
            raise ValueError(f"Invalid interval {name}:{start}-{end}")
        position = start - 1
        remaining = end - start + 1
        output = bytearray()
        while remaining:
            line_position = position % line_bases
            take = min(remaining, line_bases - line_position)
            byte_position = (
                offset + (position // line_bases) * line_width + line_position
            )
            self.handle.seek(byte_position)
            output.extend(self.handle.read(take))
            position += take
            remaining -= take
        return output.decode("ascii").upper()

    def close(self) -> None:
        self.handle.close()


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--protein-id", required=True)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-qc", required=True, type=Path)
    args = parser.parse_args()

    cds: list[dict[str, object]] = []
    with args.gff.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attributes = parse_attributes(fields[8])
            if attributes.get("protein_id") != args.protein_id:
                continue
            cds.append(
                {
                    "seqid": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "phase": fields[7],
                    "gene": attributes.get("gene", ""),
                    "product": attributes.get("product", ""),
                    "parent": attributes.get("Parent", ""),
                }
            )
    if not cds:
        raise ValueError(f"Protein {args.protein_id} was not found in {args.gff}")
    if len({str(row["seqid"]) for row in cds}) != 1:
        raise ValueError("CDS records occur on multiple sequences")
    if len({str(row["strand"]) for row in cds}) != 1:
        raise ValueError("CDS records use multiple strands")

    strand = str(cds[0]["strand"])
    cds.sort(key=lambda row: int(row["start"]), reverse=strand == "-")
    fasta = IndexedFasta(args.fasta, {str(cds[0]["seqid"])})
    try:
        parts: list[str] = []
        for row in cds:
            part = fasta.fetch(
                str(row["seqid"]), int(row["start"]), int(row["end"])
            )
            parts.append(
                part.translate(COMPLEMENT)[::-1] if strand == "-" else part
            )
    finally:
        fasta.close()

    coding_sequence = "".join(parts)
    translation = str(Seq(coding_sequence).translate())
    protein = translation.removesuffix("*")
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_fasta.write_text(
        f">{args.protein_id} gene={cds[0]['gene']} product={cds[0]['product']}\n"
        f"{wrap(protein)}\n",
        encoding="utf-8",
    )

    qc = {
        "protein_id": args.protein_id,
        "gene": str(cds[0]["gene"]),
        "product": str(cds[0]["product"]),
        "seqid": str(cds[0]["seqid"]),
        "strand": strand,
        "cds_blocks": len(cds),
        "cds_length_bp": len(coding_sequence),
        "length_divisible_by_three": "yes" if len(coding_sequence) % 3 == 0 else "no",
        "protein_length_aa": len(protein),
        "starts_with_methionine": "yes" if protein.startswith("M") else "no",
        "terminal_stop_present": "yes" if translation.endswith("*") else "no",
        "internal_stop_count": translation[:-1].count("*"),
        "gff_phases_transcript_order": ",".join(str(row["phase"]) for row in cds),
        "cds_intervals_transcript_order": ",".join(
            f"{row['seqid']}:{row['start']}-{row['end']}({strand})" for row in cds
        ),
    }
    args.output_qc.parent.mkdir(parents=True, exist_ok=True)
    with args.output_qc.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc), delimiter="\t")
        writer.writeheader()
        writer.writerow(qc)
    print(f"Reconstructed {args.protein_id}: {len(protein)} aa from {len(cds)} CDS blocks")


if __name__ == "__main__":
    main()
