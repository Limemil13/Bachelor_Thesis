#!/usr/bin/env python3
"""Create a targeted cartilage/growth-plate TF motif set from a MEME database."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

DEFAULT_TFS = (
    "SOX5",
    "SOX6",
    "SOX9",
    "RUNX2",
    "RUNX3",
    "SP7",
    "MEF2C",
    "HIF1A",
    "GLI1",
    "GLI2",
    "GLI3",
    "FOXA2",
    "ATF4",
    "JUN",
    "FOS",
    "STAT3",
    "CREB1",
    "SMAD1",
    "SMAD2",
    "SMAD3",
    "SMAD4",
    "SMAD5",
    "TCF7L2",
    "LEF1",
    "NFATC1",
    "NFKB1",
    "RELA",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tf", action="append", default=[])
    args = parser.parse_args()

    text = args.input.read_text(encoding="utf-8")
    first_motif = text.find("\nMOTIF ")
    if first_motif < 0:
        raise SystemExit("No MOTIF blocks found")
    header = text[: first_motif + 1]
    blocks = re.split(r"(?=^MOTIF )", text[first_motif + 1 :], flags=re.MULTILINE)
    wanted = {name.upper() for name in (args.tf or DEFAULT_TFS)}
    selected: list[str] = []
    selected_names: list[str] = []
    for block in blocks:
        if not block.startswith("MOTIF "):
            continue
        motif_header = block.splitlines()[0]
        fields = motif_header.split(maxsplit=2)
        alt_name = fields[2] if len(fields) > 2 else fields[1]
        tokens = {
            token.upper() for token in re.split(r"[^A-Za-z0-9]+", alt_name) if token
        }
        if tokens & wanted:
            selected.append(block.rstrip() + "\n")
            selected_names.append(motif_header)
    if not selected:
        raise SystemExit("No targeted TF motifs matched the database")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(header + "".join(selected), encoding="utf-8")
    print(f"Selected {len(selected)} motifs: " + "; ".join(selected_names))


if __name__ == "__main__":
    main()
