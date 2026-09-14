#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def extract_tip_names(treefile: Path):
    text = treefile.read_text().strip()
    # Finds Newick tip labels: after "(" or "," and before ":"
    names = re.findall(r"(?<=[(,])([^(),:;]+):", text)
    return sorted(set(names))


def make_label(name: str, gene: str):
    delimiter = f"_{gene}_"

    if delimiter in name:
        species_part, accession = name.split(delimiter, 1)
        species_label = species_part.replace("_", " ")
        return species_label

    return name.replace("_", " ")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    treefile = Path(args.tree)
    out_prefix = Path(args.out_prefix)

    tips = extract_tip_names(treefile)

    labels_file = out_prefix.with_suffix(".labels.txt")

    lines = [
        "LABELS",
        "SEPARATOR TAB",
        "DATA",
    ]

    for tip in tips:
        label = make_label(tip, args.gene)
        lines.append(f"{tip}\t{label}")

    labels_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Tips found: {len(tips)}")
    print(f"Wrote: {labels_file}")


if __name__ == "__main__":
    main()
