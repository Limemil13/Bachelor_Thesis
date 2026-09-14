#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

CLASS_COLORS = {
    "Class_I": "#d73027",
    "Class_II": "#4575b4",
    "Class_III": "#984ea3",
    "Class_IV": "#66c2a5",
    "Class_V": "#fdae61",
    "Outgroup": "#999999",
    "Unknown": "#000000",
}

THESIS_GENES = {"BGN", "FMOD", "OGN", "EPYC", "PRELP", "LUM"}


def extract_tip_names(treefile: Path):
    text = treefile.read_text().strip()
    tips = re.findall(r"(?<=[(,])([^(),:;]+):", text)
    return sorted(set(t.strip().strip("'").strip('"') for t in tips))


def parse_tip(tip):
    # After sed, header is GENE_CLASS_ACC, e.g. BGN_Class_I_NP_001702.1
    parts = tip.split("_")

    gene = parts[0]

    if "Class" in tip:
        if "Class_I_" in tip:
            cls = "Class_I"
        elif "Class_II_" in tip:
            cls = "Class_II"
        elif "Class_III_" in tip:
            cls = "Class_III"
        elif "Class_IV_" in tip:
            cls = "Class_IV"
        elif "Class_V_" in tip:
            cls = "Class_V"
        else:
            cls = "Unknown"
    elif "Outgroup" in tip:
        cls = "Outgroup"
    else:
        cls = "Unknown"

    return gene, cls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    treefile = Path(args.tree)
    out_prefix = Path(args.out_prefix)

    tips = extract_tip_names(treefile)

    labels_file = out_prefix.with_suffix(".labels.txt")
    classes_file = out_prefix.with_suffix(".classes.txt")
    thesis_file = out_prefix.with_suffix(".thesis_genes.txt")

    label_lines = [
        "LABELS",
        "SEPARATOR TAB",
        "DATA",
    ]

    class_lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tSLRP_class",
        "COLOR\t#000000",
        "STRIP_WIDTH\t35",
        "MARGIN\t5",
        "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tSLRP class",
        "LEGEND_SHAPES\t1\t1\t1\t1\t1\t1",
        "LEGEND_COLORS\t#d73027\t#4575b4\t#984ea3\t#66c2a5\t#fdae61\t#999999",
        "LEGEND_LABELS\tClass I\tClass II\tClass III\tClass IV\tClass V\tOutgroup",
        "DATA",
    ]

    thesis_lines = [
        "DATASET_BINARY",
        "SEPARATOR TAB",
        "DATASET_LABEL\tThesis_gene",
        "COLOR\t#000000",
        "FIELD_SHAPES\t1",
        "FIELD_LABELS\tThesis genes",
        "FIELD_COLORS\t#000000",
        "DATA",
    ]

    for tip in tips:
        gene, cls = parse_tip(tip)
        color = CLASS_COLORS.get(cls, CLASS_COLORS["Unknown"])

        label_lines.append(f"{tip}\t{gene}")
        class_lines.append(f"{tip}\t{color}\t{cls}")

        value = "1" if gene in THESIS_GENES else "0"
        thesis_lines.append(f"{tip}\t{value}")

    labels_file.write_text("\n".join(label_lines) + "\n", encoding="utf-8")
    classes_file.write_text("\n".join(class_lines) + "\n", encoding="utf-8")
    thesis_file.write_text("\n".join(thesis_lines) + "\n", encoding="utf-8")

    print(f"Tips found: {len(tips)}")
    print(f"Wrote: {labels_file}")
    print(f"Wrote: {classes_file}")
    print(f"Wrote: {thesis_file}")


if __name__ == "__main__":
    main()
