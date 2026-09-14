#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

GENE_COLORS = {
    "BGN": "#d73027",
    "DCN": "#8c564b",
    "FMOD": "#4575b4",
    "OGN": "#66c2a5",
    "EPYC": "#984ea3",
    "PRELP": "#fdae61",
    "LUM": "#e7298a",
    "Unknown": "#999999",
}


TAXON_MAP = {
    "human": ("Human", "Placental mammal", "#d73027"),
    "mouse": ("Mouse", "Placental mammal", "#d73027"),
    "rattus_norvegicus": ("Rat", "Placental mammal", "#d73027"),
    "canis_lupus_familiaris": ("Dog", "Placental mammal", "#d73027"),
    "bos_taurus": ("Cow", "Placental mammal", "#d73027"),
    "monodelphis_domestica": ("Opossum", "Marsupial", "#fc8d59"),
    "chicken": ("Chicken", "Bird", "#fee08b"),
    "gallus_gallus": ("Chicken", "Bird", "#fee08b"),
    "anolis_carolinensis": ("Anole", "Reptile", "#91cf60"),
    "frog": ("Frog", "Amphibian", "#66c2a5"),
    "xenopus": ("Frog", "Amphibian", "#66c2a5"),
    "zebrafish": ("Zebrafish", "Ray-finned fish", "#4575b4"),
    "danio_rerio": ("Zebrafish", "Ray-finned fish", "#4575b4"),
    "lepisosteus_oculatus": ("Spotted gar", "Ray-finned fish", "#4575b4"),
    "latimeria_chalumnae": ("Coelacanth", "Lobe-finned fish", "#3288bd"),
    "callorhinchus_milii": ("Elephant shark", "Cartilaginous fish", "#984ea3"),
    "scyliorhinus_canicula": ("Catshark", "Cartilaginous fish", "#984ea3"),
    "whale_shark": ("Whale shark", "Cartilaginous fish", "#984ea3"),
    "rhincodon_typus": ("Whale shark", "Cartilaginous fish", "#984ea3"),
    "lamprey": ("Lamprey", "Jawless vertebrate", "#a6d854"),
    "petromyzon_marinus": ("Lamprey", "Jawless vertebrate", "#a6d854"),
    "amphioxus": ("Amphioxus", "Outgroup chordate", "#bdbdbd"),
    "branchiostoma": ("Amphioxus", "Outgroup chordate", "#bdbdbd"),
}


def extract_tip_names(treefile: Path):
    text = treefile.read_text().strip()
    tips = re.findall(r"(?<=[(,])([^(),:;]+):", text)
    return sorted(set(t.strip().strip("'").strip('"') for t in tips))


def detect_gene(tip: str):
    parts = tip.split("_")
    for gene in GENE_COLORS:
        if gene != "Unknown" and gene in parts:
            return gene
    return "Unknown"


def detect_species(tip: str):
    lower = tip.lower()
    for key, value in TAXON_MAP.items():
        if key in lower:
            return value
    return (tip, "Unknown", "#999999")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    treefile = Path(args.tree)
    out_prefix = Path(args.out_prefix)

    tips = extract_tip_names(treefile)

    labels_file = out_prefix.with_suffix(".labels.txt")
    gene_file = out_prefix.with_suffix(".genes.txt")
    tax_file = out_prefix.with_suffix(".taxgroups.txt")

    # Labels
    label_lines = [
        "LABELS",
        "SEPARATOR TAB",
        "DATA",
    ]

    # Gene color strip
    gene_lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tGene",
        "COLOR\t#000000",
        "STRIP_WIDTH\t25",
        "MARGIN\t5",
        "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tGene",
        "LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1",
        "LEGEND_COLORS\t#d73027\t#8c564b\t#4575b4\t#66c2a5\t#984ea3\t#fdae61\t#e7298a",
        "LEGEND_LABELS\tBGN\tDCN\tFMOD\tOGN\tEPYC\tPRELP\tLUM",
        "DATA",
    ]

    # Taxonomic group color strip
    tax_lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tTaxonomic_group",
        "COLOR\t#000000",
        "STRIP_WIDTH\t25",
        "MARGIN\t5",
        "SHOW_INTERNAL\t0",
        "DATA",
    ]

    for tip in tips:
        gene = detect_gene(tip)
        gene_color = GENE_COLORS.get(gene, GENE_COLORS["Unknown"])

        species_label, tax_group, tax_color = detect_species(tip)

        # label format: Species_GENE
        clean_label = f"{species_label}_{gene}"

        label_lines.append(f"{tip}\t{clean_label}")
        gene_lines.append(f"{tip}\t{gene_color}\t{gene}")
        tax_lines.append(f"{tip}\t{tax_color}\t{tax_group}")

    labels_file.write_text("\n".join(label_lines) + "\n", encoding="utf-8")
    gene_file.write_text("\n".join(gene_lines) + "\n", encoding="utf-8")
    tax_file.write_text("\n".join(tax_lines) + "\n", encoding="utf-8")

    print(f"Tips found: {len(tips)}")
    print(f"Wrote: {labels_file}")
    print(f"Wrote: {gene_file}")
    print(f"Wrote: {tax_file}")


if __name__ == "__main__":
    main()
