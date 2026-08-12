from pathlib import Path
from Bio import Phylo
import sys

if len(sys.argv) != 3:
    print("Usage: python make_itol_annotations.py <treefile> <output_prefix>")
    sys.exit(1)

treefile = Path(sys.argv[1])
out_prefix = Path(sys.argv[2])

tree = Phylo.read(treefile, "newick")

# ---------- species naming ----------
label_map = {
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

    "latimeria_chalumnae": ("Coelacanth", "Lobe-finned fish", "#3288bd"),

    "lepisosteus_oculatus": ("Spotted gar", "Ray-finned fish", "#4575b4"),
    "zebrafish": ("Zebrafish", "Ray-finned fish", "#4575b4"),
    "danio_rerio": ("Zebrafish", "Ray-finned fish", "#4575b4"),

    "callorhinchus_milii": ("Elephant shark", "Cartilaginous fish", "#984ea3"),
    "scyliorhinus_canicula": ("Catshark", "Cartilaginous fish", "#984ea3"),
    "whale_shark": ("Whale shark", "Cartilaginous fish", "#984ea3"),
    "rhincodon_typus": ("Whale shark", "Cartilaginous fish", "#984ea3"),

    "lamprey": ("Lamprey", "Jawless vertebrate", "#a6d854"),
    "petromyzon_marinus": ("Lamprey", "Jawless vertebrate", "#a6d854"),

    "amphioxus": ("Amphioxus", "Outgroup chordate", "#bdbdbd"),
    "branchiostoma": ("Amphioxus", "Outgroup chordate", "#bdbdbd"),
}

def detect_species(raw_name: str):
    s = raw_name.lower()
    for key, value in label_map.items():
        if key in s:
            return value
    return (raw_name, "Unknown", "#999999")

# ---------- output files ----------
labels_file = out_prefix.with_suffix(".labels.txt")
strip_file = out_prefix.with_suffix(".groups.txt")

# iTOL LABELS dataset
label_lines = [
    "LABELS",
    "SEPARATOR TAB",
    "DATA"
]

# iTOL COLORSTRIP dataset
strip_lines = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tEvolutionary_group",
    "COLOR\t#000000",
    "STRIP_WIDTH\t25",
    "MARGIN\t5",
    "SHOW_INTERNAL\t0",
    "DATA"
]

for clade in tree.get_terminals():
    old = clade.name
    short_label, group, color = detect_species(old)
    label_lines.append(f"{old}\t{short_label}")
    strip_lines.append(f"{old}\t{color}\t{group}")

labels_file.write_text("\n".join(label_lines) + "\n", encoding="utf-8")
strip_file.write_text("\n".join(strip_lines) + "\n", encoding="utf-8")

print(f"Wrote: {labels_file}")
print(f"Wrote: {strip_file}")