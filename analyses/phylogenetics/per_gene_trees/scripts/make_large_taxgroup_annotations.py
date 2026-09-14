#!/usr/bin/env python3

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

try:
    from ete3 import NCBITaxa
except ImportError:  # Existing mapping files are a reproducible offline fallback.
    NCBITaxa = None


GROUPS = {
    "Mammal": "#d73027",
    "Bird": "#fee08b",
    "Reptile": "#91cf60",
    "Amphibian": "#66c2a5",
    "Ray-finned fish": "#4575b4",
    "Lobe-finned fish": "#3288bd",
    "Cartilaginous fish": "#984ea3",
    "Jawless vertebrate": "#a6d854",
    "Other chordate": "#bdbdbd",
    "Other vertebrate": "#999999",
    "Non-vertebrate": "#cccccc",
    "Unknown": "#000000",
}


# Species absent from the five already annotated large-gene trees. These broad
# vertebrate assignments are stable and make the LUM annotation reproducible
# even on a machine without ete3's local NCBI taxonomy database.
MANUAL_SPECIES_GROUPS = {
    "Amblyraja radiata": "Cartilaginous fish",
    "Aulonocara sp. 'chitande type north' Nkhata Bay": "Ray-finned fish",
    "Heptranchias perlo": "Cartilaginous fish",
    "Lampetra planeri": "Jawless vertebrate",
    "Sinocyclocheilus rhinocerous": "Ray-finned fish",
}


def clean_tip_name(name: str) -> str:
    return name.strip().strip("'").strip('"')


def extract_tip_names(treefile: Path):
    text = treefile.read_text().strip()
    raw_names = re.findall(r"(?<=[(,])([^(),:;]+):", text)
    return sorted(set(clean_tip_name(n) for n in raw_names))


def extract_organism(old_header: str) -> str:
    m = re.search(r"\[organism=([^\]]+)\]", old_header)
    if m:
        return m.group(1).strip()
    return "Unknown"


def tree_tip_from_report_header(new_header: str) -> str:
    # filter report has: GeneID_125260808|BGN|XP_048035309.1
    # iTOL tree has after sed: GeneID_125260808_BGN_XP_048035309.1
    return new_header.replace("|", "_")


def load_species_mapping(report_file: Path):
    mapping = {}

    with report_file.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if row["status"] != "kept":
                continue

            tip = tree_tip_from_report_header(row["new_header"])
            organism = extract_organism(row["old_header"])

            mapping[tip] = organism

    return mapping


def load_cached_taxonomy(mapping_files):
    """Load species-level classifications from prior final tree mappings."""
    cache = {}
    for path in mapping_files:
        if not path.exists():
            continue
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                species = row.get("species", "Unknown")
                group = row.get("group", "Unknown")
                lineage = row.get("lineage", "NA")
                if species != "Unknown" and group in GROUPS and group != "Unknown":
                    cache.setdefault(species, (group, lineage))
    return cache


def classify_lineage(lineage_names):
    names = set(lineage_names)

    if "Mammalia" in names:
        return "Mammal"
    if "Aves" in names:
        return "Bird"
    if (
        "Lepidosauria" in names
        or "Testudines" in names
        or "Crocodylia" in names
        or "Reptilia" in names
    ):
        return "Reptile"
    if "Amphibia" in names:
        return "Amphibian"

    if "Actinopteri" in names or "Actinopterygii" in names:
        return "Ray-finned fish"
    if "Chondrichthyes" in names:
        return "Cartilaginous fish"
    if "Petromyzontida" in names or "Myxini" in names or "Cyclostomata" in names:
        return "Jawless vertebrate"

    # Important: check Sarcopterygii only after mammals/birds/reptiles/amphibians,
    # because tetrapods are nested inside Sarcopterygii.
    if "Sarcopterygii" in names:
        return "Lobe-finned fish"

    if "Vertebrata" in names:
        return "Other vertebrate"
    if "Chordata" in names or "Cephalochordata" in names or "Tunicata" in names:
        return "Other chordate"

    return "Non-vertebrate"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--filter-report", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument(
        "--reuse-mapping",
        action="append",
        default=[],
        help="Existing taxonomy_mapping.tsv file to reuse (repeatable).",
    )
    args = parser.parse_args()

    treefile = Path(args.tree)
    report_file = Path(args.filter_report)
    out_prefix = Path(args.out_prefix)

    tips = extract_tip_names(treefile)
    tip_to_species = load_species_mapping(report_file)

    unique_species = sorted(set(tip_to_species.values()) - {"Unknown"})

    print(f"Tips in tree: {len(tips)}")
    print(f"Tips with species mapping: {sum(1 for t in tips if t in tip_to_species)}")
    print(f"Unique species names: {len(unique_species)}")
    print("Example species:")
    for s in unique_species[:10]:
        print(" ", s)

    mapping_files = [Path(path) for path in args.reuse_mapping]
    if not mapping_files:
        # Files are emitted via Path.with_suffix(".taxonomy_mapping.tsv"), so
        # the separator before "taxonomy" is a dot rather than an underscore.
        mapping_files = sorted(out_prefix.parent.glob("*.taxonomy_mapping.tsv"))
    cached_taxonomy = load_cached_taxonomy(mapping_files)

    ncbi = NCBITaxa() if NCBITaxa is not None else None

    group_by_species = {}
    lineage_by_species = {}

    for species in unique_species:
        if species in cached_taxonomy:
            group_by_species[species], lineage_by_species[species] = cached_taxonomy[
                species
            ]
            continue
        if species in MANUAL_SPECIES_GROUPS:
            group_by_species[species] = MANUAL_SPECIES_GROUPS[species]
            lineage_by_species[species] = "manual broad-group assignment"
            continue
        if ncbi is None:
            group_by_species[species] = "Unknown"
            lineage_by_species[species] = "NA"
            continue
        try:
            translated = ncbi.get_name_translator([species])
            taxids = translated.get(species, [])

            if not taxids:
                group_by_species[species] = "Unknown"
                lineage_by_species[species] = "NA"
                continue

            taxid = taxids[0]
            lineage = ncbi.get_lineage(taxid)
            lineage_names_dict = ncbi.get_taxid_translator(lineage)
            lineage_names = [
                lineage_names_dict[t] for t in lineage if t in lineage_names_dict
            ]

            group = classify_lineage(lineage_names)

            group_by_species[species] = group
            lineage_by_species[species] = ";".join(lineage_names)

        except Exception:
            group_by_species[species] = "Unknown"
            lineage_by_species[species] = "NA"

    # iTOL color strip
    colorstrip_file = out_prefix.with_suffix(".taxgroups.txt")

    legend_labels = list(GROUPS.keys())
    legend_colors = [GROUPS[g] for g in legend_labels]
    legend_shapes = ["1"] * len(legend_labels)

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tTaxonomic_group",
        "COLOR\t#000000",
        "STRIP_WIDTH\t35",
        "MARGIN\t10",
        "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tTaxonomic group",
        "LEGEND_SHAPES\t" + "\t".join(legend_shapes),
        "LEGEND_COLORS\t" + "\t".join(legend_colors),
        "LEGEND_LABELS\t" + "\t".join(legend_labels),
        "DATA",
    ]

    group_counts = Counter()

    for tip in tips:
        species = tip_to_species.get(tip, "Unknown")
        group = group_by_species.get(species, "Unknown")
        color = GROUPS.get(group, GROUPS["Unknown"])

        group_counts[group] += 1
        lines.append(f"{tip}\t{color}\t{group}")

    colorstrip_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # report
    mapping_file = out_prefix.with_suffix(".taxonomy_mapping.tsv")

    report_lines = ["tip\tspecies\tgroup\tlineage"]
    for tip in tips:
        species = tip_to_species.get(tip, "Unknown")
        group = group_by_species.get(species, "Unknown")
        lineage = lineage_by_species.get(species, "NA")
        report_lines.append(f"{tip}\t{species}\t{group}\t{lineage}")

    mapping_file.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Wrote: {colorstrip_file}")
    print(f"Wrote: {mapping_file}")
    print(f"Cached mapping files reused: {len(mapping_files)}")
    print("Group counts:")
    for group, count in group_counts.most_common():
        print(f"  {group}: {count}")


if __name__ == "__main__":
    main()
