#!/usr/bin/env python3

"""Assemble the reduced human SLRP reference/context FASTA."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]

OUTDIR = ROOT / "analyses/phylogenetics/combined_trees/slrp_family_reference_tree"
DOWNLOADS = OUTDIR / "downloads"
LOCAL_QUERIES = ROOT / "data/queries"

OUT_FASTA = OUTDIR / "SLRP_reference_context.faa"
OUT_REPORT = OUTDIR / "SLRP_reference_context_selected.tsv"

# Reduced SLRP family context: focal genes, close paralogs and one LINGO outgroup.
GENES = {
    # Class I / Clade 2 context
    "ASPN": {"class": "Class_I", "preferred": ["NP_060150.4", "NP_001180264.1"]},
    "BGN": {"class": "Class_I", "preferred": ["NP_001702.1"]},
    "DCN": {"class": "Class_I", "preferred": ["NP_001911.1"]},
    "ECM2": {"class": "Class_I", "preferred": []},
    # Class II / Clade 3 context
    "FMOD": {"class": "Class_II", "preferred": ["NP_002014.2"]},
    "LUM": {"class": "Class_II", "preferred": []},
    "KERA": {"class": "Class_II", "preferred": []},
    "OMD": {"class": "Class_II", "preferred": ["NP_005005.1"]},
    "PRELP": {"class": "Class_II", "preferred": ["NP_002716.1", "NP_958505.1"]},
    # Class III / Clade 4 context
    "OGN": {
        "class": "Class_III",
        "preferred": ["NP_148935.1", "NP_054776.1", "NP_077727.3"],
    },
    "EPYC": {"class": "Class_III", "preferred": ["NP_004941.2"]},
    "OPTC": {"class": "Class_III", "preferred": []},
    # Non-canonical / broader context
    "CHAD": {"class": "Class_IV", "preferred": []},
    "CHADL": {"class": "Class_IV", "preferred": []},
    "NYX": {"class": "Class_IV", "preferred": []},
    "TSKU": {"class": "Class_IV", "preferred": []},
    "PODN": {"class": "Class_V", "preferred": []},
    "PODNL1": {"class": "Class_V", "preferred": []},
    # Outgroup, similar to paper logic using Lingo family
    "LINGO1": {"class": "Outgroup", "preferred": []},
}


def read_fasta(path):
    # Preserve full headers so source accessions and descriptions remain auditable.
    records = []
    header = None
    seq = []

    if not path.exists():
        return records

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

    if header is not None:
        records.append((header, "".join(seq)))

    return records


def write_fasta(records, path):
    # Write the selected human reference sequences under normalized tree-tip names.
    with path.open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


def accession(header):
    # Recover the accession from common NCBI header formats.
    return header.split()[0]


def is_plausible_slrp_length(seq, gene):
    # Apply a broad length screen while allowing the longer outgroup proteins.
    # LINGO is larger than SLRPs, so allow it separately.
    if gene == "LINGO1":
        return 400 <= len(seq) <= 900
    return 150 <= len(seq) <= 700


def candidate_files_for_gene(gene):
    # Search pinned downloads first and local query FASTAs as a fallback.
    files = []

    local = LOCAL_QUERIES / f"{gene}.faa"
    if local.exists():
        files.append(local)

    dl = DOWNLOADS / gene / f"{gene}_human" / "ncbi_dataset" / "data" / "protein.faa"
    if dl.exists():
        files.append(dl)

    return files


def choose_record(gene, records, preferred):
    # Prefer listed reference accessions, otherwise choose a plausible complete model.
    if not records:
        return None

    # Keep only plausible length first.
    plausible = [(h, s) for h, s in records if is_plausible_slrp_length(s, gene)]
    if not plausible:
        plausible = records

    # Preferred accession first.
    for pref in preferred:
        for h, s in plausible:
            if accession(h) == pref:
                return h, s, "preferred_accession"

    # Then curated NP_ records.
    np_records = [(h, s) for h, s in plausible if accession(h).startswith("NP_")]
    if np_records:
        # Use longest NP_ as simple representative if no preferred accession specified.
        h, s = sorted(np_records, key=lambda x: len(x[1]), reverse=True)[0]
        return h, s, "longest_NP"

    # Otherwise longest plausible.
    h, s = sorted(plausible, key=lambda x: len(x[1]), reverse=True)[0]
    return h, s, "longest_plausible"


def main():
    # Select one human sequence per reference gene and record exactly which source
    # file and accession contributed to the combined context FASTA.
    selected = []
    report = [
        "gene\tclass\taccession\tlength\tselection_reason\told_header\tnew_header"
    ]

    for gene, meta in GENES.items():
        all_records = []

        for f in candidate_files_for_gene(gene):
            all_records.extend(read_fasta(f))

        chosen = choose_record(gene, all_records, meta["preferred"])

        if chosen is None:
            print(f"WARNING: no sequence found for {gene}")
            continue

        old_header, seq, reason = chosen
        acc = accession(old_header)
        cls = meta["class"]
        new_header = f"{gene}|{cls}|{acc}"

        selected.append((new_header, seq))
        report.append(
            f"{gene}\t{cls}\t{acc}\t{len(seq)}\t{reason}\t{old_header}\t{new_header}"
        )

    write_fasta(selected, OUT_FASTA)
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Selected proteins: {len(selected)}")
    print(f"Wrote: {OUT_FASTA}")
    print(f"Wrote: {OUT_REPORT}")


if __name__ == "__main__":
    main()
