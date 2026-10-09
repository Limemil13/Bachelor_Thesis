"""reading every msa -> measuring seq len, gaps, comparable positions, pairwise identity and flags unusual stuff and thats tuff
dependent on annotation quality but should be fine as for now it has always been, stay positive

"""

import csv
from itertools import combinations
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ALIGNMENT_DIR = BASE / "alignments" / "canonical"
OUT_FILE = ALIGNMENT_DIR / "msa_quality_summary_canonical.tsv"
ALIGNMENT_FILES = {
    "BGN": ALIGNMENT_DIR / "BGN_canonical_aligned.faa",
    "DCN": ALIGNMENT_DIR / "DCN_canonical_aligned.faa",
    "EPYC": ALIGNMENT_DIR / "EPYC_canonical_aligned.faa",
    "FMOD": ALIGNMENT_DIR / "FMOD_canonical_aligned.faa",
    "OGN": ALIGNMENT_DIR / "OGN_canonical_aligned.faa",
    "PRELP": ALIGNMENT_DIR / "PRELP_canonical_aligned.faa",
    "LUM": ALIGNMENT_DIR / "LUM_canonical_aligned.faa",
}
EXPECTED_COUNTS = {
    "BGN": 13,
    "DCN": 15,
    "EPYC": 15,
    "FMOD": 14,
    "OGN": 14,
    "PRELP": 16,
    "LUM": 16,
}


def read_fasta(path: Path):
    # Alignment gaps are measurements in this QC step, so they must be retained.
    # Complete headers are also kept to verify gene/species/accession identity.
    records = []
    header = None
    seq_chunks = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        records.append((header, "".join(seq_chunks)))

    return records


def parse_header(header: str):
    first = header.split()[0]

    if "|" in first:
        parts = first.split("|")
        if len(parts) >= 3:
            return parts[1], parts[0], parts[2]

    if "_" in first:
        # useful fallback
        parts = first.split("_")
        return "UNKNOWN", first, ""

    return "UNKNOWN", first, ""


def pairwise_identity(seq_a: str, seq_b: str) -> tuple[float, int]:
    #comparing the residues only when both have aa, gaps are not counted as mismatch
    comparable = [
        (a, b) for a, b in zip(seq_a, seq_b, strict=False) if a != "-" and b != "-"
    ]
    if not comparable:
        return 0.0, 0
    matches = sum(a == b for a, b in comparable)
    return (matches / len(comparable)) * 100.0, len(comparable)


def main():
    # Validate panel membership and rectangular alignment shape before calculating qc
    rows = []
    seen_keys = set()

    for gene_from_file, aln_file in ALIGNMENT_FILES.items():
        if not aln_file.exists():
            raise FileNotFoundError(f"Missing canonical alignment: {aln_file}")
        records = read_fasta(aln_file)

        if len(records) != EXPECTED_COUNTS[gene_from_file]:
            raise ValueError(
                f"{aln_file} contains {len(records)} sequences; "
                f"expected {EXPECTED_COUNTS[gene_from_file]}"
            )

        aln_lengths = [len(seq) for _, seq in records]
        if len(set(aln_lengths)) != 1:
            raise ValueError(f"Alignment rows have unequal lengths in {aln_file}")
        alignment_length = aln_lengths[0]

        pair_values: dict[int, list[tuple[float, int]]] = {
            index: [] for index in range(len(records))
        }
        for left, right in combinations(range(len(records)), 2):
            identity, comparable_sites = pairwise_identity(
                records[left][1], records[right][1]
            )
            pair_values[left].append((identity, comparable_sites))
            pair_values[right].append((identity, comparable_sites))

        for index, (header, seq) in enumerate(records):
            gene, species, accession = parse_header(header)

            if gene == "UNKNOWN":
                gene = gene_from_file
            if gene != gene_from_file:
                raise ValueError(
                    f"{aln_file} contains gene {gene}, expected {gene_from_file}"
                )

            key = (gene, species, accession)
            if key in seen_keys:
                raise ValueError(f"Duplicate alignment key: {key}")
            seen_keys.add(key)

            gap_count = seq.count("-")
            ungapped_length = len(seq.replace("-", ""))
            gap_percent = round((gap_count / len(seq)) * 100, 2)
            identities = [value for value, _ in pair_values[index]]
            comparable_counts = [count for _, count in pair_values[index]]
            mean_pairwise_identity = round(sum(identities) / len(identities), 2)
            mean_comparable_sites = round(
                sum(comparable_counts) / len(comparable_counts), 1
            )
            if gap_percent > 40:
                msa_status = "check_high_gap"
            elif gap_percent > 25:
                msa_status = "watch"
            else:
                msa_status = "ok"

            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "accession": accession,
                    "alignment_file": aln_file.name,
                    "alignment_length": alignment_length,
                    "ungapped_length": ungapped_length,
                    "gap_count": gap_count,
                    "gap_percent": gap_percent,
                    "mean_pairwise_identity_percent": mean_pairwise_identity,
                    "mean_pairwise_comparable_sites": mean_comparable_sites,
                    "msa_status": msa_status,
                }
            )

    expected_total = sum(EXPECTED_COUNTS.values())
    if len(rows) != expected_total:
        raise ValueError(
            f"Expected {expected_total} canonical MSA rows, found {len(rows)}"
        )

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with OUT_FILE.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "gene",
            "species",
            "accession",
            "alignment_file",
            "alignment_length",
            "ungapped_length",
            "gap_count",
            "gap_percent",
            "mean_pairwise_identity_percent",
            "mean_pairwise_comparable_sites",
            "msa_status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT_FILE}")
    print(f"Rows: {len(rows)}")


if __name__ == "__main__":
    main()
