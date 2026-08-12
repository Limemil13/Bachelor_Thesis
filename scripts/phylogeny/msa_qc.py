from pathlib import Path
import csv


ALIGNMENT_DIR = Path("analyses/slrp_phylogeny/04_alignments")
OUT_FILE = Path("analyses/slrp_phylogeny/04_alignments/msa_quality_summary.tsv")


def read_fasta(path: Path):
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
    """
    Expected header example:
    human|BGN|XP_016885213.1 forward_pident=100.00 ...
    or cleaned:
    human_BGN_XP_016885213.1
    """
    first = header.split()[0]

    if "|" in first:
        parts = first.split("|")
        if len(parts) >= 3:
            return parts[1], parts[0], parts[2]

    if "_" in first:
        # fallback, not perfect but useful
        parts = first.split("_")
        return "UNKNOWN", first, ""

    return "UNKNOWN", first, ""


def main():
    rows = []

    for aln_file in sorted(ALIGNMENT_DIR.glob("*_aligned.faa")):
        gene_from_file = aln_file.name.replace("_aligned.faa", "")
        records = read_fasta(aln_file)

        if not records:
            continue

        aln_lengths = [len(seq) for _, seq in records]
        alignment_length = max(aln_lengths)

        for header, seq in records:
            gene, target, hit_id = parse_header(header)

            if gene == "UNKNOWN":
                gene = gene_from_file

            gap_count = seq.count("-")
            ungapped_length = len(seq.replace("-", ""))
            gap_percent = round((gap_count / len(seq)) * 100, 2)

            if gap_percent > 40:
                msa_status = "check_high_gap"
            elif gap_percent > 25:
                msa_status = "watch"
            else:
                msa_status = "ok"

            rows.append({
                "gene": gene,
                "target": target,
                "hit_id": hit_id,
                "alignment_file": aln_file.name,
                "alignment_length": alignment_length,
                "ungapped_length": ungapped_length,
                "gap_count": gap_count,
                "gap_percent": gap_percent,
                "msa_status": msa_status,
            })

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with OUT_FILE.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "gene",
            "target",
            "hit_id",
            "alignment_file",
            "alignment_length",
            "ungapped_length",
            "gap_count",
            "gap_percent",
            "msa_status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT_FILE}")
    print(f"Rows: {len(rows)}")


if __name__ == "__main__":
    main()