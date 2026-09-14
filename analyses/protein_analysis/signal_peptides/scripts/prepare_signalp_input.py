import csv
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
IN_DIR = BASE / "candidates"
OUT_DIR = BASE / "signal_peptides" / "input"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FASTA = OUT_DIR / "slrp_candidates_canonical_for_signalp.faa"
OUT_MAP = OUT_DIR / "signalp_canonical_id_mapping.tsv"
OUT_LUM_FASTA = OUT_DIR / "LUM_canonical_for_signalp.faa"
OUT_LUM_MAP = OUT_DIR / "LUM_signalp_id_mapping.tsv"
OUT_DCN_FASTA = OUT_DIR / "DCN_canonical_for_signalp.faa"
OUT_DCN_MAP = OUT_DIR / "DCN_signalp_id_mapping.tsv"
OUT_PENDING_FASTA = OUT_DIR / "pending_17_canonical_for_signalp.faa"
OUT_PENDING_MAP = OUT_DIR / "pending_17_signalp_id_mapping.tsv"

files = {
    "BGN": IN_DIR / "canonical_by_gene" / "BGN_canonical.faa",
    "DCN": IN_DIR / "canonical_by_gene" / "DCN_canonical.faa",
    "FMOD": IN_DIR / "canonical_by_gene" / "FMOD_canonical.faa",
    "OGN": IN_DIR / "canonical_by_gene" / "OGN_canonical.faa",
    "EPYC": IN_DIR / "canonical_by_gene" / "EPYC_canonical.faa",
    "PRELP": IN_DIR / "canonical_by_gene" / "PRELP_canonical.faa",
    "LUM": IN_DIR / "canonical_by_gene" / "LUM_canonical.faa",
}


def read_fasta(path):
    records = []
    header = None
    seq = []
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


def clean_text(x):
    x = x.replace("|", "_")
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


rows = []
seen_ids = set()
expected_counts = {
    "BGN": 13,
    "DCN": 15,
    "FMOD": 14,
    "OGN": 14,
    "EPYC": 15,
    "PRELP": 16,
    "LUM": 16,
}
observed_counts = {}

with OUT_FASTA.open("w") as fasta_out:
    for gene, path in files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing canonical SignalP source: {path}")

        records = read_fasta(path)
        observed_counts[gene] = len(records)
        for header, seq in records:
            first = header.split()[0]
            parts = first.split("|")

            if len(parts) >= 3:
                species = parts[0]
                accession = parts[2]
            else:
                species = "unknown"
                accession = first

            simple_id = clean_text(f"{gene}__{species}__{accession}")

            # Avoid duplicate IDs if any exist
            original_simple_id = simple_id
            counter = 2
            while simple_id in seen_ids:
                simple_id = f"{original_simple_id}_{counter}"
                counter += 1
            seen_ids.add(simple_id)

            fasta_out.write(f">{simple_id}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

            rows.append(
                {
                    "signalp_id": simple_id,
                    "gene": gene,
                    "species": species,
                    "accession": accession,
                    "protein_length": len(seq),
                    "original_header": header,
                    "source_file": path.relative_to(BASE.parents[1]).as_posix(),
                }
            )

if observed_counts != expected_counts:
    raise ValueError(f"Unexpected per-gene SignalP counts: {observed_counts}")
expected_total = sum(expected_counts.values())
if len(rows) != expected_total:
    raise ValueError(
        f"Expected {expected_total} canonical SignalP proteins, found {len(rows)}"
    )

with OUT_MAP.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "signalp_id",
            "gene",
            "species",
            "accession",
            "protein_length",
            "original_header",
            "source_file",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

lum_rows = [row for row in rows if row["gene"] == "LUM"]
lum_ids = {row["signalp_id"] for row in lum_rows}
with OUT_LUM_FASTA.open("w") as fasta_out:
    for header, seq in read_fasta(OUT_FASTA):
        if header.split()[0] in lum_ids:
            fasta_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

with OUT_LUM_MAP.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "signalp_id",
            "gene",
            "species",
            "accession",
            "protein_length",
            "original_header",
            "source_file",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(lum_rows)

if len(lum_rows) != 16:
    raise ValueError(f"Expected 16 LUM SignalP proteins, found {len(lum_rows)}")

dcn_rows = [row for row in rows if row["gene"] == "DCN"]
dcn_ids = {row["signalp_id"] for row in dcn_rows}
with OUT_DCN_FASTA.open("w") as fasta_out:
    for header, seq in read_fasta(OUT_FASTA):
        if header.split()[0] in dcn_ids:
            fasta_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

with OUT_DCN_MAP.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "signalp_id",
            "gene",
            "species",
            "accession",
            "protein_length",
            "original_header",
            "source_file",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(dcn_rows)

if len(dcn_rows) != 15:
    raise ValueError(f"Expected 15 DCN SignalP proteins, found {len(dcn_rows)}")

pending_bgn_accessions = {"XP_038638277.1", "XP_048475905.1"}
pending_rows = [
    row
    for row in rows
    if row["gene"] == "DCN"
    or (row["gene"] == "BGN" and row["accession"] in pending_bgn_accessions)
]
pending_ids = {row["signalp_id"] for row in pending_rows}
with OUT_PENDING_FASTA.open("w") as fasta_out:
    for header, seq in read_fasta(OUT_FASTA):
        if header.split()[0] in pending_ids:
            fasta_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

with OUT_PENDING_MAP.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "signalp_id",
            "gene",
            "species",
            "accession",
            "protein_length",
            "original_header",
            "source_file",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(pending_rows)

if len(pending_rows) != 17:
    raise ValueError(f"Expected 17 pending SignalP proteins, found {len(pending_rows)}")

print("Wrote:")
print(OUT_FASTA)
print(OUT_MAP)
print(OUT_LUM_FASTA)
print(OUT_LUM_MAP)
print(OUT_DCN_FASTA)
print(OUT_DCN_MAP)
print(OUT_PENDING_FASTA)
print(OUT_PENDING_MAP)
print("Sequences:", len(rows))
