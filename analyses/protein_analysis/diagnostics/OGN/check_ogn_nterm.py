import csv
from pathlib import Path

fasta = Path("analyses/protein_analysis/candidates/OGN_domain_input.faa")
signalp = Path(
    "analyses/protein_analysis/signal_peptides/tables/signalp_summary_clean.tsv"
)
out = Path("analyses/protein_analysis/diagnostics/OGN/OGN_nterm_check.tsv")

# Read SignalP table
signalp_info = {}
with signalp.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["gene"] == "OGN":
            key = row["signalp_id"]
            signalp_info[key] = row

# Read FASTA
records = []
header = None
seq = []

with fasta.open() as f:
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

rows = []

for header, sequence in records:
    first = header.split()[0]
    parts = first.split("|")

    if len(parts) >= 3:
        species = parts[0]
        gene = parts[1]
        accession = parts[2]
    else:
        species = "unknown"
        gene = "OGN"
        accession = first

    signalp_id = f"OGN_{species}_{accession}"
    sp = signalp_info.get(signalp_id, {})

    rows.append(
        {
            "species": species,
            "accession": accession,
            "length": len(sequence),
            "signal_peptide": sp.get("signal_peptide", ""),
            "signalp_quality": sp.get("signalp_quality", ""),
            "cleavage_site": sp.get("cleavage_site", ""),
            "cleavage_probability": sp.get("cleavage_probability", ""),
            "nterm_80aa": sequence[:80],
        }
    )

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "species",
            "accession",
            "length",
            "signal_peptide",
            "signalp_quality",
            "cleavage_site",
            "cleavage_probability",
            "nterm_80aa",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print("Wrote:", out)
