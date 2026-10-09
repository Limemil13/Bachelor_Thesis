"""Write a diagnostic OGN FASTA with the dog start-site correction."""

from pathlib import Path

inp = Path("analyses/protein_analysis/candidates/OGN_domain_input.faa")
out = Path(
    "analyses/protein_analysis/diagnostics/OGN/OGN_candidates_dog_start_corrected.faa"
)
report = Path("analyses/protein_analysis/diagnostics/OGN/OGN_dog_correction_report.tsv")

records = []
header = None
seq = []

# Parse all original OGN records; the source FASTA remains unchanged.
with inp.open() as f:
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

# Apply the documented dog start-site correction and report the exact trimming.
with out.open("w") as fasta_out, report.open("w") as rep:
    rep.write("species\taccession\toriginal_length\tcorrected_length\tcorrection\n")

    for h, s in records:
        first = h.split()[0]
        parts = first.split("|")

        if len(parts) >= 3:
            species = parts[0]
            accession = parts[2]
        else:
            species = "unknown"
            accession = first

        correction = "none"
        new_seq = s
        new_header = h

        if species == "canis_lupus_familiaris" and accession == "XP_038383338.1":
            motif = "MKTLQST"
            pos = s.find(motif)
            if pos != -1:
                new_seq = s[pos:]
                new_header = h + " dog_start_corrected_from_MKTLQST"
                correction = f"trimmed_first_{pos}_aa"

        fasta_out.write(f">{new_header}\n")
        for i in range(0, len(new_seq), 60):
            fasta_out.write(new_seq[i : i + 60] + "\n")

        rep.write(f"{species}\t{accession}\t{len(s)}\t{len(new_seq)}\t{correction}\n")

print("Wrote:", out)
print("Wrote:", report)
