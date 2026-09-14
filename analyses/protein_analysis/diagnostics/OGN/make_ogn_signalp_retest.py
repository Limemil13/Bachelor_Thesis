from pathlib import Path

inp = Path("analyses/protein_analysis/candidates/OGN_domain_input.faa")
out = Path("analyses/protein_analysis/diagnostics/OGN/OGN_signalp_retest.faa")

wanted = {
    "human": "reference_normal",
    "mouse": "reference_normal",
    "bos_taurus": "reference_normal",
    "canis_lupus_familiaris": "dog_full",
    "lepisosteus_oculatus": "gar_full",
}

records = []
header = None
seq = []

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

with out.open("w") as f:
    for h, s in records:
        first = h.split()[0]
        parts = first.split("|")
        species = parts[0]
        accession = parts[2]

        if species not in wanted:
            continue

        simple = f"OGN_{species}_{accession}_{wanted[species]}"
        f.write(f">{simple}\n{s}\n")

        # special dog trimmed test
        if species == "canis_lupus_familiaris":
            motif = "MKTLQST"
            pos = s.find(motif)
            print("Dog MKTLQST position, 0-based:", pos)
            if pos != -1:
                trimmed = s[pos:]
                f.write(
                    f">OGN_canis_lupus_familiaris_{accession}_dog_trimmed_from_MKTLQST\n{trimmed}\n"
                )
                print("Dog full length:", len(s))
                print("Dog trimmed length:", len(trimmed))

print("Wrote:", out)
