from pathlib import Path

inp = Path("analyses/protein_analysis/candidates/OGN_domain_input.faa")
out = Path("analyses/protein_analysis/diagnostics/OGN/OGN_nterm_100aa.faa")

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
        simple = h.split()[0].replace("|", "_")
        f.write(f">{simple}\n{s[:100]}\n")

print("Wrote:", out)
