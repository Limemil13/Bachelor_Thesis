import sys
import pandas as pd
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python scripts/extract_gene_tpm.py <quant.sf> <output.tsv>")
    sys.exit(1)

genes_file = "genes_of_interest.tsv"
quant_file = sys.argv[1]
output_file = sys.argv[2]

genes = pd.read_csv(genes_file, sep="\t")
quant = pd.read_csv(quant_file, sep="\t")

rows = []

for _, row in genes.iterrows():
    gene = row["gene_symbol"]
    transcripts = row["transcript_ids"].split(",")

    sub = quant[quant["Name"].isin(transcripts)]

    rows.append({
        "gene_symbol": gene,
        "transcripts_found": len(sub),
        "total_TPM": sub["TPM"].sum(),
        "total_NumReads": sub["NumReads"].sum()
    })

out = pd.DataFrame(rows)

Path(output_file).parent.mkdir(parents=True, exist_ok=True)

out.to_csv(output_file, sep="\t", index=False)

print(out)
print(f"\nSaved to: {output_file}")
