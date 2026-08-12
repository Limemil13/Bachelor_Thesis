import re
from collections import defaultdict

fasta = "transcriptome/Mus_musculus.GRCm39.cdna.all.fa"
output = "genes_expanded.tsv"

genes_to_find = [
    # SLRPs
    "Dcn", "Fmod", "Aspn", "Bgn", "Lum", "Epyc", "Ogn", "Omd",
    "Prelp", "Kera", "Chad", "Chadl", "Tsku", "Podn", "Podnl1", "Ecm2",

    # cartilage / growth plate controls
    "Col2a1", "Acan", "Sox9", "Comp", "Matn1",
    "Col9a1", "Col9a2", "Col11a1", "Col10a1", "Runx2", "Mmp13",

    # broad ECM controls
    "Col1a1", "Col1a2", "Col3a1", "Fn1", "Vim",

    # non-cartilage tissue controls
    "Alb", "Myh6", "Mbp", "Krt14", "Acta1"
]

gene_data = defaultdict(lambda: {"gene_id": None, "transcripts": []})

with open(fasta, "r", encoding="utf-8") as f:
    for line in f:
        if not line.startswith(">"):
            continue

        transcript_match = re.match(r">(\S+)", line)
        gene_symbol_match = re.search(r"gene_symbol:([^\s]+)", line)
        gene_id_match = re.search(r"gene:(ENSMUSG[0-9]+\.[0-9]+)", line)
        biotype_match = re.search(r"transcript_biotype:([^\s]+)", line)

        if not transcript_match or not gene_symbol_match or not gene_id_match:
            continue

        transcript_id = transcript_match.group(1)
        gene_symbol = gene_symbol_match.group(1)
        gene_id = gene_id_match.group(1)
        biotype = biotype_match.group(1) if biotype_match else ""

        if gene_symbol not in genes_to_find:
            continue

        # Prefer protein-coding transcripts and avoid retained introns
        if "retained_intron" in biotype:
            continue

        gene_data[gene_symbol]["gene_id"] = gene_id
        gene_data[gene_symbol]["transcripts"].append(transcript_id)

with open(output, "w", encoding="utf-8") as out:
    out.write("gene_symbol\tgene_id\ttranscript_ids\n")
    for gene in genes_to_find:
        if gene in gene_data and gene_data[gene]["transcripts"]:
            transcripts = ",".join(gene_data[gene]["transcripts"])
            out.write(f"{gene}\t{gene_data[gene]['gene_id']}\t{transcripts}\n")
        else:
            out.write(f"{gene}\tNOT_FOUND\tNOT_FOUND\n")

print(f"Saved: {output}")
