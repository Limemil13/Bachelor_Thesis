from pathlib import Path
import csv


BASE = Path("analyses/slrp_phylogeny")

MSA_IN = BASE / "04_alignments" / "msa_quality_summary.tsv"
DOMAIN_IN = BASE / "03_domains" / "domain_scan" / "domain_summary_by_protein.tsv"

OUT_DIR = BASE / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MSA_OUT = OUT_DIR / "msa_quality_table_clean.tsv"
DOMAIN_OUT = OUT_DIR / "domain_quality_table_clean.tsv"
COMBINED_OUT = OUT_DIR / "combined_domain_msa_quality_table.tsv"


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def msa_quality_label(gap_percent: float):
    if gap_percent >= 40:
        return "Needs inspection"
    if gap_percent >= 25:
        return "Watch"
    return "Good"


def domain_quality_label(domain_status: str, n_lrr: int):
    if domain_status == "SLRP_like" and n_lrr >= 2:
        return "Good SLRP-like domain architecture"
    if domain_status == "SLRP_like":
        return "Weak but SLRP-like"
    return "Needs inspection"


def main():
    msa_rows = read_tsv(MSA_IN)
    domain_rows = read_tsv(DOMAIN_IN)

    # ---------- clean MSA table ----------
    clean_msa = []

    for r in msa_rows:
        gap_percent = float(r["gap_percent"])

        clean_msa.append({
            "Gene": r["gene"],
            "Species": r["target"],
            "Protein accession": r["hit_id"],
            "Alignment file": r["alignment_file"],
            "Protein length": r["ungapped_length"],
            "Alignment length": r["alignment_length"],
            "Gap percentage": r["gap_percent"],
            "MSA quality": msa_quality_label(gap_percent),
            "Raw MSA status": r["msa_status"],
        })

    msa_fields = [
        "Gene",
        "Species",
        "Protein accession",
        "Alignment file",
        "Protein length",
        "Alignment length",
        "Gap percentage",
        "MSA quality",
        "Raw MSA status",
    ]

    write_tsv(MSA_OUT, clean_msa, msa_fields)

    # ---------- clean domain table ----------
    clean_domain = []

    for r in domain_rows:
        n_lrr = int(r["n_lrr_domains"])

        clean_domain.append({
            "Gene": r["gene"],
            "Species": r["target"],
            "Protein accession": r["hit_id"],
            "Protein length": r["protein_len"],
            "Number of significant Pfam domains": r["n_significant_domains"],
            "Number of LRR-related domains": r["n_lrr_domains"],
            "Detected domain families": r["domain_names"],
            "Domain quality": domain_quality_label(r["domain_status"], n_lrr),
            "Raw domain status": r["domain_status"],
        })

    domain_fields = [
        "Gene",
        "Species",
        "Protein accession",
        "Protein length",
        "Number of significant Pfam domains",
        "Number of LRR-related domains",
        "Detected domain families",
        "Domain quality",
        "Raw domain status",
    ]

    write_tsv(DOMAIN_OUT, clean_domain, domain_fields)

    # ---------- combined domain + MSA table ----------
    domain_by_key = {}

    for r in clean_domain:
        key = (r["Gene"], r["Species"], r["Protein accession"])
        domain_by_key[key] = r

    combined = []

    for m in clean_msa:
        key = (m["Gene"], m["Species"], m["Protein accession"])
        d = domain_by_key.get(key)

        if d is None:
            continue

        combined.append({
            "Gene": m["Gene"],
            "Species": m["Species"],
            "Protein accession": m["Protein accession"],
            "Protein length": m["Protein length"],
            "Gap percentage": m["Gap percentage"],
            "MSA quality": m["MSA quality"],
            "Number of LRR-related domains": d["Number of LRR-related domains"],
            "Domain quality": d["Domain quality"],
            "Detected domain families": d["Detected domain families"],
        })

    combined_fields = [
        "Gene",
        "Species",
        "Protein accession",
        "Protein length",
        "Gap percentage",
        "MSA quality",
        "Number of LRR-related domains",
        "Domain quality",
        "Detected domain families",
    ]

    write_tsv(COMBINED_OUT, combined, combined_fields)

    print("Written:")
    print(MSA_OUT)
    print(DOMAIN_OUT)
    print(COMBINED_OUT)


if __name__ == "__main__":
    main()