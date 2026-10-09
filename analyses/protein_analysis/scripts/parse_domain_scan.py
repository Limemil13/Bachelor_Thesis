"""Parse the canonical Pfam hmmscan domtblout into thesis-facing TSV files."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from canonical_dataset import MANIFEST

BASE = Path(__file__).resolve().parents[1]
SCAN_DIR = BASE / "domains" / "domain_scan_canonical"
DOMTBLOUT = SCAN_DIR / "domain_candidates_canonical.pfam.domtblout"
HITS_OUT = SCAN_DIR / "domain_hits_canonical.tsv"
SUMMARY_OUT = SCAN_DIR / "domain_summary_by_protein_canonical.tsv"
I_EVALUE_THRESHOLD = 1e-3


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_lrr(domain_name: str, description: str) -> bool:
    # Pfam we catching all the names, actually only 3lrrrr :)
    text = f"{domain_name} {description}".lower()
    return (
        "lrr" in domain_name.lower()
        or "leucine rich repeat" in text
        or "leucine-rich repeat" in text
    )


def parse_domtblout() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with DOMTBLOUT.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip().split(maxsplit=22)
            if len(parts) < 23:
                raise ValueError(f"Malformed domtblout line: {raw_line.rstrip()}")

            domain_name = parts[0]
            domain_acc = parts[1]
            protein_id = parts[3]
            protein_len = int(parts[5])
            i_evalue = float(parts[12])
            if i_evalue > I_EVALUE_THRESHOLD:
                continue

            id_parts = protein_id.split("|")
            if len(id_parts) != 3:
                raise ValueError(
                    f"Expected gene|species|accession query ID: {protein_id}"
                )
            gene, species, accession = id_parts
            description = parts[22]

            rows.append(
                {
                    "gene": gene,
                    "species": species,
                    "accession": accession,
                    "protein_id": protein_id,
                    "protein_len": protein_len,
                    "domain_name": domain_name,
                    "domain_acc": domain_acc,
                    "description": description,
                    "i_evalue": f"{i_evalue:.6g}",
                    "domain_score": parts[13],
                    "ali_from": parts[17],
                    "ali_to": parts[18],
                    "env_from": parts[19],
                    "env_to": parts[20],
                    "is_lrr": "yes" if is_lrr(domain_name, description) else "no",
                }
            )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    #if < 2 LRR then its flagged
    if not DOMTBLOUT.exists():
        raise SystemExit(f"Missing hmmscan output: {DOMTBLOUT}")

    manifest = read_manifest()
    hits = parse_domtblout()
    canonical_ids = {record["canonical_id"] for record in manifest}
    hits = [hit for hit in hits if hit["protein_id"] in canonical_ids]
    hit_fields = [
        "gene",
        "species",
        "accession",
        "protein_id",
        "protein_len",
        "domain_name",
        "domain_acc",
        "description",
        "i_evalue",
        "domain_score",
        "ali_from",
        "ali_to",
        "env_from",
        "env_to",
        "is_lrr",
    ]
    hits.sort(
        key=lambda row: (
            str(row["gene"]),
            str(row["species"]),
            int(row["ali_from"]),
            str(row["domain_name"]),
        )
    )
    write_tsv(HITS_OUT, hits, hit_fields)

    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        by_id[str(hit["protein_id"])].append(hit)

    summaries: list[dict[str, object]] = []
    for record in manifest:
        protein_hits = by_id.get(record["canonical_id"], [])
        lrr_hits = [hit for hit in protein_hits if hit["is_lrr"] == "yes"]
        names = sorted({str(hit["domain_name"]) for hit in protein_hits})
        lrr_names = sorted({str(hit["domain_name"]) for hit in lrr_hits})
        summaries.append(
            {
                "gene": record["gene"],
                "species": record["species"],
                "accession": record["accession"],
                "protein_len": record["analyzed_length"],
                "n_significant_domains": len(protein_hits),
                "n_lrr_domains": len(lrr_hits),
                "domain_names": ",".join(names),
                "lrr_domain_names": ",".join(lrr_names),
                "domain_status": "SLRP_like" if len(lrr_hits) >= 2 else "check",
            }
        )

    summary_fields = [
        "gene",
        "species",
        "accession",
        "protein_len",
        "n_significant_domains",
        "n_lrr_domains",
        "domain_names",
        "lrr_domain_names",
        "domain_status",
    ]
    write_tsv(SUMMARY_OUT, summaries, summary_fields)
    print(f"Wrote {HITS_OUT} ({len(hits)} significant domains)")
    print(f"Wrote {SUMMARY_OUT} ({len(summaries)} proteins)")


if __name__ == "__main__":
    main()
