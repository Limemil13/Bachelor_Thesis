"""Fetch reproducible Bgee 16 expression evidence for thesis SLRP candidates.

The script keeps anatomy-only calls (all supported data types) separate from
stage-resolved calls supported by bulk or single-cell RNA-seq.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(__file__).resolve().parent
RAW_DIR = BASE / "bgee" / "raw"
TABLE_DIR = BASE / "tables"
API = "https://www.bgee.org/api/"

SPECIES = {
    "mouse": {
        "taxon_id": 10090,
        "symbols": {
            "BGN": "Bgn",
            "EPYC": "Epyc",
            "FMOD": "Fmod",
            "OGN": "Ogn",
            "OMD": "Omd",
            "PRELP": "Prelp",
            "DCN": "Dcn",
            "LUM": "Lum",
            "ASPN": "Aspn",
        },
    },
    "chicken": {
        "taxon_id": 9031,
        "symbols": {
            "BGN": "BGN",
            "EPYC": "EPYC",
            "FMOD": "FMOD",
            "OGN": "OGN",
            "OMD": "OMD",
            "PRELP": "PRELP",
            "DCN": "DCN",
            "LUM": "LUM",
            "ASPN": "ASPN",
        },
    },
    "zebrafish": {
        "taxon_id": 7955,
        "symbols": {
            "BGN": "bgna",
            "EPYC": "epyc",
            "FMOD": "fmodb",
            "OGN": "ogna",
            "OMD": "omd",
            "PRELP": "prelp",
            "DCN": "dcn",
            "LUM": "lum",
            "ASPN": "aspn",
        },
    },
}
# Terms used to flag potentially relevant musculoskeletal anatomy calls.
SKELETAL_PATTERN = re.compile(
    r"cartilage|chondr|growth plate|epiphys|joint|bone|skelet|tendon|ligament|"
    r"intervertebral|notochord|sclerotome|somite|femur|tibia|humerus|vertebr",
    re.IGNORECASE,
)


def request_json(params: list[tuple[str, str]]) -> tuple[dict, str]:
    url = API + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(
        url, headers={"User-Agent": "SLRP-thesis-reproducibility/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = json.load(response)
    if payload.get("code") != 200 or payload.get("status") != "SUCCESS":
        raise RuntimeError(f"Bgee request failed: {url}: {payload}")
    return payload, url


def exact_gene_match(
    species: str, gene: str, symbol: str, taxon_id: int
) -> dict[str, str]:
    payload, url = request_json(
        [
            ("display_type", "json"),
            ("page", "gene"),
            ("query", symbol),
            ("species_id", str(taxon_id)),
            ("limit", "20"),
        ]
    )
    matches = payload["data"]["result"]["geneMatches"]
    exact = [
        item["gene"]
        for item in matches
        if item["gene"]["name"].casefold() == symbol.casefold()
        and int(item["gene"]["species"]["id"]) == taxon_id
    ]
    if len(exact) > 1:
        raise ValueError(f"Ambiguous exact Bgee matches for {species} {gene}/{symbol}")
    selected = exact[0] if exact else None
    return {
        "species": species,
        "taxon_id": str(taxon_id),
        "thesis_gene": gene,
        "queried_symbol": symbol,
        "bgee_gene_id": selected["geneId"] if selected else "",
        "bgee_gene_name": selected["name"] if selected else "",
        "bgee_description": selected.get("description", "") if selected else "",
        "mapping_status": "exact_symbol_match"
        if selected
        else "not_found_exact_symbol",
        "nonexact_search_hits": ";".join(
            f"{item['gene']['geneId']}:{item['gene']['name']}:{item.get('matchSource', '')}"
            for item in matches
        ),
        "gene_search_url": url,
    }


def expression_query(gene_row: dict[str, str], mode: str) -> tuple[dict, str]:
    params = [
        ("display_type", "json"),
        ("page", "gene"),
        ("action", "expression"),
        ("gene_id", gene_row["bgee_gene_id"]),
        ("species_id", gene_row["taxon_id"]),
        ("cond_param", "anat_entity"),
        ("expr_type", "expressed"),
    ]
    if mode == "stage_rnaseq":
        params.extend(
            [
                ("cond_param", "dev_stage"),
                ("data_type", "RNA_SEQ"),
                ("data_type", "SC_RNA_SEQ"),
            ]
        )
    payload, url = request_json(params)
    return payload, url


def flatten_call(
    call: dict, gene_row: dict[str, str], mode: str, url: str
) -> dict[str, str]:
    # Flatten the nested response while retaining IDs and request provenance.
    condition = call.get("condition", {})
    anatomy = condition.get("anatEntity") or {}
    stage = condition.get("devStage") or {}
    score = call.get("expressionScore") or {}
    data_types = call.get("dataTypesWithData") or []
    anatomy_name = anatomy.get("name", "")
    return {
        "species": gene_row["species"],
        "taxon_id": gene_row["taxon_id"],
        "thesis_gene": gene_row["thesis_gene"],
        "bgee_gene_id": gene_row["bgee_gene_id"],
        "bgee_gene_name": gene_row["bgee_gene_name"],
        "query_mode": mode,
        "anatomical_entity_id": anatomy.get("id", ""),
        "anatomical_entity_name": anatomy_name,
        "developmental_stage_id": stage.get("id", ""),
        "developmental_stage_name": stage.get("name", ""),
        "expression_score": str(score.get("expressionScore", "")),
        "expression_score_confidence": score.get("expressionScoreConfidence", ""),
        "fdr": str(call.get("fdr", "")),
        "expression_quality": call.get("expressionQuality", ""),
        "supporting_data_types": ";".join(data_types),
        "skeletal_cartilage_term_match": "yes"
        if SKELETAL_PATTERN.search(anatomy_name)
        else "no",
        "request_url": url,
    }


def write_tsv(
    path: Path, rows: list[dict[str, str]], fields: list[str] | None = None
) -> None:
    # Keep a stable TSV column order.
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    # Resolve exact Bgee gene IDs before requesting expression calls.
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    # Run independent lookups concurrently, then sort the results.
    lookup_tasks = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        for species, info in SPECIES.items():
            for gene, symbol in info["symbols"].items():
                lookup_tasks.append(
                    executor.submit(
                        exact_gene_match, species, gene, symbol, int(info["taxon_id"])
                    )
                )
        gene_rows = [future.result() for future in as_completed(lookup_tasks)]
    gene_rows.sort(key=lambda row: (row["species"], row["thesis_gene"]))
    write_tsv(TABLE_DIR / "bgee_slrp_gene_mapping.tsv", gene_rows)

    # Retain raw JSON, flattened calls, and query-level provenance.
    calls: list[dict[str, str]] = []
    provenance: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(expression_query, gene_row, mode): (gene_row, mode)
            for gene_row in gene_rows
            if gene_row["bgee_gene_id"]
            for mode in ("anatomy_all_data", "stage_rnaseq")
        }
        for future in as_completed(futures):
            gene_row, mode = futures[future]
            payload, url = future.result()
            raw_path = (
                RAW_DIR / f"{gene_row['species']}_{gene_row['thesis_gene']}_{mode}.json"
            )
            raw_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            result_calls = payload.get("data", {}).get("calls", [])
            calls.extend(
                flatten_call(call, gene_row, mode, url) for call in result_calls
            )
            provenance.append(
                {
                    "species": gene_row["species"],
                    "thesis_gene": gene_row["thesis_gene"],
                    "query_mode": mode,
                    "call_count": str(len(result_calls)),
                    "raw_json": str(raw_path.relative_to(BASE.parent.parent)).replace(
                        "\\", "/"
                    ),
                    "request_url": url,
                    "retrieved_date": "2026-08-22",
                    "bgee_release": "16.0",
                }
            )
            time.sleep(0.05)

    # Sort after concurrent retrieval so output does not depend on completion time.
    calls.sort(
        key=lambda row: (
            row["species"],
            row["thesis_gene"],
            row["query_mode"],
            -float(row["expression_score"] or 0),
            row["anatomical_entity_name"],
        )
    )
    call_fields = [
        "species",
        "taxon_id",
        "thesis_gene",
        "bgee_gene_id",
        "bgee_gene_name",
        "query_mode",
        "anatomical_entity_id",
        "anatomical_entity_name",
        "developmental_stage_id",
        "developmental_stage_name",
        "expression_score",
        "expression_score_confidence",
        "fdr",
        "expression_quality",
        "supporting_data_types",
        "skeletal_cartilage_term_match",
        "request_url",
    ]
    write_tsv(TABLE_DIR / "bgee_slrp_expression_calls_all.tsv", calls, call_fields)
    # Keyword filtering is a convenience subset, not proof of biological absence.
    relevant = [row for row in calls if row["skeletal_cartilage_term_match"] == "yes"]
    write_tsv(
        TABLE_DIR / "bgee_slrp_skeletal_cartilage_calls.tsv", relevant, call_fields
    )
    provenance.sort(
        key=lambda row: (row["species"], row["thesis_gene"], row["query_mode"])
    )
    write_tsv(TABLE_DIR / "bgee_slrp_query_provenance.tsv", provenance)

    # Summarize mappings and call counts while keeping missing mappings explicit.
    summary: list[dict[str, str]] = []
    for gene_row in gene_rows:
        for mode in ("anatomy_all_data", "stage_rnaseq"):
            if not gene_row["bgee_gene_id"]:
                summary.append(
                    {
                        "species": gene_row["species"],
                        "thesis_gene": gene_row["thesis_gene"],
                        "bgee_gene_id": "",
                        "query_mode": mode,
                        "total_expressed_call_count": "",
                        "skeletal_cartilage_match_count": "",
                        "best_skeletal_cartilage_expression_score": "",
                        "best_skeletal_cartilage_term": "",
                        "best_skeletal_cartilage_stage": "",
                        "best_skeletal_cartilage_quality": "",
                        "best_skeletal_cartilage_data_types": "",
                        "interpretation": (
                            "gene not resolved by exact symbol in Bgee; missing database evidence, "
                            "not evidence of biological absence"
                        ),
                    }
                )
                continue
            subset = [
                row
                for row in calls
                if row["species"] == gene_row["species"]
                and row["thesis_gene"] == gene_row["thesis_gene"]
                and row["query_mode"] == mode
            ]
            matched = [
                row for row in subset if row["skeletal_cartilage_term_match"] == "yes"
            ]
            best = max(
                matched,
                key=lambda row: float(row["expression_score"] or 0),
                default=None,
            )
            summary.append(
                {
                    "species": gene_row["species"],
                    "thesis_gene": gene_row["thesis_gene"],
                    "bgee_gene_id": gene_row["bgee_gene_id"],
                    "query_mode": mode,
                    "total_expressed_call_count": str(len(subset)),
                    "skeletal_cartilage_match_count": str(len(matched)),
                    "best_skeletal_cartilage_expression_score": best["expression_score"]
                    if best
                    else "",
                    "best_skeletal_cartilage_term": best["anatomical_entity_name"]
                    if best
                    else "",
                    "best_skeletal_cartilage_stage": best["developmental_stage_name"]
                    if best
                    else "",
                    "best_skeletal_cartilage_quality": best["expression_quality"]
                    if best
                    else "",
                    "best_skeletal_cartilage_data_types": best["supporting_data_types"]
                    if best
                    else "",
                    "interpretation": (
                        "curated expression support found"
                        if matched
                        else "no matching Bgee call; this is missing evidence, not evidence of absence"
                    ),
                }
            )
    write_tsv(TABLE_DIR / "bgee_slrp_expression_summary.tsv", summary)

    print(f"Mapped genes: {len(gene_rows)}")
    print(f"Expression queries: {len(provenance)}")
    print(f"Expression calls: {len(calls)}")
    print(f"Skeletal/cartilage term matches: {len(relevant)}")


if __name__ == "__main__":
    main()
