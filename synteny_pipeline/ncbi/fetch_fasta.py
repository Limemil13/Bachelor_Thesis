from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError


class FastaFetchError(RuntimeError):
    pass


@dataclass(frozen=True)
class FastaFetchResult:
    accession: str
    fasta_path: Path


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _http_get(
    url: str, timeout: int = 60, retries: int = 4, backoff: float = 1.5
) -> str:
    last_error = None

    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except HTTPError as e:
            last_error = e
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(backoff**i)
                continue
            raise
        except URLError as e:
            last_error = e
            time.sleep(backoff**i)

    raise FastaFetchError(f"HTTP request failed after retries: {last_error}")


def _get_json(url: str) -> dict:
    text = _http_get(url)
    try:
        return json.loads(text)
    except Exception:
        raise FastaFetchError(f"Could not parse JSON. First 200 chars:\n{text[:200]}")


def _esearch_uid(db: str, term: str) -> str:
    params = {
        "db": db,
        "term": term,
        "retmode": "json",
        "retmax": "1",
    }
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    data = _get_json(url)
    ids = data.get("esearchresult", {}).get("idlist", [])

    if not ids:
        raise FastaFetchError(f"No UID found for {db} term={term}")

    return str(ids[0])


def _elink_uids(dbfrom: str, db: str, uid: str) -> list[str]:
    params = {
        "dbfrom": dbfrom,
        "db": db,
        "id": uid,
        "retmode": "json",
    }
    url = f"{EUTILS}/elink.fcgi?" + urllib.parse.urlencode(params)
    data = _get_json(url)

    out = []
    for linkset in data.get("linksets", []):
        for linkdb in linkset.get("linksetdbs", []):
            for linked_uid in linkdb.get("links", []):
                out.append(str(linked_uid))

    return out


def _esummary_accessionversion(db: str, uid: str) -> str:
    params = {"db": db, "id": uid, "retmode": "json"}
    url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
    data = _get_json(url)

    record = data.get("result", {}).get(str(uid), {})
    accver = record.get("accessionversion")
    if not accver:
        raise FastaFetchError(f"No accessionversion returned for {db} uid={uid}")

    return accver


def _resolve_to_protein_accession(accession: str) -> str:
    accession = accession.strip()

    if accession.startswith(("NP_", "XP_", "YP_", "WP_")):
        return accession

    if accession.startswith(("NM_", "XM_", "NR_", "XR_")):
        nuccore_uid = _esearch_uid("nuccore", f"{accession}[Accession]")
        protein_uids = _elink_uids("nuccore", "protein", nuccore_uid)
        if not protein_uids:
            raise FastaFetchError(f"{accession} has no linked protein in NCBI.")
        return _esummary_accessionversion("protein", protein_uids[0])

    return accession


def fetch_protein_fasta(accession: str, out_fasta: Path) -> FastaFetchResult:
    """Fetch a protein FASTA from NCBI E-utilities."""
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    protein_accession = _resolve_to_protein_accession(accession)

    params = {
        "db": "protein",
        "id": protein_accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params)

    try:
        text = _http_get(url)
    except Exception as e:
        raise FastaFetchError(f"Could not fetch FASTA for {protein_accession}: {e}")

    if not text.strip().startswith(">"):
        raise FastaFetchError(
            f"NCBI did not return FASTA for {protein_accession}. First 200 chars:\n{text[:200]}"
        )

    out_fasta.write_text(text, encoding="utf-8")
    return FastaFetchResult(accession=protein_accession, fasta_path=out_fasta)
