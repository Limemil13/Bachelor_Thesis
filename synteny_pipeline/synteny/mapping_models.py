from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class NeighborHit:
    neighbor_name: str
    method: str
    hit_found: bool

    sseqid: str | None = None
    start_1based: int | None = None
    end_1based: int | None = None
    strand: str | None = None
    pident: float | None = None
    length: int | None = None
    evalue: float | None = None
    bitscore: float | None = None
    query_protein_length: int | None = None


@dataclass(frozen=True)
class NeighborMappingResult:
    ref_assembly: str
    target_assembly: str
    gene_name: str
    flank: int
    n: int
    target_contig: str
    target_region_start_1based: int
    target_region_end_1based: int
    hits: list[NeighborHit]

    def to_json(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path

    @staticmethod
    def from_json(path: Path) -> NeighborMappingResult:
        data = json.loads(path.read_text(encoding="utf-8"))
        hits = [NeighborHit(**h) for h in data["hits"]]
        data["hits"] = hits
        return NeighborMappingResult(**data)
