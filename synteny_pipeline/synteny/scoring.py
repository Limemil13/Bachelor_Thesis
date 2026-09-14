from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ScoreConfig:
    min_pident: float = 30.0
    max_evalue: float = 1e-5
    min_align_len: int = 60


@dataclass(frozen=True)
class NeighborScoreRow:
    neighbor_name: str
    hit_found: bool
    pident: float | None
    evalue: float | None
    length: int | None
    pass_thresholds: bool


@dataclass(frozen=True)
class MappingScoreResult:
    mapping_json: Path
    ref_assembly: str
    target_assembly: str
    gene_name: str
    thresholds: ScoreConfig
    total_neighbors: int
    hits_found: int
    hits_passing: int
    score_fraction_passing: float
    rows: list[NeighborScoreRow]

    def to_json(self, out_path: Path) -> Path:
        out_path.parent.mkdir(parents=True, exist_ok=True)

        payload: dict[str, Any] = {
            "mapping_json": str(self.mapping_json),
            "ref_assembly": self.ref_assembly,
            "target_assembly": self.target_assembly,
            "gene_name": self.gene_name,
            "thresholds": {
                "min_pident": self.thresholds.min_pident,
                "max_evalue": self.thresholds.max_evalue,
                "min_align_len": self.thresholds.min_align_len,
            },
            "total_neighbors": self.total_neighbors,
            "hits_found": self.hits_found,
            "hits_passing": self.hits_passing,
            "score_fraction_passing": self.score_fraction_passing,
            "rows": [
                {
                    "neighbor_name": row.neighbor_name,
                    "hit_found": row.hit_found,
                    "pident": row.pident,
                    "evalue": row.evalue,
                    "length": row.length,
                    "pass_thresholds": row.pass_thresholds,
                }
                for row in self.rows
            ],
        }

        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return out_path

    def to_tsv(self, out_path: Path) -> Path:
        out_path.parent.mkdir(parents=True, exist_ok=True)

        lines = ["neighbor_name\thit_found\tpident\tevalue\tlength\tpass_thresholds\n"]
        for row in self.rows:
            lines.append(
                f"{row.neighbor_name}\t{int(row.hit_found)}\t"
                f"{'' if row.pident is None else row.pident}\t"
                f"{'' if row.evalue is None else row.evalue}\t"
                f"{'' if row.length is None else row.length}\t"
                f"{int(row.pass_thresholds)}\n"
            )

        out_path.write_text("".join(lines), encoding="utf-8")
        return out_path


def passes_thresholds(
    cfg: ScoreConfig, pident: float, evalue: float, length: int
) -> bool:
    return (
        pident >= cfg.min_pident
        and evalue <= cfg.max_evalue
        and length >= cfg.min_align_len
    )


def score_mapping_json(mapping_json: Path, cfg: ScoreConfig) -> MappingScoreResult:
    mapping_json = mapping_json.resolve()
    data = json.loads(mapping_json.read_text(encoding="utf-8"))

    rows: list[NeighborScoreRow] = []
    hits_found = 0
    hits_passing = 0

    for hit in data.get("hits", []):
        neighbor_name = hit.get("neighbor_name", "unknown")
        hit_found = bool(hit.get("hit_found", False))
        pident = hit.get("pident")
        evalue = hit.get("evalue")
        length = hit.get("length")

        if hit_found:
            hits_found += 1

        passed = False
        if (
            hit_found
            and pident is not None
            and evalue is not None
            and length is not None
        ):
            passed = passes_thresholds(cfg, float(pident), float(evalue), int(length))
            if passed:
                hits_passing += 1

        rows.append(
            NeighborScoreRow(
                neighbor_name=neighbor_name,
                hit_found=hit_found,
                pident=None if pident is None else float(pident),
                evalue=None if evalue is None else float(evalue),
                length=None if length is None else int(length),
                pass_thresholds=passed,
            )
        )

    total_neighbors = len(rows)
    score_fraction = hits_passing / total_neighbors if total_neighbors else 0.0

    return MappingScoreResult(
        mapping_json=mapping_json,
        ref_assembly=data["ref_assembly"],
        target_assembly=data["target_assembly"],
        gene_name=data["gene_name"],
        thresholds=cfg,
        total_neighbors=total_neighbors,
        hits_found=hits_found,
        hits_passing=hits_passing,
        score_fraction_passing=score_fraction,
        rows=rows,
    )
