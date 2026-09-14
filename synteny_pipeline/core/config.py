from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PipelineConfig:
    flank_bp: int = 50_000
    neighbors_each_side: int = 10
    blast_mode: str = "local"
    evolution_mode: str = "auto"


def default_project_root() -> Path:
    return Path.cwd()
