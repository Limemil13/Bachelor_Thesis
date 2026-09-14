from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @property
    def data_dir(self) -> Path:
        return self.root / "data"

    @property
    def ncbi_dir(self) -> Path:
        return self.data_dir / "ncbi"

    @property
    def genomes_dir(self) -> Path:
        return self.ncbi_dir / "genomes"

    @property
    def genes_dir(self) -> Path:
        return self.ncbi_dir / "genes"

    @property
    def blast_dir(self) -> Path:
        return self.data_dir / "blast"

    @property
    def runs_dir(self) -> Path:
        return self.root / "runs"

    def ensure(self) -> None:
        for path in [
            self.data_dir,
            self.ncbi_dir,
            self.genomes_dir,
            self.genes_dir,
            self.blast_dir,
            self.runs_dir,
        ]:
            path.mkdir(parents=True, exist_ok=True)

    def new_run_dir(self) -> Path:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = self.runs_dir / timestamp
        run_dir.mkdir(parents=True, exist_ok=False)
        return run_dir
