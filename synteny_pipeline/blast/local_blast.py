from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class BlastNotFound(RuntimeError):
    pass


@dataclass(frozen=True)
class BlastHit:
    sseqid: str
    sstart: int
    send: int
    pident: float
    length: int
    evalue: float
    bitscore: float

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def start_1based(self) -> int:
        return min(self.sstart, self.send)

    @property
    def end_1based(self) -> int:
        return max(self.sstart, self.send)


def require_blast_tools() -> None:
    for exe in ("makeblastdb", "blastn", "tblastn", "blastp"):
        if shutil.which(exe) is None:
            raise BlastNotFound(
                f"NCBI BLAST+ tool '{exe}' not found on PATH.\n"
                f"Install BLAST+ and restart PyCharm so the terminal sees it."
            )


def make_nucl_db(genome_fasta: Path, db_prefix: Path) -> None:
    require_blast_tools()
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "makeblastdb",
        "-in",
        str(genome_fasta),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
    ]
    subprocess.run(cmd, check=True)


def run_tblastn_top1(query_faa: Path, db_prefix: Path) -> BlastHit | None:

    # Protein query vs nucleotide genome (typical for finding a gene locus when only genome exists).

    require_blast_tools()
    outfmt = "6 sseqid sstart send pident length evalue bitscore"
    cmd = [
        "tblastn",
        "-query",
        str(query_faa),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-outfmt",
        outfmt,
    ]
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    line = p.stdout.strip().splitlines()
    if not line:
        return None
    parts = line[0].split()
    return BlastHit(
        sseqid=parts[0],
        sstart=int(parts[1]),
        send=int(parts[2]),
        pident=float(parts[3]),
        length=int(parts[4]),
        evalue=float(parts[5]),
        bitscore=float(parts[6]),
    )
