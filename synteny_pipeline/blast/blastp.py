from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


class BlastNotFound(RuntimeError):
    pass


@dataclass(frozen=True)
class BlastpHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    evalue: float
    bitscore: float


def _require_exe(name: str) -> str:
    try:
        subprocess.run([name, "-version"], check=True, capture_output=True, text=True)
        return name
    except FileNotFoundError:
        raise BlastNotFound(f"{name} not found on PATH")


def make_prot_db(protein_faa: Path, db_prefix: Path):
    makeblastdb = _require_exe("makeblastdb")
    db_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        makeblastdb,
        "-in",
        str(protein_faa),
        "-dbtype",
        "prot",
        "-out",
        str(db_prefix),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed:\n{result.stderr}")


def run_blastp_top1(query_faa: Path, db_prefix: Path) -> BlastpHit | None:
    # Return the top blastp hit, or None if there is no hit
    blastp = _require_exe("blastp")
    outfmt = "6 qseqid sseqid pident length evalue bitscore"

    cmd = [
        blastp,
        "-query",
        str(query_faa),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-outfmt",
        outfmt,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"blastp failed:\n{result.stderr}")

    line = result.stdout.strip()
    if not line:
        return None

    parts = line.split("\t")
    return BlastpHit(
        qseqid=parts[0],
        sseqid=parts[1],
        pident=float(parts[2]),
        length=int(parts[3]),
        evalue=float(parts[4]),
        bitscore=float(parts[5]),
    )
