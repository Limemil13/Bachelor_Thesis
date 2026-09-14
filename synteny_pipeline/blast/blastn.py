from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


class BlastNotFound(RuntimeError):
    pass


@dataclass(frozen=True)
class BlastnHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    evalue: float
    bitscore: float
    sstart: int
    send: int

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def start_1based(self) -> int:
        return min(self.sstart, self.send)

    @property
    def end_1based(self) -> int:
        return max(self.sstart, self.send)


def _require_blastn() -> str:
    try:
        subprocess.run(
            ["blastn", "-version"], check=True, capture_output=True, text=True
        )
        return "blastn"
    except FileNotFoundError:
        raise BlastNotFound("blastn not found on PATH")


def run_blastn_top1(query_fna: Path, db_prefix: Path) -> BlastnHit | None:
    # Return the best blastn hit, or None if no hit is found
    blastn = _require_blastn()
    outfmt = "6 qseqid sseqid pident length evalue bitscore sstart send"

    cmd = [
        blastn,
        "-query",
        str(query_fna),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-max_hsps",
        "1",
        "-evalue",
        "1e-10",
        "-outfmt",
        outfmt,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{result.stderr}")

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        return None

    best_parts = None
    best_bitscore = None

    for line in lines:
        parts = line.split("\t")
        if len(parts) != 8:
            raise RuntimeError(f"Unexpected blastn output:\n{line}")

        bitscore = float(parts[5])
        if best_bitscore is None or bitscore > best_bitscore:
            best_bitscore = bitscore
            best_parts = parts

    return BlastnHit(
        qseqid=best_parts[0],
        sseqid=best_parts[1],
        pident=float(best_parts[2]),
        length=int(best_parts[3]),
        evalue=float(best_parts[4]),
        bitscore=float(best_parts[5]),
        sstart=int(best_parts[6]),
        send=int(best_parts[7]),
    )
