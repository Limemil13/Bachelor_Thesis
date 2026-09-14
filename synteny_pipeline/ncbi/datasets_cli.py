from __future__ import annotations

import shutil
import subprocess
import zipfile
from dataclasses import dataclass
from pathlib import Path


class DatasetsCliNotFound(RuntimeError):
    pass


@dataclass(frozen=True)
class GenomeDownloadResult:
    assembly: str
    zip_path: Path
    unzip_dir: Path


def require_datasets_cli() -> str:
    exe = shutil.which("datasets")
    if exe:
        return exe

    candidates = [
        Path.home() / "miniconda3" / "bin" / "datasets.exe",
        Path.home() / "miniconda3" / "Scripts" / "datasets.exe",
        Path.home() / "miniconda3" / "Library" / "bin" / "datasets.exe",
    ]

    for path in candidates:
        if path.exists():
            return str(path)

    raise DatasetsCliNotFound(
        "NCBI 'datasets' CLI not found. Checked PATH and common miniconda locations."
    )


def download_genome_package(
    assembly: str,
    out_dir: Path,
    include_gff3: bool = True,
    include_genome_fasta: bool = True,
    force: bool = False,
) -> GenomeDownloadResult:
    out_dir = out_dir.resolve()
    raw_dir = out_dir / "raw"
    unzip_dir = out_dir / "unzipped"

    raw_dir.mkdir(parents=True, exist_ok=True)
    unzip_dir.mkdir(parents=True, exist_ok=True)

    zip_path = raw_dir / f"{assembly}.zip"

    if zip_path.exists() and not force:
        if not any(unzip_dir.iterdir()):
            unzip_file(zip_path, unzip_dir)
        return GenomeDownloadResult(
            assembly=assembly, zip_path=zip_path, unzip_dir=unzip_dir
        )

    datasets_exe = require_datasets_cli()

    cmd = [
        datasets_exe,
        "download",
        "genome",
        "accession",
        assembly,
        "--filename",
        str(zip_path),
    ]

    if include_gff3:
        cmd += ["--include", "gff3"]
    if include_genome_fasta:
        cmd += ["--include", "genome"]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "NCBI datasets download failed.\n\n"
            f"Command:\n  {' '.join(cmd)}\n\n"
            f"STDOUT:\n{proc.stdout}\n\n"
            f"STDERR:\n{proc.stderr}\n"
        )

    unzip_file(zip_path, unzip_dir)
    return GenomeDownloadResult(
        assembly=assembly, zip_path=zip_path, unzip_dir=unzip_dir
    )


def unzip_file(zip_path: Path, unzip_dir: Path) -> None:
    with zipfile.ZipFile(zip_path, "r") as handle:
        handle.extractall(unzip_dir)
