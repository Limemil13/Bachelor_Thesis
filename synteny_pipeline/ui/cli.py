from __future__ import annotations

import csv
import gzip
import json
import shutil
import statistics
import subprocess
import textwrap
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import typer
from pyfaidx import Fasta
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from synteny_pipeline.blast.blastn import BlastNotFound as BlastnNotFound
from synteny_pipeline.blast.blastn import run_blastn_top1
from synteny_pipeline.blast.blastp import BlastNotFound as BlastpNotFound
from synteny_pipeline.blast.blastp import make_prot_db, run_blastp_top1
from synteny_pipeline.blast.local_blast import (
    BlastNotFound,
    make_nucl_db,
    run_tblastn_top1,
)
from synteny_pipeline.core.config import PipelineConfig, default_project_root
from synteny_pipeline.core.paths import ProjectPaths
from synteny_pipeline.genome.find_gene import find_gene_by_name_fast
from synteny_pipeline.genome.gff_region import genes_in_region_fast
from synteny_pipeline.genome.protein_ids import protein_id_for_gene_symbol_first
from synteny_pipeline.genome.region_fasta import (
    ContigNotFoundError,
    InvalidCoordinatesError,
    RegionSpec,
    extract_region_to_fasta,
)
from synteny_pipeline.ncbi.datasets_cli import (
    DatasetsCliNotFound,
    download_genome_package,
)
from synteny_pipeline.ncbi.fetch_fasta import FastaFetchError, fetch_protein_fasta
from synteny_pipeline.ncbi.find_files import find_genome_files
from synteny_pipeline.phylogeny.project import create_phylogeny_run
from synteny_pipeline.synteny.compare import compare_profiles, profile_from_block
from synteny_pipeline.synteny.mapping_models import NeighborHit, NeighborMappingResult
from synteny_pipeline.synteny.neighbors import neighbors_around_focal, pick_middle_gene
from synteny_pipeline.synteny.scoring import ScoreConfig, score_mapping_json

matplotlib.use("Agg")
app = typer.Typer(add_completion=False)
console = Console()


def _require_exe(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def read_fasta_records_multi(path: Path) -> list[tuple[str, str]]:
    """Return ``(header, sequence)`` records from a FASTA file."""
    records: list[tuple[str, str]] = []
    header = None
    seq_chunks: list[str] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks).upper()))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        records.append((header, "".join(seq_chunks).upper()))

    return records


def parse_candidate_header(header: str) -> dict[str, str]:
    """Parse the metadata fields in a candidate FASTA header."""
    parts = header.split()
    main = parts[0]

    fields = {
        "target": "",
        "query": "",
        "hit_id": "",
        "forward_pident": "",
        "forward_qcovs": "",
        "reverse_hit": "",
    }

    main_parts = main.split("|")
    if len(main_parts) >= 3:
        fields["target"] = main_parts[0]
        fields["query"] = main_parts[1]
        fields["hit_id"] = main_parts[2]

    for p in parts[1:]:
        if "=" not in p:
            continue
        key, value = p.split("=", 1)
        if key in fields:
            fields[key] = value

    return fields


def fasta_read_one(path: Path) -> tuple[str, str]:
    header = None
    seq_chunks = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].split()[0]
            else:
                seq_chunks.append(line)
    if header is None:
        raise ValueError(f"No FASTA header found in {path}")
    return header, "".join(seq_chunks).upper()


def fasta_sequence_length(path: Path) -> int:
    _, seq = fasta_read_one(path)
    return len(seq)


def fasta_write_one(path: Path, header: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = seq.replace(" ", "").replace("\n", "")
    wrapped = "\n".join(textwrap.wrap(seq, 80))
    path.write_text(f">{header}\n{wrapped}\n", encoding="utf-8")


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def extract_subseq_from_single_fasta(
    fa_path: Path, start_1based: int, end_1based: int, strand: str
) -> str:
    """
    Extract subsequence from a single-record FASTA.
    start/end are 1-based inclusive.
    strand is '+' or '-'.
    """
    _, seq = fasta_read_one(fa_path)
    if start_1based < 1 or end_1based > len(seq) or start_1based > end_1based:
        raise ValueError(
            f"Invalid slice {start_1based}-{end_1based} for length {len(seq)} in {fa_path}"
        )
    sub = seq[start_1based - 1 : end_1based]
    if strand == "-":
        sub = revcomp(sub)
    return sub


def run_blastx_top1(query_fna: Path, db_prefix: Path) -> dict[str, object] | None:
    """Return the top BLASTX hit, or ``None`` when no hit passes."""
    _require_exe("blastx")
    cmd = [
        "blastx",
        "-query",
        str(query_fna),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-outfmt",
        "6 qseqid sseqid pident length evalue bitscore sstart send",
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip() or "blastx failed")
    out = p.stdout.strip()
    if not out:
        return None
    cols = out.splitlines()[0].split("\t")
    return {
        "qseqid": cols[0],
        "sseqid": cols[1],
        "pident": float(cols[2]),
        "length": int(cols[3]),
        "evalue": float(cols[4]),
        "bitscore": float(cols[5]),
        "sstart": int(cols[6]),
        "send": int(cols[7]),
    }


def keep_neighbor_gene(symbol: str | None) -> bool:
    if not symbol:
        return False

    symbol = symbol.upper()
    if symbol.startswith("LOC"):
        return False
    if symbol in {"CCDST", "CRCT1"}:
        return False

    return True


def fasta_get_full_header_by_id(fasta: Path, wanted_id: str) -> str:
    """Return a protein's full FASTA header without the leading ``>``."""
    wanted_id = normalize_acc(wanted_id)

    with open_text_maybe_gzip(fasta) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                continue

            full_header = line[1:]
            first_id = normalize_acc(full_header.split()[0])

            if (
                first_id == wanted_id
                or first_id.split(".")[0] == wanted_id.split(".")[0]
            ):
                return full_header

    return ""


def find_protein_fasta(unzip_dir: Path) -> Path:
    """Find a plain or gzipped protein FASTA in an NCBI dataset folder."""
    candidates = list(unzip_dir.rglob("*.faa")) + list(unzip_dir.rglob("*.faa.gz"))
    if not candidates:
        raise FileNotFoundError(f"No protein FASTA (.faa/.faa.gz) found in {unzip_dir}")
    return candidates[0]


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def normalize_acc(x: str) -> str:
    return (
        x.split()[0]
        .replace("ref|", "")
        .replace("gb|", "")
        .replace("sp|", "")
        .split("|")[-1]
    )


def fasta_get_record_by_id(fasta: Path, wanted_id: str) -> tuple[str, str] | None:
    wanted_id = normalize_acc(wanted_id)
    header = None
    seq = []

    with open_text_maybe_gzip(fasta) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    current_id = normalize_acc(header)
                    if current_id == wanted_id:
                        return header, "".join(seq)

                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)

    if header is not None and normalize_acc(header) == wanted_id:
        return header, "".join(seq)

    return None


def run_blastp_table_top1(
    query_faa: Path,
    db_prefix: Path,
    evalue: float = 1e-5,
) -> dict[str, object] | None:
    """Return the top BLASTP hit with query coverage for tabular reports."""
    _require_exe("blastp")

    cmd = [
        "blastp",
        "-query",
        str(query_faa),
        "-db",
        str(db_prefix),
        "-max_target_seqs",
        "1",
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length evalue bitscore qcovs",
    ]

    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip() or "blastp failed")

    out = p.stdout.strip()
    if not out:
        return None

    c = out.splitlines()[0].split("\t")
    return {
        "qseqid": c[0],
        "sseqid": c[1],
        "pident": float(c[2]),
        "length": int(c[3]),
        "evalue": float(c[4]),
        "bitscore": float(c[5]),
        "qcovs": float(c[6]),
    }


def simplify_reverse_description(desc: str) -> str:
    """Remove the accession and organism suffix from a FASTA description."""
    if not desc:
        return ""

    parts = desc.split(maxsplit=1)
    if len(parts) == 2:
        desc = parts[1]

    if "[" in desc:
        desc = desc.split("[", 1)[0].strip()

    return desc


def interpret_blastp_row(row: dict[str, object]) -> str:
    """Classify a reciprocal BLASTP row using its reference description."""
    reverse_desc = row.get("reverse_description", "").lower()
    query = row.get("query", "").lower()

    if not row.get("forward_best_hit"):
        return "no forward hit"

    if not row.get("reverse_best_hit"):
        return "no reverse hit"

    if query and query in reverse_desc:
        return "likely ortholog"

    # Product names cover reference headers that omit the gene symbol.
    known_names = {
        "DCN": ["decorin"],
        "BGN": ["biglycan"],
        "FMOD": ["fibromodulin"],
        "ASPN": ["asporin"],
    }

    query_upper = row.get("query", "").upper()
    for synonym in known_names.get(query_upper, []):
        if synonym in reverse_desc:
            return "likely ortholog"

    return "unclear / possible paralog"


@app.command("extract-blastp-candidates")
def extract_blastp_candidates(
    blastp_tsv: Path = typer.Option(
        ..., "--blastp-tsv", exists=True, help="Path to blastp_ortholog.tsv"
    ),
    out_faa: Path | None = typer.Option(
        None, "--out-faa", help="Output candidate FASTA"
    ),
    root: Path = typer.Option(default_project_root(), "--root"),
):
    """Extract forward BLASTP best hits for downstream protein analyses."""
    paths = ProjectPaths(root=root)
    paths.ensure()

    rows = []
    with blastp_tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    if not rows:
        console.print("[bold red]No rows found in BLASTP TSV.[/bold red]")
        raise typer.Exit(code=1)

    query = rows[0].get("query", "QUERY")

    if out_faa is None:
        out_faa = blastp_tsv.parent / f"{query}_candidates.faa"

    written = 0
    missing = []

    with out_faa.open("w", encoding="utf-8") as out:
        for row in rows:
            target = row.get("target", "")
            hit_id = row.get("forward_best_hit", "")

            if not target or not hit_id:
                missing.append((target, hit_id, "missing target or forward_best_hit"))
                continue

            proteome_faa = (
                paths.data_dir / "ncbi" / "proteomes" / target / "protein.faa"
            )

            if not proteome_faa.exists():
                missing.append((target, hit_id, f"proteome missing: {proteome_faa}"))
                continue

            record = fasta_get_record_by_id(proteome_faa, hit_id)

            if record is None:
                missing.append((target, hit_id, "sequence not found in proteome"))
                continue

            header, seq = record

            clean_header = (
                f"{target}|{query}|{hit_id} "
                f"forward_pident={row.get('forward_pident', '')} "
                f"forward_qcovs={row.get('forward_qcovs', '')} "
                f"reverse_hit={row.get('reverse_best_hit', '')}"
            )

            out.write(f">{clean_header}\n")
            out.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

            written += 1

    console.print(
        Panel.fit(
            "[bold green]BLASTP candidates extracted[/bold green]\n\n"
            f"[bold]Query:[/bold] {query}\n"
            f"[bold]Input TSV:[/bold] {blastp_tsv}\n"
            f"[bold]Output FASTA:[/bold] {out_faa}\n"
            f"[bold]Sequences written:[/bold] {written}\n"
            f"[bold]Missing/problem rows:[/bold] {len(missing)}"
        )
    )

    if missing:
        console.print("[bold yellow]Missing/problem rows:[/bold yellow]")
        for target, hit_id, reason in missing[:20]:
            console.print(f"  {target} | {hit_id} | {reason}")


@app.command("summarize-candidates")
def summarize_candidates(
    candidates_faa: Path = typer.Option(
        ..., "--candidates-faa", exists=True, help="Candidate FASTA file"
    ),
):
    """Summarize candidate lengths and reciprocal-BLAST metadata."""
    records = read_fasta_records_multi(candidates_faa)

    if not records:
        console.print("[bold red]No FASTA records found.[/bold red]")
        raise typer.Exit(code=1)

    rows = []

    lengths = [len(seq) for _, seq in records]
    median_len = statistics.median(lengths)

    for header, seq in records:
        info = parse_candidate_header(header)
        length = len(seq)

        if length < 0.7 * median_len:
            length_status = "short"
        elif length > 1.3 * median_len:
            length_status = "long"
        else:
            length_status = "ok"

        rows.append(
            {
                "query": info["query"],
                "target": info["target"],
                "hit_id": info["hit_id"],
                "length": length,
                "median_length": median_len,
                "length_status": length_status,
                "forward_pident": info["forward_pident"],
                "forward_qcovs": info["forward_qcovs"],
                "reverse_hit": info["reverse_hit"],
            }
        )

    out_dir = candidates_faa.parent
    out_tsv = out_dir / "candidates_summary.tsv"
    out_txt = out_dir / "candidates_summary_readable.txt"

    headers = [
        "query",
        "target",
        "hit_id",
        "length",
        "median_length",
        "length_status",
        "forward_pident",
        "forward_qcovs",
        "reverse_hit",
    ]

    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    txt_lines = []
    txt_lines.append(
        f"{'query':<8} {'target':<25} {'hit_id':<18} "
        f"{'len':>6} {'status':<8} {'pid':>7} {'qcov':>7} {'reverse_hit':<18}\n"
    )
    txt_lines.append("-" * 110 + "\n")

    for r in rows:
        txt_lines.append(
            f"{r['query']:<8} "
            f"{r['target']:<25} "
            f"{r['hit_id']:<18} "
            f"{r['length']:>6} "
            f"{r['length_status']:<8} "
            f"{r['forward_pident']:>7} "
            f"{r['forward_qcovs']:>7} "
            f"{r['reverse_hit']:<18}\n"
        )

    out_txt.write_text("".join(txt_lines), encoding="utf-8")

    console.print(
        Panel.fit(
            "[bold green]Candidate FASTA summary written[/bold green]\n\n"
            f"[bold]Input:[/bold] {candidates_faa}\n"
            f"[bold]Sequences:[/bold] {len(records)}\n"
            f"[bold]Median length:[/bold] {median_len}\n"
            f"[bold]TSV:[/bold] {out_tsv}\n"
            f"[bold]Readable TXT:[/bold] {out_txt}"
        )
    )


@app.command("filter-candidates-by-length")
def filter_candidates_by_length(
    candidates_faa: Path = typer.Option(..., "--candidates-faa", exists=True),
    min_ratio: float = typer.Option(0.7, "--min-ratio"),
    max_ratio: float = typer.Option(1.3, "--max-ratio"),
):
    """Filter candidates by their length relative to the panel median."""
    records = read_fasta_records_multi(candidates_faa)

    if not records:
        console.print("[bold red]No FASTA records found.[/bold red]")
        raise typer.Exit(code=1)

    lengths = [len(seq) for _, seq in records]
    median_len = statistics.median(lengths)

    out_dir = candidates_faa.parent
    stem = candidates_faa.stem.replace("_candidates", "")
    out_faa = out_dir / f"{stem}_candidates_filtered.faa"
    out_excluded = out_dir / f"{stem}_candidates_excluded.tsv"

    kept = []
    excluded = []

    for header, seq in records:
        length = len(seq)
        low = min_ratio * median_len
        high = max_ratio * median_len

        info = parse_candidate_header(header)

        if low <= length <= high:
            kept.append((header, seq))
        else:
            excluded.append(
                {
                    "query": info.get("query", ""),
                    "target": info.get("target", ""),
                    "hit_id": info.get("hit_id", ""),
                    "length": length,
                    "median_length": median_len,
                    "reason": f"outside {min_ratio}-{max_ratio}x median length",
                }
            )

    with out_faa.open("w", encoding="utf-8") as f:
        for header, seq in kept:
            f.write(f">{header}\n")
            f.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

    with out_excluded.open("w", encoding="utf-8", newline="") as f:
        fieldnames = ["query", "target", "hit_id", "length", "median_length", "reason"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded)

    console.print(
        Panel.fit(
            "[bold green]Candidate FASTA filtered[/bold green]\n\n"
            f"[bold]Input:[/bold] {candidates_faa}\n"
            f"[bold]Median length:[/bold] {median_len}\n"
            f"[bold]Kept:[/bold] {len(kept)}\n"
            f"[bold]Excluded:[/bold] {len(excluded)}\n"
            f"[bold]Filtered FASTA:[/bold] {out_faa}\n"
            f"[bold]Excluded TSV:[/bold] {out_excluded}"
        )
    )


@app.command("prepare-domain-inputs")
def prepare_domain_inputs(
    root: Path = typer.Option(default_project_root(), "--root"),
):
    """Collect filtered candidate FASTAs for domain analysis."""
    paths = ProjectPaths(root=root)
    paths.ensure()

    gene_to_fasta = {
        "BGN": paths.runs_dir
        / "BGN_run"
        / "blastp_ortholog"
        / "BGN"
        / "BGN_candidates_filtered.faa",
        "FMOD": paths.runs_dir
        / "FMOD_run"
        / "blastp_ortholog"
        / "FMOD"
        / "FMOD_candidates_filtered.faa",
        "OGN": paths.runs_dir
        / "ogn_run"
        / "blastp_ortholog"
        / "OGN"
        / "OGN_candidates_filtered.faa",
        "EPYC": paths.runs_dir
        / "epyc_run"
        / "blastp_ortholog"
        / "EPYC"
        / "EPYC_candidates_filtered.faa",
        "PRELP": paths.runs_dir
        / "prelp_run"
        / "blastp_ortholog"
        / "PRELP"
        / "PRELP_candidates_filtered.faa",
    }

    out_dir = paths.runs_dir / "domain_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_faa = out_dir / "domain_candidates_all.faa"
    metadata_tsv = out_dir / "domain_candidates_metadata.tsv"

    metadata_rows = []
    total_written = 0
    missing_files = []

    with combined_faa.open("w", encoding="utf-8") as combined:
        for gene, fasta_path in gene_to_fasta.items():
            if not fasta_path.exists():
                missing_files.append(str(fasta_path))
                continue

            gene_out = out_dir / f"{gene}_domain_input.faa"
            shutil.copyfile(fasta_path, gene_out)

            records = read_fasta_records_multi(fasta_path)

            for header, seq in records:
                info = parse_candidate_header(header)

                target = info.get("target", "")
                hit_id = info.get("hit_id", "")

                new_header = f"{gene}|{target}|{hit_id}"

                combined.write(f">{new_header}\n")
                combined.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

                metadata_rows.append(
                    {
                        "gene": gene,
                        "target": target,
                        "hit_id": hit_id,
                        "length": len(seq),
                        "forward_pident": info.get("forward_pident", ""),
                        "forward_qcovs": info.get("forward_qcovs", ""),
                        "reverse_hit": info.get("reverse_hit", ""),
                        "original_header": header,
                    }
                )

                total_written += 1

    headers = [
        "gene",
        "target",
        "hit_id",
        "length",
        "forward_pident",
        "forward_qcovs",
        "reverse_hit",
        "original_header",
    ]

    with metadata_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata_rows)

    lines = [
        f"[bold]Output folder:[/bold] {out_dir}",
        f"[bold]Combined FASTA:[/bold] {combined_faa}",
        f"[bold]Metadata TSV:[/bold] {metadata_tsv}",
        f"[bold]Sequences written:[/bold] {total_written}",
    ]

    if missing_files:
        lines.append("\n[bold yellow]Missing filtered FASTA files:[/bold yellow]")
        for m in missing_files:
            lines.append(f"  {m}")

    console.print(
        Panel.fit(
            "[bold green]Domain-analysis input prepared[/bold green]\n\n"
            + "\n".join(lines)
        )
    )


@app.command()
def init(
    root: Path = typer.Option(
        default_project_root(),
        "--root",
        help="Project root folder where cache/output folders will be created.",
    ),
):
    paths = ProjectPaths(root=root)
    paths.ensure()

    cfg = PipelineConfig()
    console.print(
        Panel.fit(
            f"[green]Mapping + scoring done[/green]\n\n"
            f"[bold]Root:[/bold] {paths.root}\n"
            f"[bold]Genomes cache:[/bold] {paths.genomes_dir}\n"
            f"[bold]Runs output:[/bold] {paths.runs_dir}\n\n"
            f"[bold]Default settings:[/bold]\n"
            f"  flank_bp={cfg.flank_bp}\n"
            f"  neighbors_each_side={cfg.neighbors_each_side}\n"
            f"  blast_mode={cfg.blast_mode}\n"
            f"  evolution_mode={cfg.evolution_mode}\n"
        )
    )


@app.command("download-genome")
def download_genome(
    assembly: str = typer.Option(
        ..., "--assembly", help="Assembly accession like GCF_... or GCA_..."
    ),
    root: Path = typer.Option(
        default_project_root(),
        "--root",
        help="Project root folder (same as used in init).",
    ),
    force: bool = typer.Option(False, "--force", help="Redownload even if cached."),
):

    paths = ProjectPaths(root=root)
    paths.ensure()

    assembly_dir = paths.genomes_dir / assembly

    try:
        res = download_genome_package(
            assembly=assembly,
            out_dir=assembly_dir,
            include_gff3=True,
            include_genome_fasta=True,
            force=force,
        )
    except DatasetsCliNotFound as e:
        console.print(f"[bold red]{e}[/bold red]")
        raise typer.Exit(code=1)

    console.print(
        Panel.fit(
            f"[bold green]Downloaded genome package[/bold green]\n\n"
            f"[bold]Assembly:[/bold] {res.assembly}\n"
            f"[bold]Zip:[/bold] {res.zip_path}\n"
            f"[bold]Unzipped:[/bold] {res.unzip_dir}\n"
        )
    )


@app.command("download-proteome")
def download_proteome(
    assembly: str = typer.Option(
        ..., "--assembly", help="Assembly accession like GCF_..."
    ),
    name: str = typer.Option(
        ..., "--name", help="Clean species name, e.g. human, mouse"
    ),
    root: Path = typer.Option(default_project_root(), "--root"),
    force: bool = typer.Option(False, "--force"),
):
    """
    Download only the protein FASTA from NCBI Datasets and store it cleanly in:
    data/ncbi/proteomes/<name>/protein.faa
    """
    paths = ProjectPaths(root=root)
    paths.ensure()

    proteome_dir = paths.data_dir / "ncbi" / "proteomes" / name
    final_faa = proteome_dir / "protein.faa"

    if final_faa.exists() and not force:
        console.print(
            Panel.fit(
                "[bold green]Proteome already exists[/bold green]\n\n"
                f"[bold]Name:[/bold] {name}\n"
                f"[bold]Protein FASTA:[/bold] {final_faa}\n\n"
                f"Use --force to redownload."
            )
        )
        return

    tmp_dir = proteome_dir / "_tmp"
    zip_path = proteome_dir / f"{name}.zip"

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    proteome_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        assembly,
        "--include",
        "protein",
        "--filename",
        str(zip_path),
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        console.print(
            f"[bold red]NCBI datasets download failed for {assembly}[/bold red]"
        )
        raise typer.Exit(code=1)

    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(tmp_dir)
    except Exception as e:
        console.print(f"[bold red]Failed to unzip proteome:[/bold red] {e}")
        raise typer.Exit(code=1)

    faa_files = list(tmp_dir.rglob("*.faa")) + list(tmp_dir.rglob("*.faa.gz"))

    if not faa_files:
        console.print(
            f"[bold red]No protein FASTA found in downloaded package for {assembly}[/bold red]"
        )
        raise typer.Exit(code=1)

    source_faa = faa_files[0]

    if final_faa.exists():
        final_faa.unlink()

    if source_faa.suffix == ".gz":
        import gzip

        with gzip.open(source_faa, "rt", encoding="utf-8", errors="replace") as fin:
            final_faa.write_text(fin.read(), encoding="utf-8")
    else:
        shutil.copy2(source_faa, final_faa)

    shutil.rmtree(tmp_dir)

    console.print(
        Panel.fit(
            "[bold green]Proteome downloaded[/bold green]\n\n"
            f"[bold]Species name:[/bold] {name}\n"
            f"[bold]Assembly:[/bold] {assembly}\n"
            f"[bold]Protein FASTA:[/bold] {final_faa}"
        )
    )


@app.command("init-phylogeny-run")
def init_phylogeny_run(
    run_name: str = typer.Option("slrp_phylogeny", "--run-name"),
    root: Path = typer.Option(default_project_root(), "--root"),
):
    """Create the directory structure for an SLRP phylogeny run."""
    folders = create_phylogeny_run(root=root, run_name=run_name)

    lines = ["[bold green]Phylogeny run folder created[/bold green]\n"]
    for name, path in folders.items():
        lines.append(f"[bold]{name}:[/bold] {path}")

    console.print(Panel.fit("\n".join(lines)))


@app.command("inspect-genome")
def inspect_genome(
    assembly: str = typer.Option(
        ..., "--assembly", help="Assembly accession like GCF_... or GCA_..."
    ),
    root: Path = typer.Option(
        default_project_root(),
        "--root",
        help="Project root folder (same as used in init).",
    ),
):

    paths = ProjectPaths(root=root)
    paths.ensure()

    assembly_dir = paths.genomes_dir / assembly
    unzip_dir = assembly_dir / "unzipped"

    files = find_genome_files(unzip_dir)

    gff_line = (
        f"[bold]GFF3:[/bold] {files.gff3}"
        if files.gff3
        else "[bold]GFF3:[/bold] (not found)"
    )
    ar_line = (
        f"[bold]Assembly report:[/bold] {files.assembly_report}"
        if files.assembly_report
        else "[bold]Assembly report:[/bold] (not found)"
    )
    dc_line = (
        f"[bold]Dataset catalog:[/bold] {files.dataset_catalog}"
        if files.dataset_catalog
        else "[bold]Dataset catalog:[/bold] (not found)"
    )

    console.print(
        Panel.fit(
            f"[bold green]Genome files found[/bold green]\n\n"
            f"[bold]Assembly:[/bold] {assembly}\n"
            f"[bold]Genome FASTA:[/bold] {files.genome_fasta}\n"
            f"{gff_line}\n"
            f"{ar_line}\n"
            f"{dc_line}\n"
        )
    )


@app.command("extract-region")
def extract_region(
    assembly: str = typer.Option(
        ..., "--assembly", help="Assembly accession like GCF_... or GCA_..."
    ),
    contig: str = typer.Option(
        ...,
        "--contig",
        help="Contig/chromosome name as in the FASTA (e.g., NC_..., chr1, etc.)",
    ),
    start: int = typer.Option(
        ..., "--start", help="1-based start coordinate (inclusive)"
    ),
    end: int = typer.Option(..., "--end", help="1-based end coordinate (inclusive)"),
    root: Path = typer.Option(
        default_project_root(), "--root", help="Project root folder."
    ),
):
    paths = ProjectPaths(root=root)
    paths.ensure()

    # Locate FASTA inside this assembly cache
    unzip_dir = paths.genomes_dir / assembly / "unzipped"
    files = find_genome_files(unzip_dir)
    genome_fasta = files.genome_fasta

    # Output into a new run folder
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = paths.runs_dir / run_id / "regions" / assembly
    out_fasta = out_dir / f"{assembly}_{contig}_{start}_{end}.fasta"

    region = RegionSpec(contig=contig, start_1based=start, end_1based=end)

    try:
        written = extract_region_to_fasta(genome_fasta, region, out_fasta)
    except (ContigNotFoundError, InvalidCoordinatesError) as e:
        console.print(f"[bold red]{e}[/bold red]")
        raise typer.Exit(code=1)

    console.print(
        Panel.fit(
            f"[bold green]Region extracted[/bold green]\n\n"
            f"[bold]Assembly:[/bold] {assembly}\n"
            f"[bold]Genome FASTA:[/bold] {genome_fasta}\n"
            f"[bold]Region:[/bold] {contig}:{start}-{end} (1-based)\n"
            f"[bold]Output FASTA:[/bold] {written}\n"
        )
    )


@app.command("list-genes")
def list_genes(
    assembly: str = typer.Option(
        ..., "--assembly", help="Assembly accession like GCF_... or GCA_..."
    ),
    contig: str = typer.Option(
        ..., "--contig", help="Contig/chromosome name as in the FASTA/GFF"
    ),
    start: int = typer.Option(
        ..., "--start", help="1-based start coordinate (inclusive)"
    ),
    end: int = typer.Option(..., "--end", help="1-based end coordinate (inclusive)"),
    root: Path = typer.Option(
        default_project_root(), "--root", help="Project root folder."
    ),
):
    paths = ProjectPaths(root=root)
    paths.ensure()

    unzip_dir = paths.genomes_dir / assembly / "unzipped"
    files = find_genome_files(unzip_dir)

    if not files.gff3:
        console.print(
            "[bold red]No GFF3 annotation found for this assembly.[/bold red]"
        )
        raise typer.Exit(code=1)

    try:
        genes = genes_in_region_fast(
            gff_path=files.gff3,
            contig=contig,
            start_1based=start,
            end_1based=end,
        )
    except Exception as e:
        console.print(f"[bold red]Failed to parse GFF region:[/bold red] {e}")
        raise typer.Exit(code=1)

    if not genes:
        console.print(
            Panel.fit(
                f"[yellow]No gene features found in region[/yellow]\n\n"
                f"[bold]Assembly:[/bold] {assembly}\n"
                f"[bold]Region:[/bold] {contig}:{start}-{end}"
            )
        )
        return

    lines = []
    lines.append(f"Found {len(genes)} gene(s) in {contig}:{start}-{end}\n")

    for g in genes[:80]:
        name = g.gene_name or "-"
        gid = g.gene_id or "-"
        lines.append(
            f"{g.start_1based:>10}-{g.end_1based:<10} {g.strand}  {name}  (ID={gid})"
        )

    if len(genes) > 80:
        lines.append(f"\n... showing first 80 of {len(genes)}")

    console.print(
        Panel.fit("[bold green]Genes in region[/bold green]\n\n" + "\n".join(lines))
    )


@app.command("neighbors-in-region")
def neighbors_in_region(
    assembly: str = typer.Option(..., "--assembly"),
    contig: str = typer.Option(..., "--contig"),
    start: int = typer.Option(..., "--start"),
    end: int = typer.Option(..., "--end"),
    n: int = typer.Option(2, "--n", help="Number of neighbors upstream/downstream"),
    root: Path = typer.Option(default_project_root(), "--root"),
):

    # Pick a focal gene automatically and show its upstream/downstream neighbors.

    paths = ProjectPaths(root=root)
    paths.ensure()

    unzip_dir = paths.genomes_dir / assembly / "unzipped"
    files = find_genome_files(unzip_dir)

    genes = genes_in_region_fast(
        gff_path=files.gff3,
        contig=contig,
        start_1based=start,
        end_1based=end,
    )

    if len(genes) < 2:
        console.print("[red]Not enough genes in region.[/red]")
        raise typer.Exit(code=1)

    focal = pick_middle_gene(genes)
    block = neighbors_around_focal(genes, focal, n_each_side=n)

    lines = []
    lines.append(f"[bold]Focal gene:[/bold] {focal.gene_name or '-'}\n")

    lines.append("[bold]Upstream:[/bold]")
    for g in block.upstream:
        lines.append(
            f"  {g.start_1based}-{g.end_1based} {g.strand} {g.gene_name or '-'}"
        )

    lines.append("\n[bold]Downstream:[/bold]")
    for g in block.downstream:
        lines.append(
            f"  {g.start_1based}-{g.end_1based} {g.strand} {g.gene_name or '-'}"
        )

    console.print(
        Panel.fit("[bold green]Neighbor block[/bold green]\n\n" + "\n".join(lines))
    )


@app.command("neighbors-around-gene")
def neighbors_around_gene(
    assembly: str = typer.Option(..., "--assembly"),
    gene_name: str = typer.Option(
        ...,
        "--gene-name",
        help="Gene name as it appears in the GFF (e.g., foxo3, tp53)",
    ),
    flank: int = typer.Option(
        50000, "--flank", help="Bases upstream and downstream to include"
    ),
    n: int = typer.Option(10, "--n", help="Neighbors upstream and downstream"),
    root: Path = typer.Option(default_project_root(), "--root"),
):

    # Find a gene by name in the GFF, extract ±flank region, and return N upstream/downstream neighbors.

    paths = ProjectPaths(root=root)
    paths.ensure()

    unzip_dir = paths.genomes_dir / assembly / "unzipped"
    files = find_genome_files(unzip_dir)

    if not files.gff3:
        console.print(
            "[bold red]No GFF3 annotation found for this assembly.[/bold red]"
        )
        raise typer.Exit(code=1)

    found = find_gene_by_name_fast(files.gff3, gene_name=gene_name)
    if not found:
        console.print(
            Panel.fit(
                f"[bold red]Gene not found in GFF[/bold red]\n\n"
                f"[bold]Assembly:[/bold] {assembly}\n"
                f"[bold]Gene name:[/bold] {gene_name}\n\n"
                "Tip: gene names differ across species; use accession or BLAST-based mapping when needed."
            )
        )
        raise typer.Exit(code=1)

    # Compute region around gene
    region_start = max(1, found.start_1based - flank)
    region_end = found.end_1based + flank

    genes = genes_in_region_fast(
        gff_path=files.gff3,
        contig=found.contig,
        start_1based=region_start,
        end_1based=region_end,
    )

    # Find the focal gene within that region list (match by gene_id if possible, else by name)
    focal = None
    for g in genes:
        if found.gene_id and g.gene_id == found.gene_id:
            focal = g
            break
    if focal is None:
        for g in genes:
            if g.gene_name and g.gene_name.lower() == gene_name.lower():
                focal = g
                break
    if focal is None:
        console.print(
            "[bold red]Internal error: focal gene not found inside its own region list.[/bold red]"
        )
        raise typer.Exit(code=1)

    block = neighbors_around_focal(genes, focal, n_each_side=n)

    lines = []
    lines.append(f"[bold]Assembly:[/bold] {assembly}")
    lines.append(
        f"[bold]Focal gene:[/bold] {focal.gene_name or '-'} (ID={focal.gene_id or '-'})"
    )
    lines.append(
        f"[bold]Locus:[/bold] {found.contig}:{found.start_1based}-{found.end_1based} ({found.strand})"
    )
    lines.append(
        f"[bold]Region:[/bold] {found.contig}:{region_start}-{region_end} (±{flank} bp)\n"
    )

    lines.append("[bold]Upstream:[/bold]")
    for g in block.upstream:
        lines.append(
            f"  {g.start_1based}-{g.end_1based} {g.strand}  {g.gene_name or '-'} (ID={g.gene_id or '-'})"
        )

    lines.append("\n[bold]Downstream:[/bold]")
    for g in block.downstream:
        lines.append(
            f"  {g.start_1based}-{g.end_1based} {g.strand}  {g.gene_name or '-'} (ID={g.gene_id or '-'})"
        )

    console.print(
        Panel.fit(
            "[bold green]Neighbors around gene[/bold green]\n\n" + "\n".join(lines)
        )
    )


@app.command("compare-synteny")
def compare_synteny(
    assemblies: str = typer.Option(
        ..., "--assemblies", help="Comma-separated assemblies, e.g. GCF_...,..."
    ),
    gene_name: str = typer.Option(
        ..., "--gene-name", help="Gene name to search in each assembly GFF"
    ),
    flank: int = typer.Option(50000, "--flank"),
    n: int = typer.Option(10, "--n"),
    root: Path = typer.Option(default_project_root(), "--root"),
):

    # For multiple assemblies, find the focal gene by name, extract ±flank, pick N neighbors, and compare neighbor overlap across assemblies.

    paths = ProjectPaths(root=root)
    paths.ensure()

    assembly_list = [a.strip() for a in assemblies.split(",") if a.strip()]
    if len(assembly_list) < 2:
        console.print("[red]Provide at least 2 assemblies in --assemblies[/red]")
        raise typer.Exit(code=1)

    profiles = []
    failed = []

    for assembly in assembly_list:
        unzip_dir = paths.genomes_dir / assembly / "unzipped"

        if not unzip_dir.exists():
            failed.append(
                (
                    assembly,
                    "Not downloaded (missing unzipped folder). Run download-genome first.",
                )
            )
            continue

        try:
            files = find_genome_files(unzip_dir)
        except Exception as e:
            failed.append((assembly, f"Could not locate genome files: {e}"))
            continue

        if not files.gff3:
            failed.append((assembly, "No GFF3"))
            continue

        found = find_gene_by_name_fast(files.gff3, gene_name=gene_name)
        if not found:
            failed.append((assembly, f"Gene '{gene_name}' not found"))
            continue

        region_start = max(1, found.start_1based - flank)
        region_end = found.end_1based + flank

        genes = genes_in_region_fast(
            gff_path=files.gff3,
            contig=found.contig,
            start_1based=region_start,
            end_1based=region_end,
        )

        # Identify focal gene record inside region list
        focal = None
        for g in genes:
            if found.gene_id and g.gene_id == found.gene_id:
                focal = g
                break
        if focal is None:
            for g in genes:
                if g.gene_name and g.gene_name.lower() == gene_name.lower():
                    focal = g
                    break
        if focal is None:
            failed.append((assembly, "Focal gene not found inside region list"))
            continue

        block = neighbors_around_focal(genes, focal, n_each_side=n)
        profiles.append(profile_from_block(assembly, block))

    if len(profiles) < 2:
        console.print(
            Panel.fit(
                "[bold red]Not enough successful assemblies to compare[/bold red]\n\n"
                + "\n".join([f"{a}: {msg}" for a, msg in failed])
                if failed
                else ""
            )
        )
        raise typer.Exit(code=1)

    comp = compare_profiles(profiles)

    # Print a readable report
    lines = []
    lines.append(f"[bold]Gene:[/bold] {gene_name}")
    lines.append(
        f"[bold]flank:[/bold] {flank}  [bold]neighbors each side:[/bold] {n}\n"
    )

    lines.append("[bold]Per-assembly neighbor sets:[/bold]")
    for p in profiles:
        lines.append(f"\n[bold]{p.assembly}[/bold]")
        lines.append(
            f"  upstream:   {', '.join(p.upstream) if p.upstream else '(none)'}"
        )
        lines.append(
            f"  downstream: {', '.join(p.downstream) if p.downstream else '(none)'}"
        )

    lines.append("\n[bold]Pairwise overlap counts (neighbors shared):[/bold]")
    for (a, b), k in comp.pairwise_overlap.items():
        lines.append(f"  {a} vs {b}: {k}")

    if failed:
        lines.append("\n[bold yellow]Failed assemblies:[/bold yellow]")
        for a, msg in failed:
            lines.append(f"  {a}: {msg}")

    console.print(
        Panel.fit("[bold green]Synteny comparison[/bold green]\n\n" + "\n".join(lines))
    )


@app.command("map-by-blast")
def map_by_blast(
    assembly: str = typer.Option(
        ..., "--assembly", help="Target assembly (genome to search)"
    ),
    query_accession: str | None = typer.Option(
        None, "--query-accession", help="Protein accession (NP_/XP_/...)"
    ),
    query_faa: Path | None = typer.Option(
        None, "--query-faa", exists=True, help="Protein FASTA of focal gene"
    ),
    flank: int = typer.Option(
        50000, "--flank", help="Extract +/- flank bp around best hit"
    ),
    n: int = typer.Option(10, "--n", help="Neighbors each side (if annotation exists)"),
    root: Path = typer.Option(default_project_root(), "--root"),
):

    # Use tblastn (protein->genome) to locate the best hit, then extract +/- flank and list genes in that region.

    paths = ProjectPaths(root=root)
    paths.ensure()

    # Choose query source (accession OR fasta)
    if (query_accession is None) == (query_faa is None):
        console.print(
            "[bold red]Provide exactly one of --query-accession OR --query-faa[/bold red]"
        )
        raise typer.Exit(code=1)

    if query_accession is not None:
        try:
            query_faa = fetch_protein_fasta(
                accession=query_accession,
                out_fasta=paths.data_dir / "queries" / f"{query_accession}.faa",
            ).fasta_path
        except FastaFetchError as e:
            console.print(f"[bold red]{e}[/bold red]")
            raise typer.Exit(code=1)

    unzip_dir = paths.genomes_dir / assembly / "unzipped"
    files = find_genome_files(unzip_dir)

    if not files.genome_fasta:
        console.print("[bold red]No genome FASTA found for this assembly.[/bold red]")
        raise typer.Exit(code=1)

    # Build BLAST DB (cached)
    db_dir = paths.genomes_dir / assembly / "blastdb"
    db_prefix = db_dir / "genome_nucl"

    try:
        if not (db_dir.exists() and any(db_dir.glob("genome_nucl.*"))):
            make_nucl_db(files.genome_fasta, db_prefix)
    except BlastNotFound as e:
        console.print(f"[bold red]{e}[/bold red]")
        raise typer.Exit(code=1)

    hit = run_tblastn_top1(query_faa=query_faa, db_prefix=db_prefix)
    if hit is None:
        console.print(
            Panel.fit(
                f"[bold red]No BLAST hit found[/bold red]\n\n"
                f"Assembly: {assembly}\nQuery: {query_faa}"
            )
        )
        raise typer.Exit(code=1)

    region_start = max(1, hit.start_1based - flank)
    region_end = hit.end_1based + flank

    # Extract the hit-centered target region.
    run_dir = paths.new_run_dir()
    out_fa = (
        run_dir
        / "regions"
        / assembly
        / f"{assembly}_{hit.sseqid}_{region_start}_{region_end}.fasta"
    )
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    region = RegionSpec(
        contig=hit.sseqid, start_1based=region_start, end_1based=region_end
    )

    extract_region_to_fasta(
        genome_fasta=files.genome_fasta,
        region=region,
        out_fasta=out_fa,
    )

    # Annotate the region and select the feature nearest the hit midpoint.
    msg_lines = [
        f"[bold]Assembly:[/bold] {assembly}",
        f"[bold]Best hit:[/bold] {hit.sseqid}:{hit.start_1based}-{hit.end_1based} ({hit.strand})",
        f"[bold]Region:[/bold] {hit.sseqid}:{region_start}-{region_end} (±{flank})",
        f"[bold]Extracted FASTA:[/bold] {out_fa}",
    ]

    if files.gff3:
        genes = genes_in_region_fast(files.gff3, hit.sseqid, region_start, region_end)
        msg_lines.append(f"\n[bold]Genes in region:[/bold] {len(genes)}")

        # pick focal = closest gene to hit midpoint
        mid = (hit.start_1based + hit.end_1based) // 2
        focal = (
            min(genes, key=lambda g: abs(((g.start_1based + g.end_1based) // 2) - mid))
            if genes
            else None
        )

        if focal:
            block = neighbors_around_focal(genes, focal, n_each_side=n)
            msg_lines.append(
                f"\n[bold]Focal (closest gene):[/bold] {focal.gene_name or focal.gene_id}"
            )

            msg_lines.append("[bold]Upstream:[/bold]")
            for g in block.upstream:
                msg_lines.append(
                    f"  {g.start_1based}-{g.end_1based} {g.strand} {g.gene_name or g.gene_id}"
                )

            msg_lines.append("[bold]Downstream:[/bold]")
            for g in block.downstream:
                msg_lines.append(
                    f"  {g.start_1based}-{g.end_1based} {g.strand} {g.gene_name or g.gene_id}"
                )
        else:
            msg_lines.append(
                "[yellow]No gene features found in this region GFF.[/yellow]"
            )
    else:
        msg_lines.append(
            "[yellow]No GFF available: only the region FASTA was extracted.[/yellow]"
        )

    console.print(
        Panel.fit(
            "[bold green]BLAST mapping result[/bold green]\n\n" + "\n".join(msg_lines)
        )
    )


@app.command("map-neighbors-by-blastn")
def map_neighbors_by_blastn(
    ref_assembly: str = typer.Option(..., "--ref-assembly", help="Reference assembly"),
    target_assembly: str = typer.Option(
        ..., "--target-assembly", help="Target assembly"
    ),
    gene_name: str = typer.Option(..., "--gene-name", help="Focal gene name"),
    flank: int = typer.Option(50000, "--flank", help="Flank size around focal gene"),
    n: int = typer.Option(10, "--n", help="Neighbors on each side"),
    root: Path = typer.Option(default_project_root(), "--root"),
):
    # Map reference neighbor gene DNA into the target focal region using BLASTN

    paths = ProjectPaths(root=root)
    paths.ensure()

    # reference files
    ref_unzip = paths.genomes_dir / ref_assembly / "unzipped"
    ref_files = find_genome_files(ref_unzip)

    if not ref_files.gff3 or not ref_files.genome_fasta:
        console.print("[red]Reference assembly needs both GFF3 and genome FASTA.[/red]")
        raise typer.Exit(code=1)

    ref_found = find_gene_by_name_fast(ref_files.gff3, gene_name=gene_name)
    if not ref_found:
        console.print(
            Panel.fit(
                f"[red]Reference focal gene not found[/red]\n\n"
                f"Assembly: {ref_assembly}\n"
                f"Gene: {gene_name}"
            )
        )
        raise typer.Exit(code=1)

    ref_start = max(1, ref_found.start_1based - flank)
    ref_end = ref_found.end_1based + flank

    ref_genes = genes_in_region_fast(
        gff_path=ref_files.gff3,
        contig=ref_found.contig,
        start_1based=ref_start,
        end_1based=ref_end,
    )

    ref_focal = None
    for gene in ref_genes:
        if ref_found.gene_id and gene.gene_id == ref_found.gene_id:
            ref_focal = gene
            break

    if ref_focal is None:
        for gene in ref_genes:
            if gene.gene_name and gene.gene_name.lower() == gene_name.lower():
                ref_focal = gene
                break

    if ref_focal is None:
        console.print(
            "[red]Could not recover focal gene inside reference region.[/red]"
        )
        raise typer.Exit(code=1)

    ref_block = neighbors_around_focal(ref_genes, ref_focal, n_each_side=n)
    ref_neighbors = list(ref_block.upstream) + list(ref_block.downstream)

    if not ref_neighbors:
        console.print("[red]No neighbors found in reference block.[/red]")
        raise typer.Exit(code=1)

    # target files
    tgt_unzip = paths.genomes_dir / target_assembly / "unzipped"
    tgt_files = find_genome_files(tgt_unzip)

    if not tgt_files.genome_fasta:
        console.print("[red]Target assembly needs a genome FASTA.[/red]")
        raise typer.Exit(code=1)

    if not tgt_files.gff3:
        console.print("[red]Target assembly has no GFF3. Use map-by-blast first.[/red]")
        raise typer.Exit(code=1)

    tgt_found = find_gene_by_name_fast(tgt_files.gff3, gene_name=gene_name)
    if not tgt_found:
        console.print(
            Panel.fit(
                f"[red]Target focal gene not found[/red]\n\n"
                f"Assembly: {target_assembly}\n"
                f"Gene: {gene_name}"
            )
        )
        raise typer.Exit(code=1)

    tgt_start = max(1, tgt_found.start_1based - flank)
    tgt_end = tgt_found.end_1based + flank

    run_dir = paths.new_run_dir()
    out_dir = (
        run_dir
        / "neighbor_mapping"
        / f"{ref_assembly}_to_{target_assembly}"
        / gene_name
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    target_region_fa = out_dir / (
        f"target_region_{target_assembly}_{tgt_found.contig}_{tgt_start}_{tgt_end}.fna"
    )

    extract_region_to_fasta(
        genome_fasta=tgt_files.genome_fasta,
        region=RegionSpec(
            contig=tgt_found.contig,
            start_1based=tgt_start,
            end_1based=tgt_end,
        ),
        out_fasta=target_region_fa,
    )

    # BLAST db for target region
    db_dir = out_dir / "blastdb"
    db_prefix = db_dir / "target_region_db"

    try:
        if not (db_dir.exists() and any(db_dir.glob("target_region_db.*"))):
            make_nucl_db(target_region_fa, db_prefix)
    except (BlastNotFound, BlastnNotFound) as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)

    lines = [
        f"[bold]Reference:[/bold] {ref_assembly} (focal={gene_name})",
        f"[bold]Target:[/bold] {target_assembly}",
        f"[bold]Target region:[/bold] {tgt_found.contig}:{tgt_start}-{tgt_end} (±{flank})",
        "",
    ]

    hit_count = 0

    for gene in ref_neighbors:
        label = gene.gene_name or gene.gene_id or "unknown"
        query_fa = (
            out_dir
            / "queries"
            / f"{label}_{gene.contig}_{gene.start_1based}_{gene.end_1based}.fna"
        )
        query_fa.parent.mkdir(parents=True, exist_ok=True)

        extract_region_to_fasta(
            genome_fasta=ref_files.genome_fasta,
            region=RegionSpec(
                contig=gene.contig,
                start_1based=gene.start_1based,
                end_1based=gene.end_1based,
            ),
            out_fasta=query_fa,
        )

        hit = run_blastn_top1(query_fna=query_fa, db_prefix=db_prefix)

        if hit is None:
            lines.append(f"[red]{label:<20} -> no hit[/red]")
            continue

        hit_count += 1
        lines.append(
            f"[green]{label:<20} -> {hit.sseqid}:{hit.start_1based}-{hit.end_1based} ({hit.strand}) "
            f"pident={hit.pident:.1f} len={hit.length} e={hit.evalue:.2e}[/green]"
        )

    lines.append(f"\n[bold]Hits found:[/bold] {hit_count}/{len(ref_neighbors)}")
    lines.append(f"[bold]Output folder:[/bold] {out_dir}")

    console.print(
        Panel.fit(
            "[bold green]Neighbor mapping by BLASTN[/bold green]\n\n" + "\n".join(lines)
        )
    )

    if not ref_files.gff3 or not ref_files.genome_fasta:
        console.print(
            "[bold red]Reference must have both GFF3 and genome FASTA.[/bold red]"
        )
        raise typer.Exit(code=1)

    ref_found = find_gene_by_name_fast(ref_files.gff3, gene_name=gene_name)
    if not ref_found:
        console.print(
            Panel.fit(
                f"[bold red]Reference focal gene not found[/bold red]\n\n"
                f"Assembly: {ref_assembly}\nGene: {gene_name}"
            )
        )
        raise typer.Exit(code=1)

    ref_region_start = max(1, ref_found.start_1based - flank)
    ref_region_end = ref_found.end_1based + flank

    ref_genes = genes_in_region_fast(
        gff_path=ref_files.gff3,
        contig=ref_found.contig,
        start_1based=ref_region_start,
        end_1based=ref_region_end,
    )

    # Identify focal gene record inside the region list (same logic as compare_synteny)
    ref_focal = None
    for g in ref_genes:
        if ref_found.gene_id and g.gene_id == ref_found.gene_id:
            ref_focal = g
            break
    if ref_focal is None:
        for g in ref_genes:
            if g.gene_name and g.gene_name.lower() == gene_name.lower():
                ref_focal = g
                break
    if ref_focal is None:
        console.print(
            "[bold red]Internal error: reference focal gene not found inside region list.[/bold red]"
        )
        raise typer.Exit(code=1)

    ref_block = neighbors_around_focal(ref_genes, ref_focal, n_each_side=n)

    # Exclude the focal gene from neighbor mapping.
    ref_neighbors = list(ref_block.upstream) + list(ref_block.downstream)

    if not ref_neighbors:
        console.print("[bold red]No neighbors found in reference block.[/bold red]")
        raise typer.Exit(code=1)

    # Load target genome + GFF

    tgt_unzip = paths.genomes_dir / target_assembly / "unzipped"
    tgt_files = find_genome_files(tgt_unzip)

    if not tgt_files.genome_fasta:
        console.print("[bold red]Target must have genome FASTA.[/bold red]")
        raise typer.Exit(code=1)

    if not tgt_files.gff3:
        console.print(
            "[bold red]Target has no GFF3. Locate its region with map-by-blast first.[/bold red]"
        )
        raise typer.Exit(code=1)

    tgt_found = find_gene_by_name_fast(tgt_files.gff3, gene_name=gene_name)
    if not tgt_found:
        console.print(
            Panel.fit(
                f"[bold red]Target focal gene not found by name[/bold red]\n\n"
                f"Target assembly: {target_assembly}\nGene: {gene_name}\n\n"
                "Use map-by-blast to locate the target region before mapping neighbors."
            )
        )
        raise typer.Exit(code=1)

    tgt_region_start = max(1, tgt_found.start_1based - flank)
    tgt_region_end = tgt_found.end_1based + flank

    # Create run dir + extract target region fasta

    run_dir = paths.new_run_dir()
    out_dir = (
        run_dir
        / "neighbor_mapping"
        / f"{ref_assembly}_to_{target_assembly}"
        / gene_name
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    target_region_fa = (
        out_dir
        / f"target_region_{target_assembly}_{tgt_found.contig}_{tgt_region_start}_{tgt_region_end}.fna"
    )
    extract_region_to_fasta(
        genome_fasta=tgt_files.genome_fasta,
        region=RegionSpec(
            contig=tgt_found.contig,
            start_1based=tgt_region_start,
            end_1based=tgt_region_end,
        ),
        out_fasta=target_region_fa,
    )

    # Build BLAST DB on target region

    db_dir = out_dir / "blastdb"
    db_prefix = db_dir / "target_region_db"

    try:
        if not (db_dir.exists() and any(db_dir.glob("target_region_db.*"))):
            make_nucl_db(target_region_fa, db_prefix)
    except (BlastNotFound, BlastnNotFound) as e:
        console.print(f"[bold red]{e}[/bold red]")
        raise typer.Exit(code=1)

    # Extract each reference neighbor DNA and BLASTN into target region

    results_lines = []
    results_lines.append(f"[bold]Reference:[/bold] {ref_assembly}  (focal={gene_name})")
    results_lines.append(f"[bold]Target:[/bold] {target_assembly}")
    results_lines.append(
        f"[bold]Target region:[/bold] {tgt_found.contig}:{tgt_region_start}-{tgt_region_end} (±{flank})\n"
    )

    hits_found = 0

    for g in ref_neighbors:
        q_fa = (
            out_dir
            / "queries"
            / f"{g.gene_name or g.gene_id}_{g.contig}_{g.start_1based}_{g.end_1based}.fna"
        )
        q_fa.parent.mkdir(parents=True, exist_ok=True)

        extract_region_to_fasta(
            genome_fasta=ref_files.genome_fasta,
            region=RegionSpec(
                contig=g.contig, start_1based=g.start_1based, end_1based=g.end_1based
            ),
            out_fasta=q_fa,
        )

        hit = run_blastn_top1(query_fna=q_fa, db_prefix=db_prefix)

        label = g.gene_name or g.gene_id or "unknown"
        if hit is None:
            results_lines.append(f"[red]{label:<20} -> no hit[/red]")
        else:
            hits_found += 1
            results_lines.append(
                f"[green]{label:<20} -> {hit.sseqid}:{hit.start_1based}-{hit.end_1based} ({hit.strand}) "
                f"pident={hit.pident:.1f} len={hit.length} e={hit.evalue:.2e}[/green]"
            )

    results_lines.append(
        f"\n[bold]Hits found:[/bold] {hits_found}/{len(ref_neighbors)}"
    )
    results_lines.append(f"[bold]Output folder:[/bold] {out_dir}")

    console.print(
        Panel.fit(
            "[bold green]Neighbor mapping by BLASTN[/bold green]\n\n"
            + "\n".join(results_lines)
        )
    )


def resolve_species_focal_region(
    *,
    paths: ProjectPaths,
    assembly: str,
    species_files,
    gene_name: str,
    flank: int,
    focal_query_candidates: list[tuple[str, Path]],
    use_canonical_focal_cache: bool,
):
    # Find the focal region in one species, optionally reusing a cached region

    genome_db_dir = paths.genomes_dir / assembly / "blastdb"
    genome_db = genome_db_dir / "genome_nucl"

    if not (genome_db_dir.exists() and any(genome_db_dir.glob("genome_nucl.*"))):
        make_nucl_db(species_files.genome_fasta, genome_db)

    class FocalHit:
        def __init__(self, sseqid, start_1based, end_1based, strand):
            self.sseqid = sseqid
            self.start_1based = start_1based
            self.end_1based = end_1based
            self.strand = strand

    if use_canonical_focal_cache:
        for source_label, _query_faa in focal_query_candidates:
            query_label = _query_label_for_cache(source_label)
            cache_dir = get_canonical_focal_cache_dir(
                paths=paths,
                target_assembly=assembly,
                gene_name=gene_name,
                query_label=query_label,
                flank=flank,
            )
            cached = load_cached_focal_hit(cache_dir)
            region_db = cache_dir / "blastdb" / "target_region_db"

            if cached is None:
                continue

            cached_region_fa = cache_dir / cached["region_fna_name"]
            has_db = any(region_db.parent.glob(region_db.name + ".*"))

            if cached_region_fa.exists() and has_db:
                hit = FocalHit(
                    sseqid=cached["hit"]["sseqid"],
                    start_1based=int(cached["hit"]["start_1based"]),
                    end_1based=int(cached["hit"]["end_1based"]),
                    strand=cached["hit"]["strand"],
                )
                region_start = int(cached["region"]["start_1based"])
                region_end = int(cached["region"]["end_1based"])

                return (
                    hit,
                    region_start,
                    region_end,
                    cached_region_fa,
                    region_db,
                    source_label,
                    True,
                )

    hit = None
    used_label = None

    for source_label, query_faa in focal_query_candidates:
        try:
            hit = run_tblastn_top1(query_faa=query_faa, db_prefix=genome_db)
        except Exception:
            hit = None

        if hit is not None:
            used_label = source_label
            break

    if hit is None:
        raise RuntimeError(
            f"No focal hit found in {assembly} for any supplied focal query."
        )

    genome_fa = Fasta(
        str(species_files.genome_fasta), as_raw=True, sequence_always_upper=True
    )
    contig_len = len(genome_fa[hit.sseqid])

    region_start = max(1, hit.start_1based - flank)
    region_end = min(hit.end_1based + flank, contig_len)

    if use_canonical_focal_cache:
        query_label = _query_label_for_cache(used_label)
        cache_dir = get_canonical_focal_cache_dir(
            paths=paths,
            target_assembly=assembly,
            gene_name=gene_name,
            query_label=query_label,
            flank=flank,
        )
        cache_dir.mkdir(parents=True, exist_ok=True)

        region_fa = (
            cache_dir
            / f"target_region_{assembly}_{hit.sseqid}_{region_start}_{region_end}.fna"
        )
        if not region_fa.exists():
            extract_region_to_fasta(
                genome_fasta=species_files.genome_fasta,
                region=RegionSpec(
                    contig=hit.sseqid,
                    start_1based=region_start,
                    end_1based=region_end,
                ),
                out_fasta=region_fa,
            )

        region_db = cache_dir / "blastdb" / "target_region_db"
        ensure_region_db(region_fa, region_db)

        save_cached_focal_hit(
            cache_dir,
            {
                "assembly": assembly,
                "gene_name": gene_name,
                "source_label": used_label,
                "flank": flank,
                "hit": {
                    "sseqid": hit.sseqid,
                    "start_1based": hit.start_1based,
                    "end_1based": hit.end_1based,
                    "strand": hit.strand,
                },
                "region": {
                    "start_1based": region_start,
                    "end_1based": region_end,
                },
                "region_fna_name": region_fa.name,
            },
        )

        return hit, region_start, region_end, region_fa, region_db, used_label, False

    return hit, region_start, region_end, None, None, used_label, False


@app.command("map-neighbors-by-tblastn")
def map_neighbors_by_tblastn(
    ref_assembly: str = typer.Option(..., "--ref-assembly", help="Reference assembly"),
    target_assembly: str = typer.Option(
        ..., "--target-assembly", help="Target assembly"
    ),
    gene_name: str = typer.Option(..., "--gene-name", help="Focal gene"),
    flank: int = typer.Option(50000, "--flank", help="Flank size around focal region"),
    n: int = typer.Option(10, "--n", help="Neighbors on each side"),
    root: Path = typer.Option(default_project_root(), "--root"),
    run_dir_override: Path | None = None,
    reciprocal: bool = typer.Option(
        True, "--reciprocal/--no-reciprocal", help="Run reciprocal BLAST checks"
    ),
    reciprocal_trigger_pident: float = typer.Option(
        45.0,
        "--reciprocal-trigger-pident",
        help="Only check hits below this percent identity",
    ),
    focal_accession: str | None = typer.Option(
        None,
        "--focal-accession",
        help="Protein accession used to anchor the focal locus",
    ),
    focal_faa: Path | None = typer.Option(
        None,
        "--focal-faa",
        help="Protein FASTA used to anchor the focal locus",
    ),
    focal_fallback_accession: list[str] = typer.Option(
        [],
        "--focal-fallback-accession",
        help="Optional fallback focal accessions",
    ),
    use_canonical_focal_cache: bool = typer.Option(
        False,
        "--use-canonical-focal-cache/--no-use-canonical-focal-cache",
        help="Reuse one cached focal region per target species",
    ),
):
    # Map reference neighbors into a target region using BLASTP/TBLASTN/BLASTN

    paths = ProjectPaths(root=root)
    paths.ensure()

    if focal_faa is not None and focal_accession is not None:
        console.print("[red]Use only one of --focal-accession or --focal-faa[/red]")
        raise typer.Exit(code=1)

    if focal_faa is not None and not focal_faa.exists():
        console.print(f"[red]Missing focal FASTA: {focal_faa}[/red]")
        raise typer.Exit(code=1)

    if use_canonical_focal_cache and focal_faa is None and focal_accession is None:
        console.print(
            "[red]Canonical focal cache needs --focal-accession or --focal-faa[/red]"
        )
        raise typer.Exit(code=1)

    ref_unzip = paths.genomes_dir / ref_assembly / "unzipped"
    ref_files = find_genome_files(ref_unzip)

    if not ref_files.gff3:
        console.print("[red]Reference assembly has no GFF3[/red]")
        raise typer.Exit(code=1)

    try:
        focal_queries = resolve_focal_query_candidates(
            paths=paths,
            ref_files=ref_files,
            gene_name=gene_name,
            focal_accession=focal_accession,
            focal_faa=focal_faa,
            focal_fallback_accessions=focal_fallback_accession,
        )
    except (FastaFetchError, RuntimeError) as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)

    try:
        block, ref_focal, _ref_start, _ref_end, focal_mode = resolve_reference_block(
            paths=paths,
            ref_assembly=ref_assembly,
            ref_files=ref_files,
            gene_name=gene_name,
            flank=flank,
            n=n,
            focal_query_candidates=focal_queries,
            use_canonical_focal_cache=use_canonical_focal_cache,
        )
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)

    neighbors = [
        gene.gene_name
        for gene in (block.upstream + block.downstream)
        if gene.gene_name and keep_neighbor_gene(gene.gene_name)
    ]

    if not neighbors:
        console.print("[red]No usable neighbors found in reference block[/red]")
        raise typer.Exit(code=1)

    tgt_unzip = paths.genomes_dir / target_assembly / "unzipped"
    tgt_files = find_genome_files(tgt_unzip)

    if not tgt_files.genome_fasta:
        console.print("[red]Target assembly has no genome FASTA[/red]")
        raise typer.Exit(code=1)

    target_prot_db = None
    if tgt_files.protein_faa:
        prot_db_dir = paths.genomes_dir / target_assembly / "blastdb"
        target_prot_db = prot_db_dir / "proteome_prot"
        try:
            if not any(prot_db_dir.glob("proteome_prot.*")):
                make_prot_db(tgt_files.protein_faa, target_prot_db)
        except (BlastpNotFound, RuntimeError) as e:
            console.print(
                f"[yellow]Could not build target proteome DB, falling back to genome only: {e}[/yellow]"
            )
            target_prot_db = None

    run_dir = run_dir_override if run_dir_override is not None else paths.new_run_dir()
    out_dir = (
        run_dir
        / "neighbor_mapping"
        / f"{ref_assembly}_to_{target_assembly}"
        / gene_name
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        (
            focal_hit,
            tgt_start,
            tgt_end,
            cached_region_fa,
            cached_region_db,
            used_query,
            used_cache,
        ) = resolve_species_focal_region(
            paths=paths,
            assembly=target_assembly,
            species_files=tgt_files,
            gene_name=gene_name,
            flank=flank,
            focal_query_candidates=focal_queries,
            use_canonical_focal_cache=use_canonical_focal_cache,
        )
    except BlastNotFound as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)

    if use_canonical_focal_cache:
        region_fa = cached_region_fa
        region_db = cached_region_db
        region_db_dir = region_db.parent

        local_region_fa = out_dir / region_fa.name
        if not local_region_fa.exists():
            shutil.copy2(region_fa, local_region_fa)
    else:
        region_fa = (
            out_dir
            / f"target_region_{target_assembly}_{focal_hit.sseqid}_{tgt_start}_{tgt_end}.fna"
        )

        extract_region_to_fasta(
            genome_fasta=tgt_files.genome_fasta,
            region=RegionSpec(
                contig=focal_hit.sseqid,
                start_1based=tgt_start,
                end_1based=tgt_end,
            ),
            out_fasta=region_fa,
        )

        region_db_dir = out_dir / "blastdb"
        region_db_dir.mkdir(parents=True, exist_ok=True)
        region_db = region_db_dir / "target_region_db"

    try:
        if not (
            region_db_dir.exists() and any(region_db_dir.glob("target_region_db.*"))
        ):
            make_nucl_db(region_fa, region_db)
    except BlastNotFound as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=1)

    hits: list[NeighborHit] = []
    lines = [
        f"[bold]Reference:[/bold] {ref_assembly} (focal={gene_name}, mode={focal_mode})",
        f"[bold]Target:[/bold] {target_assembly}",
        f"[bold]Target region:[/bold] {focal_hit.sseqid}:{tgt_start}-{tgt_end} (±{flank})",
    ]

    if use_canonical_focal_cache:
        cache_status = "reused" if used_cache else "built"
        lines.append(f"[bold]Canonical cache:[/bold] {cache_status} ({used_query})")

    lines.append("")
    hit_count = 0

    for symbol in neighbors:
        protein_id = protein_id_for_gene_symbol_first(ref_files.gff3, symbol)

        if not protein_id:
            ref_gene = find_gene_by_name_fast(ref_files.gff3, gene_name=symbol)
            if not ref_gene:
                lines.append(
                    f"{symbol:20s} -> [yellow]missing in reference GFF[/yellow]"
                )
                hits.append(
                    NeighborHit(
                        neighbor_name=symbol,
                        method="blastn",
                        hit_found=False,
                        query_protein_length=None,
                    )
                )
                continue

            if not ref_files.genome_fasta:
                lines.append(
                    f"{symbol:20s} -> [yellow]reference genome FASTA missing[/yellow]"
                )
                hits.append(
                    NeighborHit(
                        neighbor_name=symbol,
                        method="blastn",
                        hit_found=False,
                        query_protein_length=None,
                    )
                )
                continue

            query_fa = (
                out_dir
                / "queries_fallback_dna"
                / f"{symbol}_{ref_gene.contig}_{ref_gene.start_1based}_{ref_gene.end_1based}.fna"
            )
            query_fa.parent.mkdir(parents=True, exist_ok=True)

            extract_region_to_fasta(
                genome_fasta=ref_files.genome_fasta,
                region=RegionSpec(
                    contig=ref_gene.contig,
                    start_1based=ref_gene.start_1based,
                    end_1based=ref_gene.end_1based,
                ),
                out_fasta=query_fa,
            )

            try:
                hit = run_blastn_top1(query_fna=query_fa, db_prefix=region_db)
            except Exception as e:
                lines.append(f"{symbol:20s} -> [yellow]blastn failed: {e}[/yellow]")
                hits.append(
                    NeighborHit(
                        neighbor_name=symbol,
                        method="blastn",
                        hit_found=False,
                        query_protein_length=None,
                    )
                )
                continue

            if hit is None:
                lines.append(f"{symbol:20s} -> no hit [blastn]")
                hits.append(
                    NeighborHit(
                        neighbor_name=symbol,
                        method="blastn",
                        hit_found=False,
                        query_protein_length=None,
                    )
                )
                continue

            hit_count += 1
            lines.append(
                f"{symbol:20s} -> {hit.sseqid}:{hit.start_1based}-{hit.end_1based} ({hit.strand}) "
                f"pident={hit.pident:.1f} len={hit.length} e={hit.evalue:.2g} [blastn]"
            )
            hits.append(
                NeighborHit(
                    neighbor_name=symbol,
                    method="blastn",
                    hit_found=True,
                    sseqid=hit.sseqid,
                    start_1based=hit.start_1based,
                    end_1based=hit.end_1based,
                    strand=hit.strand,
                    pident=hit.pident,
                    length=hit.length,
                    evalue=hit.evalue,
                    bitscore=hit.bitscore,
                    query_protein_length=None,
                )
            )
            continue

        try:
            query_faa = fetch_protein_fasta(
                accession=protein_id,
                out_fasta=paths.data_dir / "queries" / f"{protein_id}.faa",
            ).fasta_path
            query_len = fasta_sequence_length(query_faa)
        except FastaFetchError as e:
            lines.append(f"{symbol:20s} -> [yellow]protein fetch failed: {e}[/yellow]")
            hits.append(
                NeighborHit(
                    neighbor_name=symbol,
                    method="tblastn",
                    hit_found=False,
                    query_protein_length=None,
                )
            )
            continue

        hit = None
        method = None

        if target_prot_db is not None:
            try:
                prot_hit = run_blastp_top1(
                    query_faa=query_faa, db_prefix=target_prot_db
                )
            except Exception as e:
                prot_hit = None
                lines.append(
                    f"{symbol:20s} -> [yellow]blastp failed, using tblastn: {e}[/yellow]"
                )

            if prot_hit is not None and tgt_files.gff3:
                target_gene = find_gene_by_name_fast(tgt_files.gff3, gene_name=symbol)
                if target_gene:

                    class SimpleHit:
                        def __init__(
                            self,
                            sseqid,
                            start_1based,
                            end_1based,
                            strand,
                            pident,
                            length,
                            evalue,
                            bitscore,
                        ):
                            self.sseqid = sseqid
                            self.start_1based = start_1based
                            self.end_1based = end_1based
                            self.strand = strand
                            self.pident = pident
                            self.length = length
                            self.evalue = evalue
                            self.bitscore = bitscore

                    hit = SimpleHit(
                        sseqid=target_gene.contig,
                        start_1based=target_gene.start_1based,
                        end_1based=target_gene.end_1based,
                        strand=target_gene.strand,
                        pident=prot_hit.pident,
                        length=prot_hit.length,
                        evalue=prot_hit.evalue,
                        bitscore=prot_hit.bitscore,
                    )
                    method = "blastp"

        if hit is None:
            hit = run_tblastn_top1(query_faa=query_faa, db_prefix=region_db)
            method = "tblastn"

        if hit is None:
            lines.append(f"{symbol:20s} -> no hit")
            hits.append(
                NeighborHit(
                    neighbor_name=symbol,
                    method=method or "tblastn",
                    hit_found=False,
                    query_protein_length=query_len,
                )
            )
            continue

        hit_count += 1
        lines.append(
            f"{symbol:20s} -> {hit.sseqid}:{hit.start_1based}-{hit.end_1based} ({hit.strand}) "
            f"pident={hit.pident:.1f} len={hit.length} e={hit.evalue:.2g} [{method}]"
        )
        hits.append(
            NeighborHit(
                neighbor_name=symbol,
                method=method or "tblastn",
                hit_found=True,
                sseqid=hit.sseqid,
                start_1based=hit.start_1based,
                end_1based=hit.end_1based,
                strand=hit.strand,
                pident=hit.pident,
                length=hit.length,
                evalue=hit.evalue,
                bitscore=hit.bitscore,
                query_protein_length=query_len,
            )
        )

    reciprocal_rows = []

    if reciprocal:
        reciprocal_dir = out_dir / "reciprocal"
        reciprocal_dir.mkdir(parents=True, exist_ok=True)

        ref_prot_fa = reciprocal_dir / "ref_neighbors_proteins.faa"
        protein_records = []

        for symbol in neighbors:
            protein_id = protein_id_for_gene_symbol_first(ref_files.gff3, symbol)
            if not protein_id:
                continue

            protein_fa = paths.data_dir / "queries" / f"{protein_id}.faa"
            if protein_fa.exists():
                header, seq = fasta_read_one(protein_fa)
                protein_records.append((protein_id, seq))

        if protein_records:
            with ref_prot_fa.open("w", encoding="utf-8") as handle:
                for accession, seq in protein_records:
                    handle.write(f">{accession}\n")
                    handle.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

            prot_db = reciprocal_dir / "ref_neighbors_prot_db"
            if not any(prot_db.parent.glob(prot_db.name + ".*")):
                make_prot_db(ref_prot_fa, prot_db)
        else:
            prot_db = None

        ref_dna_fa = reciprocal_dir / "ref_neighbors_dna.fna"
        fallback_dir = out_dir / "queries_fallback_dna"

        if fallback_dir.exists():
            with ref_dna_fa.open("w", encoding="utf-8") as handle:
                for file in fallback_dir.glob("*.fna"):
                    header, seq = fasta_read_one(file)
                    symbol = file.name.split("_")[0]
                    handle.write(f">{symbol}\n")
                    handle.write("\n".join(textwrap.wrap(seq, 80)) + "\n")

            dna_db = reciprocal_dir / "ref_neighbors_nucl_db"
            if not any(dna_db.parent.glob(dna_db.name + ".*")):
                make_nucl_db(ref_dna_fa, dna_db)
        else:
            dna_db = None

        for hit in hits:
            if not hit.hit_found:
                continue
            if hit.pident is None or hit.pident >= reciprocal_trigger_pident:
                continue

            try:
                target_seq = extract_subseq_from_single_fasta(
                    region_fa,
                    hit.start_1based,
                    hit.end_1based,
                    hit.strand,
                )
            except Exception as e:
                reciprocal_rows.append(
                    {
                        "neighbor": hit.neighbor_name,
                        "method": hit.method,
                        "reciprocal_pass": False,
                        "reason": f"extract_failed: {e}",
                    }
                )
                continue

            back_query = reciprocal_dir / f"{hit.neighbor_name}_targetHit.fna"
            fasta_write_one(
                back_query, header=f"{hit.neighbor_name}_targetHit", seq=target_seq
            )

            if hit.method == "tblastn":
                if prot_db is None:
                    reciprocal_rows.append(
                        {
                            "neighbor": hit.neighbor_name,
                            "method": "blastx",
                            "reciprocal_pass": False,
                            "reason": "no_ref_protein_db",
                        }
                    )
                    continue

                back_hit = run_blastx_top1(query_fna=back_query, db_prefix=prot_db)
                if back_hit is None:
                    reciprocal_rows.append(
                        {
                            "neighbor": hit.neighbor_name,
                            "method": "blastx",
                            "reciprocal_pass": False,
                            "reason": "no_back_hit",
                        }
                    )
                    continue

                expected = protein_id_for_gene_symbol_first(
                    ref_files.gff3, hit.neighbor_name
                )

                reciprocal_rows.append(
                    {
                        "neighbor": hit.neighbor_name,
                        "method": "blastx",
                        "back_best": back_hit["sseqid"],
                        "expected": expected,
                        "back_pident": back_hit["pident"],
                        "back_len": back_hit["length"],
                        "back_evalue": back_hit["evalue"],
                        "reciprocal_pass": expected is not None
                        and back_hit["sseqid"] == expected,
                    }
                )

            elif hit.method == "blastn":
                if dna_db is None:
                    reciprocal_rows.append(
                        {
                            "neighbor": hit.neighbor_name,
                            "method": "blastn_back",
                            "reciprocal_pass": False,
                            "reason": "no_ref_nucl_db",
                        }
                    )
                    continue

                back_hit = run_blastn_top1(query_fna=back_query, db_prefix=dna_db)
                if back_hit is None:
                    reciprocal_rows.append(
                        {
                            "neighbor": hit.neighbor_name,
                            "method": "blastn_back",
                            "reciprocal_pass": False,
                            "reason": "no_back_hit",
                        }
                    )
                    continue

                reciprocal_rows.append(
                    {
                        "neighbor": hit.neighbor_name,
                        "method": "blastn_back",
                        "back_best": back_hit.sseqid,
                        "expected": hit.neighbor_name,
                        "back_pident": back_hit.pident,
                        "back_len": back_hit.length,
                        "back_evalue": back_hit.evalue,
                        "reciprocal_pass": back_hit.sseqid == hit.neighbor_name,
                    }
                )

        if reciprocal_rows:
            rec_json = reciprocal_dir / "reciprocal.json"
            rec_tsv = reciprocal_dir / "reciprocal.tsv"

            rec_json.write_text(json.dumps(reciprocal_rows, indent=2), encoding="utf-8")

            header = [
                "neighbor",
                "method",
                "reciprocal_pass",
                "back_best",
                "expected",
                "back_pident",
                "back_len",
                "back_evalue",
                "reason",
            ]
            tsv_lines = ["\t".join(header) + "\n"]
            for row in reciprocal_rows:
                tsv_lines.append(
                    "\t".join(str(row.get(col, "")) for col in header) + "\n"
                )
            rec_tsv.write_text("".join(tsv_lines), encoding="utf-8")

            lines.append(f"\n[bold]Reciprocal checks:[/bold] {rec_json} | {rec_tsv}")

    mapping = NeighborMappingResult(
        ref_assembly=ref_assembly,
        target_assembly=target_assembly,
        gene_name=gene_name,
        flank=flank,
        n=n,
        target_contig=focal_hit.sseqid,
        target_region_start_1based=tgt_start,
        target_region_end_1based=tgt_end,
        hits=hits,
    )
    json_path = mapping.to_json(out_dir / "mapping.json")

    lines.append(f"\n[bold]Hits found:[/bold] {hit_count}/{len(neighbors)}")
    lines.append(f"[bold]Output folder:[/bold] {out_dir}")
    lines.append(f"[bold]JSON:[/bold] {json_path}")

    console.print(
        Panel.fit(
            "[bold green]Neighbor mapping by TBLASTN[/bold green]\n\n"
            + "\n".join(lines)
        )
    )


@app.command("map-and-score")
def map_and_score(
    ref_assembly: str = typer.Option(..., "--ref-assembly"),
    target_assembly: str = typer.Option(..., "--target-assembly"),
    gene_name: str = typer.Option(..., "--gene-name"),
    flank: int = typer.Option(50000, "--flank"),
    n: int = typer.Option(10, "--n"),
    min_pident: float = typer.Option(30.0, "--min-pident"),
    max_evalue: float = typer.Option(1e-5, "--max-evalue"),
    min_len: int = typer.Option(60, "--min-len"),
    root: Path = typer.Option(default_project_root(), "--root"),
    run_dir_override: Path | None = None,
    focal_accession: str | None = typer.Option(None, "--focal-accession"),
    focal_faa: Path | None = typer.Option(None, "--focal-faa", exists=True),
    focal_fallback_accession: list[str] = typer.Option(
        [], "--focal-fallback-accession"
    ),
    reciprocal: bool = typer.Option(True, "--reciprocal/--no-reciprocal"),
    use_canonical_focal_cache: bool = typer.Option(
        False,
        "--use-canonical-focal-cache/--no-use-canonical-focal-cache",
    ),
):
    # Run mapping and write a simple score summary

    paths = ProjectPaths(root=root)
    paths.ensure()

    run_dir = run_dir_override if run_dir_override is not None else paths.new_run_dir()
    out_dir = (
        run_dir
        / "neighbor_mapping"
        / f"{ref_assembly}_to_{target_assembly}"
        / gene_name
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    map_neighbors_by_tblastn(
        ref_assembly=ref_assembly,
        target_assembly=target_assembly,
        gene_name=gene_name,
        flank=flank,
        n=n,
        root=root,
        run_dir_override=run_dir,
        reciprocal=reciprocal,
        reciprocal_trigger_pident=45.0,
        focal_accession=focal_accession,
        focal_faa=focal_faa,
        focal_fallback_accession=focal_fallback_accession,
        use_canonical_focal_cache=use_canonical_focal_cache,
    )

    mapping_path = out_dir / "mapping.json"
    if not mapping_path.exists():
        console.print(f"[red]Missing mapping.json: {mapping_path}[/red]")
        raise typer.Exit(code=1)

    cfg = ScoreConfig(
        min_pident=min_pident,
        max_evalue=max_evalue,
        min_align_len=min_len,
    )
    score = score_mapping_json(mapping_path, cfg)

    score_json = out_dir / "score.json"
    score_tsv = out_dir / "score.tsv"
    score.to_json(score_json)
    score.to_tsv(score_tsv)

    console.print(
        Panel.fit(
            "[bold green]Mapping + scoring done[/bold green]\n\n"
            f"[bold]Run folder:[/bold] {run_dir}\n"
            f"[bold]Mapping:[/bold] {mapping_path}\n"
            f"[bold]Score JSON:[/bold] {score_json}\n"
            f"[bold]Score TSV:[/bold] {score_tsv}\n"
            f"[bold]Thresholds:[/bold] min_pident={min_pident}, max_evalue={max_evalue}, min_len={min_len}\n"
            f"[bold]Score:[/bold] {score.score_fraction_passing:.3f} "
            f"({score.hits_passing}/{score.total_neighbors} passing)"
        )
    )


@app.command("ref-vs-many")
def ref_vs_many(
    ref_assembly: str = typer.Option(..., "--ref-assembly"),
    targets: str = typer.Option(
        ..., "--targets", help="Comma-separated target assemblies"
    ),
    gene_name: str = typer.Option(..., "--gene-name"),
    flank: int = typer.Option(50000, "--flank"),
    n: int = typer.Option(10, "--n"),
    min_pident: float = typer.Option(30.0, "--min-pident"),
    max_evalue: float = typer.Option(1e-5, "--max-evalue"),
    min_len: int = typer.Option(60, "--min-len"),
    root: Path = typer.Option(default_project_root(), "--root"),
    focal_accession: str | None = typer.Option(None, "--focal-accession"),
    focal_faa: Path | None = typer.Option(None, "--focal-faa", exists=True),
    focal_fallback_accession: list[str] = typer.Option(
        [], "--focal-fallback-accession"
    ),
    reciprocal: bool = typer.Option(True, "--reciprocal/--no-reciprocal"),
    use_canonical_focal_cache: bool = typer.Option(
        False,
        "--use-canonical-focal-cache/--no-use-canonical-focal-cache",
    ),
):
    """Run one reference against a set of targets."""

    target_list = [target.strip() for target in targets.split(",") if target.strip()]
    if not target_list:
        console.print("[red]No targets provided[/red]")
        raise typer.Exit(code=1)

    paths = ProjectPaths(root=root)
    paths.ensure()

    run_dir = paths.new_run_dir()
    out_tsv = run_dir / "scores.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "ref_assembly\ttarget_assembly\tgene_name\tscore_fraction_passing\thits_found\thits_passing\ttotal_neighbors\tmapping_json\n"
    ]

    for target in target_list:
        try:
            map_and_score(
                ref_assembly=ref_assembly,
                target_assembly=target,
                gene_name=gene_name,
                flank=flank,
                n=n,
                min_pident=min_pident,
                max_evalue=max_evalue,
                min_len=min_len,
                root=root,
                run_dir_override=run_dir,
                focal_accession=focal_accession,
                focal_faa=focal_faa,
                focal_fallback_accession=focal_fallback_accession,
                reciprocal=reciprocal,
                use_canonical_focal_cache=use_canonical_focal_cache,
            )
        except Exception as e:
            console.print(f"[red]{target} failed:[/red] {e}")
            lines.append(
                f"{ref_assembly}\t{target}\t{gene_name}\t\t\t\t\tFAILED: {e}\n"
            )
            continue

        mapping_json = (
            run_dir
            / "neighbor_mapping"
            / f"{ref_assembly}_to_{target}"
            / gene_name
            / "mapping.json"
        )
        score_json = mapping_json.parent / "score.json"

        if not score_json.exists():
            lines.append(
                f"{ref_assembly}\t{target}\t{gene_name}\t\t\t\t\tMISSING_SCORE\n"
            )
            continue

        score_data = json.loads(score_json.read_text(encoding="utf-8"))
        lines.append(
            f"{ref_assembly}\t{target}\t{gene_name}\t"
            f"{score_data.get('score_fraction_passing', '')}\t"
            f"{score_data.get('hits_found', '')}\t"
            f"{score_data.get('hits_passing', '')}\t"
            f"{score_data.get('total_neighbors', '')}\t"
            f"{mapping_json}\n"
        )

    out_tsv.write_text("".join(lines), encoding="utf-8")

    console.print(
        Panel.fit(
            "[bold green]Ref vs many done[/bold green]\n\n"
            f"[bold]Output:[/bold] {out_tsv}\n"
            f"[bold]Targets:[/bold] {len(target_list)}"
        )
    )


def find_existing_mapping(
    run_dirs: list[Path],
    ref: str,
    tgt: str,
    gene_name: str,
) -> Path | None:
    for run_dir in run_dirs:
        mapping_path = (
            run_dir
            / "neighbor_mapping"
            / f"{ref}_to_{tgt}"
            / gene_name
            / "mapping.json"
        )
        if mapping_path.exists():
            return mapping_path
    return None


@app.command("matrix")
def matrix(
    assemblies: str = typer.Option(
        ..., "--assemblies", help="Comma-separated assemblies"
    ),
    gene_name: str = typer.Option(..., "--gene-name"),
    flank: int = typer.Option(50000, "--flank"),
    n: int = typer.Option(10, "--n"),
    min_pident: float = typer.Option(30.0, "--min-pident"),
    max_evalue: float = typer.Option(1e-5, "--max-evalue"),
    min_len: int = typer.Option(60, "--min-len"),
    root: Path = typer.Option(default_project_root(), "--root"),
    focal_accession: str | None = typer.Option(None, "--focal-accession"),
    focal_faa: Path | None = typer.Option(None, "--focal-faa", exists=True),
    focal_fallback_accession: list[str] = typer.Option(
        [], "--focal-fallback-accession"
    ),
    reciprocal: bool = typer.Option(True, "--reciprocal/--no-reciprocal"),
    reuse_run_folder: list[Path] = typer.Option(
        [],
        "--reuse-run-folder",
        help="Previous run folders to reuse if mappings already exist",
    ),
    use_canonical_focal_cache: bool = typer.Option(
        False,
        "--use-canonical-focal-cache/--no-use-canonical-focal-cache",
        help="Reuse one cached focal region per target species",
    ),
):
    """Build an all-vs-all score matrix from pairwise mappings."""

    paths = ProjectPaths(root=root)
    paths.ensure()

    if use_canonical_focal_cache and focal_faa is None and focal_accession is None:
        console.print(
            "[red]Canonical focal cache needs --focal-accession or --focal-faa[/red]"
        )
        raise typer.Exit(code=1)

    assembly_list = [
        assembly.strip() for assembly in assemblies.split(",") if assembly.strip()
    ]
    if len(assembly_list) < 2:
        console.print("[red]Need at least two assemblies[/red]")
        raise typer.Exit(code=1)

    run_dir = paths.new_run_dir()
    out_dir = run_dir / "matrix" / gene_name
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = ScoreConfig(
        min_pident=min_pident,
        max_evalue=max_evalue,
        min_align_len=min_len,
    )

    scores: dict[tuple[str, str], float] = {}
    search_dirs = [run_dir] + reuse_run_folder

    for ref in assembly_list:
        for target in assembly_list:
            if ref == target:
                scores[(ref, target)] = 1.0
                continue

            console.print(f"[cyan]{ref} -> {target}[/cyan]")

            new_mapping = (
                run_dir
                / "neighbor_mapping"
                / f"{ref}_to_{target}"
                / gene_name
                / "mapping.json"
            )

            existing_mapping = find_existing_mapping(
                run_dirs=search_dirs,
                ref=ref,
                tgt=target,
                gene_name=gene_name,
            )

            if existing_mapping is not None:
                try:
                    score = score_mapping_json(existing_mapping, cfg)
                    score.to_json(existing_mapping.parent / "score.json")
                    score.to_tsv(existing_mapping.parent / "score.tsv")
                    scores[(ref, target)] = float(score.score_fraction_passing)

                    console.print(
                        f"[green]Reused existing mapping[/green] {ref} -> {target} "
                        f"(score={score.score_fraction_passing:.3f})"
                    )
                    continue
                except Exception as e:
                    console.print(
                        f"[yellow]Could not score existing mapping, rerunning: {e}[/yellow]"
                    )

            try:
                map_neighbors_by_tblastn(
                    ref_assembly=ref,
                    target_assembly=target,
                    gene_name=gene_name,
                    flank=flank,
                    n=n,
                    root=root,
                    run_dir_override=run_dir,
                    reciprocal=reciprocal,
                    reciprocal_trigger_pident=45.0,
                    focal_accession=focal_accession,
                    focal_faa=focal_faa,
                    focal_fallback_accession=focal_fallback_accession,
                    use_canonical_focal_cache=use_canonical_focal_cache,
                )
            except Exception as e:
                console.print(f"[red]{ref} -> {target} failed:[/red] {e}")
                scores[(ref, target)] = 0.0
                continue

            if not new_mapping.exists():
                console.print(
                    f"[yellow]No mapping produced for {ref} -> {target}[/yellow]"
                )
                scores[(ref, target)] = 0.0
                continue

            try:
                score = score_mapping_json(new_mapping, cfg)
                score.to_json(new_mapping.parent / "score.json")
                score.to_tsv(new_mapping.parent / "score.tsv")
                scores[(ref, target)] = float(score.score_fraction_passing)
                console.print(
                    f"[green]Done[/green] {ref} -> {target} "
                    f"(score={score.score_fraction_passing:.3f})"
                )
            except Exception as e:
                console.print(f"[red]Scoring failed for {ref} -> {target}:[/red] {e}")
                scores[(ref, target)] = 0.0

    matrix_tsv = out_dir / "matrix.tsv"
    with matrix_tsv.open("w", encoding="utf-8") as handle:
        handle.write("ref\\target\t" + "\t".join(assembly_list) + "\n")
        for ref in assembly_list:
            row = [ref]
            for target in assembly_list:
                row.append(f"{scores.get((ref, target), 0.0):.3f}")
            handle.write("\t".join(row) + "\n")

    matrix_png = out_dir / "matrix.png"
    matrix_values = [
        [scores.get((ref, target), 0.0) for target in assembly_list]
        for ref in assembly_list
    ]

    plt.figure()
    plt.imshow(matrix_values, aspect="auto", vmin=0.0, vmax=1.0)
    plt.colorbar(label="score_fraction_passing")
    plt.xticks(range(len(assembly_list)), assembly_list, rotation=90)
    plt.yticks(range(len(assembly_list)), assembly_list)
    plt.tight_layout()
    plt.savefig(matrix_png, dpi=200)
    plt.close()

    console.print(
        Panel.fit(
            "[bold green]Matrix complete[/bold green]\n\n"
            f"[bold]Gene:[/bold] {gene_name}\n"
            f"[bold]Assemblies:[/bold] {len(assembly_list)}\n"
            f"[bold]Run folder:[/bold] {run_dir}\n"
            f"[bold]TSV:[/bold] {matrix_tsv}\n"
            f"[bold]PNG:[/bold] {matrix_png}\n"
            f"[bold]Reciprocal:[/bold] {'on' if reciprocal else 'off'}\n"
            f"[bold]Canonical cache:[/bold] {'on' if use_canonical_focal_cache else 'off'}\n"
            f"[bold]Reuse folders:[/bold] {', '.join(str(path) for path in reuse_run_folder) if reuse_run_folder else '(none)'}"
        )
    )


@app.command("reciprocal-blastp")
def reciprocal_blastp(
    ref_assembly: str = typer.Option(..., "--ref-assembly"),
    targets: str = typer.Option(
        ..., "--targets", help="Comma-separated target assemblies"
    ),
    query_accession: str = typer.Option(
        ..., "--query-accession", help="Reference protein accession, e.g. BGN NP_..."
    ),
    query_name: str = typer.Option(
        ..., "--query-name", help="Gene name label, e.g. BGN"
    ),
    evalue: float = typer.Option(1e-5, "--evalue"),
    root: Path = typer.Option(default_project_root(), "--root"),
):
    """
    BLASTP-only reciprocal ortholog check:
    query protein -> target proteome -> best target hit -> reference proteome.
    """
    paths = ProjectPaths(root=root)
    paths.ensure()

    run_dir = paths.new_run_dir()
    out_dir = run_dir / "reciprocal_blastp" / query_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Fetch query protein
    query_faa = fetch_protein_fasta(
        accession=query_accession,
        out_fasta=paths.data_dir / "queries" / f"{query_accession}.faa",
    ).fasta_path

    # Reference proteome DB
    ref_proteome = find_protein_fasta(paths.genomes_dir / ref_assembly / "unzipped")

    ref_db_prefix = paths.genomes_dir / ref_assembly / "blastdb" / "proteome_prot"
    if not any(ref_db_prefix.parent.glob(ref_db_prefix.name + ".*")):
        make_prot_db(ref_proteome, ref_db_prefix)

    target_list = [t.strip() for t in targets.split(",") if t.strip()]

    rows = []
    tsv_lines = [
        "query_name\tquery_accession\tref_assembly\ttarget_assembly\t"
        "target_best_hit\tpident\tqcovs\tlength\tevalue\tbitscore\t"
        "reverse_best_hit\treverse_pident\treverse_qcovs\treverse_evalue\treciprocal_pass\n"
    ]

    for target in target_list:
        try:
            target_unzip = paths.genomes_dir / target / "unzipped"
            target_proteome = find_protein_fasta(target_unzip)

            target_db_prefix = paths.genomes_dir / target / "blastdb" / "proteome_prot"
            if not any(target_db_prefix.parent.glob(target_db_prefix.name + ".*")):
                make_prot_db(target_proteome, target_db_prefix)

            forward = run_blastp_table_top1(query_faa, target_db_prefix, evalue=evalue)

            if forward is None:
                row = {
                    "target_assembly": target,
                    "target_best_hit": "",
                    "status": "NO_FORWARD_HIT",
                }
                rows.append(row)
                tsv_lines.append(
                    f"{query_name}\t{query_accession}\t{ref_assembly}\t{target}\t"
                    f"\t\t\t\t\t\t\t\t\t\t0\n"
                )
                continue

            # Extract target best protein sequence
            rec = fasta_get_record_by_id(target_proteome, forward["sseqid"])
            if rec is None:
                row = {
                    "target_assembly": target,
                    "target_best_hit": forward["sseqid"],
                    "status": "TARGET_SEQUENCE_NOT_FOUND",
                }
                rows.append(row)
                tsv_lines.append(
                    f"{query_name}\t{query_accession}\t{ref_assembly}\t{target}\t"
                    f"{forward['sseqid']}\t{forward['pident']}\t{forward['qcovs']}\t"
                    f"{forward['length']}\t{forward['evalue']}\t{forward['bitscore']}\t"
                    f"\t\t\t\t0\n"
                )
                continue

            target_header, target_seq = rec
            back_query = (
                out_dir
                / "target_best_hits"
                / f"{target}_{normalize_acc(forward['sseqid'])}.faa"
            )
            fasta_write_one(back_query, target_header, target_seq)

            reverse = run_blastp_table_top1(back_query, ref_db_prefix, evalue=evalue)

            reverse_best = reverse["sseqid"] if reverse else ""
            reciprocal_pass = False

            if reverse:
                reciprocal_pass = (
                    normalize_acc(reverse_best) == normalize_acc(query_accession)
                    or normalize_acc(query_accession).split(".")[0]
                    in normalize_acc(reverse_best)
                    or normalize_acc(reverse_best).split(".")[0]
                    == normalize_acc(query_accession).split(".")[0]
                )

            row = {
                "query_name": query_name,
                "query_accession": query_accession,
                "ref_assembly": ref_assembly,
                "target_assembly": target,
                "target_best_hit": forward["sseqid"],
                "pident": forward["pident"],
                "qcovs": forward["qcovs"],
                "length": forward["length"],
                "evalue": forward["evalue"],
                "bitscore": forward["bitscore"],
                "reverse_best_hit": reverse_best,
                "reverse_pident": reverse["pident"] if reverse else None,
                "reverse_qcovs": reverse["qcovs"] if reverse else None,
                "reverse_evalue": reverse["evalue"] if reverse else None,
                "reciprocal_pass": reciprocal_pass,
            }
            rows.append(row)

            reverse_pident_text = "" if reverse is None else f"{reverse['pident']:.2f}"
            reverse_qcovs_text = "" if reverse is None else f"{reverse['qcovs']:.1f}"
            reverse_evalue_text = "" if reverse is None else f"{reverse['evalue']:.2g}"

            tsv_lines.append(
                f"{query_name}\t{query_accession}\t{ref_assembly}\t{target}\t"
                f"{forward['sseqid']}\t{forward['pident']:.2f}\t{forward['qcovs']:.1f}\t"
                f"{forward['length']}\t{forward['evalue']:.2g}\t{forward['bitscore']:.1f}\t"
                f"{reverse_best}\t"
                f"{reverse_pident_text}\t"
                f"{reverse_qcovs_text}\t"
                f"{reverse_evalue_text}\t"
                f"{1 if reciprocal_pass else 0}\n"
            )

        except Exception as e:
            rows.append(
                {
                    "target_assembly": target,
                    "status": "FAILED",
                    "error": str(e),
                }
            )
            tsv_lines.append(
                f"{query_name}\t{query_accession}\t{ref_assembly}\t{target}\t"
                f"FAILED\t\t\t\t\t\t\t\t\t\t0\n"
            )

    out_json = out_dir / "reciprocal_blastp.json"
    out_tsv = out_dir / "reciprocal_blastp.tsv"

    out_json.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    out_tsv.write_text("".join(tsv_lines), encoding="utf-8")

    console.print(
        Panel.fit(
            "[bold green]Reciprocal BLASTP finished[/bold green]\n\n"
            f"[bold]Query:[/bold] {query_name} ({query_accession})\n"
            f"[bold]Reference:[/bold] {ref_assembly}\n"
            f"[bold]Targets:[/bold] {len(target_list)}\n"
            f"[bold]TSV:[/bold] {out_tsv}\n"
            f"[bold]JSON:[/bold] {out_json}"
        )
    )


@app.command("show-mapping")
def show_mapping(
    mapping_json: Path = typer.Option(
        ..., "--mapping-json", exists=True, help="Path to mapping.json"
    ),
    min_pident: float = typer.Option(30.0, "--min-pident"),
    max_evalue: float = typer.Option(1e-5, "--max-evalue"),
    min_len: int = typer.Option(60, "--min-len"),
    write_files: bool = typer.Option(
        True, "--write-files/--no-write-files", help="Write score.json + score.tsv"
    ),
):

    cfg = ScoreConfig(
        min_pident=min_pident, max_evalue=max_evalue, min_align_len=min_len
    )
    res = score_mapping_json(mapping_json, cfg)

    out_dir = mapping_json.parent
    out_json = out_dir / "score.json"
    out_tsv = out_dir / "score.tsv"

    if write_files:
        res.to_json(out_json)
        res.to_tsv(out_tsv)

    console.print(
        Panel.fit(
            "[bold green]Mapping summary[/bold green]\n\n"
            f"[bold]Reference:[/bold] {res.ref_assembly}\n"
            f"[bold]Target:[/bold] {res.target_assembly}\n"
            f"[bold]Focal gene:[/bold] {res.gene_name}\n\n"
            f"[bold]Thresholds:[/bold] min_pident={cfg.min_pident}, max_evalue={cfg.max_evalue}, min_len={cfg.min_align_len}\n"
            f"[bold]Score:[/bold] {res.score_fraction_passing:.3f}  ({res.hits_passing}/{res.total_neighbors} passing)\n"
            + (
                f"\n[bold]Wrote:[/bold] {out_json}\n[bold]Wrote:[/bold] {out_tsv}"
                if write_files
                else ""
            )
        )
    )

    table = Table(title="Neighbor hits")
    table.add_column("Neighbor", style="bold")
    table.add_column("Status")
    table.add_column("pident", justify="right")
    table.add_column("len", justify="right")
    table.add_column("evalue", justify="right")

    for r in res.rows:
        if r.pass_thresholds:
            status = "[green]PASS[/green]"
        elif r.hit_found:
            status = "[yellow]HIT[/yellow]"
        else:
            status = "[red]NOHIT[/red]"

        table.add_row(
            r.neighbor_name,
            status,
            "" if r.pident is None else f"{r.pident:.1f}",
            "" if r.length is None else str(r.length),
            "" if r.evalue is None else f"{r.evalue:.2g}",
        )

    console.print(table)


@app.command("format-blastp-table")
def format_blastp_table(
    blastp_tsv: Path = typer.Option(
        ..., "--blastp-tsv", exists=True, help="Path to blastp_ortholog.tsv"
    ),
):
    """
    Convert raw blastp_ortholog.tsv into a clean reciprocal BLASTP result table.
    Keeps raw BLAST values only. No biological interpretation.
    Writes:
      - blastp_ortholog_clean.tsv
      - blastp_ortholog_readable.txt
    """
    out_dir = blastp_tsv.parent
    out_clean_tsv = out_dir / "blastp_ortholog_clean.tsv"
    out_txt = out_dir / "blastp_ortholog_readable.txt"

    rows = []

    with blastp_tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            clean = {
                "query": row.get("query", ""),
                "reference": row.get("reference", ""),
                "target": row.get("target", ""),
                "forward_best_hit": row.get("forward_best_hit", ""),
                "forward_pident": row.get("forward_pident", ""),
                "forward_qcovs": row.get("forward_qcovs", ""),
                "forward_length": row.get("forward_length", ""),
                "forward_evalue": row.get("forward_evalue", ""),
                "forward_bitscore": row.get("forward_bitscore", ""),
                "reverse_best_hit": row.get("reverse_best_hit", ""),
                "reverse_description": row.get("reverse_description", ""),
                "reverse_pident": row.get("reverse_pident", ""),
                "reverse_qcovs": row.get("reverse_qcovs", ""),
                "reverse_evalue": row.get("reverse_evalue", ""),
                "reciprocal_pass": row.get("reciprocal_pass", ""),
                "status": row.get("status", ""),
            }
            rows.append(clean)

    headers = [
        "query",
        "reference",
        "target",
        "forward_best_hit",
        "forward_pident",
        "forward_qcovs",
        "forward_length",
        "forward_evalue",
        "forward_bitscore",
        "reverse_best_hit",
        "reverse_description",
        "reverse_pident",
        "reverse_qcovs",
        "reverse_evalue",
        "reciprocal_pass",
        "status",
    ]

    # Clean TSV for Excel / R / Python
    with out_clean_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Readable aligned TXT for quick viewing
    txt_lines = []
    txt_lines.append(
        f"{'query':<8} "
        f"{'target':<15} "
        f"{'forward_hit':<18} "
        f"{'f_pid':>7} "
        f"{'f_qcov':>7} "
        f"{'f_eval':>10} "
        f"{'reverse_hit':<18} "
        f"{'r_pid':>7} "
        f"{'r_qcov':>7} "
        f"{'r_eval':>10} "
        f"{'status'}\n"
    )
    txt_lines.append("-" * 130 + "\n")

    for r in rows:
        txt_lines.append(
            f"{r['query']:<8} "
            f"{r['target']:<15} "
            f"{r['forward_best_hit']:<18} "
            f"{r['forward_pident']:>7} "
            f"{r['forward_qcovs']:>7} "
            f"{r['forward_evalue']:>10} "
            f"{r['reverse_best_hit']:<18} "
            f"{r['reverse_pident']:>7} "
            f"{r['reverse_qcovs']:>7} "
            f"{r['reverse_evalue']:>10} "
            f"{r['status']}\n"
        )

    out_txt.write_text("".join(txt_lines), encoding="utf-8")

    console.print(
        Panel.fit(
            "[bold green]Formatted raw reciprocal BLASTP table[/bold green]\n\n"
            f"[bold]Input:[/bold] {blastp_tsv}\n"
            f"[bold]Clean TSV:[/bold] {out_clean_tsv}\n"
            f"[bold]Readable TXT:[/bold] {out_txt}"
        )
    )


@app.command("blastp-ortholog")
def blastp_ortholog(
    query: str = typer.Option(
        ..., "--query", help="Query name in data/queries, e.g. DCN, FMOD, BGN, ASPN"
    ),
    ref: str = typer.Option(
        "human", "--ref", help="Reference proteome name, default human"
    ),
    targets: str = typer.Option(
        ..., "--targets", help="Comma-separated target proteome names"
    ),
    evalue: float = typer.Option(1e-5, "--evalue"),
    root: Path = typer.Option(default_project_root(), "--root"),
):
    """
    BLASTP-only ortholog screen.

    Human query protein -> target proteome
    best target protein -> reference proteome
    reciprocal pass if reverse best hit returns the original query protein.
    """
    paths = ProjectPaths(root=root)
    paths.ensure()

    query_name = query.upper()
    query_faa = paths.data_dir / "queries" / f"{query_name}.faa"

    if not query_faa.exists():
        console.print(f"[bold red]Query FASTA not found:[/bold red] {query_faa}")
        raise typer.Exit(code=1)

    run_dir = paths.new_run_dir()
    out_dir = run_dir / "blastp_ortholog" / query_name
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_proteome = paths.data_dir / "ncbi" / "proteomes" / ref / "protein.faa"
    if not ref_proteome.exists():
        console.print(
            f"[bold red]Reference proteome not found:[/bold red] {ref_proteome}"
        )
        raise typer.Exit(code=1)

    ref_db_prefix = (
        paths.data_dir / "ncbi" / "proteomes" / ref / "blastdb" / "proteome_prot"
    )
    if not any(ref_db_prefix.parent.glob(ref_db_prefix.name + ".*")):
        make_prot_db(ref_proteome, ref_db_prefix)

    target_list = [t.strip() for t in targets.split(",") if t.strip()]

    rows = []
    tsv_lines = [
        "query\treference\ttarget\tforward_best_hit\tforward_pident\tforward_qcovs\tforward_length\tforward_evalue\tforward_bitscore\t"
        "reverse_best_hit\treverse_description\treverse_pident\treverse_qcovs\treverse_evalue\treciprocal_pass\tstatus\n"
    ]

    # Reference accession used to validate the reverse hit.
    query_header, _ = fasta_read_one(query_faa)
    query_id = normalize_acc(query_header)

    accepted_ref_ids = {
        normalize_acc(query_id),
        normalize_acc(query_id).split(".")[0],
    }

    for target in target_list:
        try:
            target_proteome = (
                paths.data_dir / "ncbi" / "proteomes" / target / "protein.faa"
            )

            if not target_proteome.exists():
                status = "TARGET_PROTEOME_MISSING"
                rows.append({"target": target, "status": status})
                tsv_lines.append(
                    f"{query_name}\t{ref}\t{target}\t\t\t\t\t\t\t\t\t\t\t0\t{status}\n"
                )
                continue

            target_db_prefix = (
                paths.data_dir
                / "ncbi"
                / "proteomes"
                / target
                / "blastdb"
                / "proteome_prot"
            )
            if not any(target_db_prefix.parent.glob(target_db_prefix.name + ".*")):
                make_prot_db(target_proteome, target_db_prefix)

            forward = run_blastp_table_top1(query_faa, target_db_prefix, evalue=evalue)

            if forward is None:
                status = "NO_FORWARD_HIT"
                rows.append({"target": target, "status": status})
                tsv_lines.append(
                    f"{query_name}\t{ref}\t{target}\t\t\t\t\t\t\t\t\t\t\t0\t{status}\n"
                )
                continue

            target_hit_id = forward["sseqid"]
            target_record = fasta_get_record_by_id(target_proteome, target_hit_id)

            if target_record is None:
                status = "TARGET_HIT_SEQUENCE_NOT_FOUND"
                rows.append(
                    {
                        "target": target,
                        "forward_best_hit": target_hit_id,
                        "status": status,
                    }
                )
                tsv_lines.append(
                    f"{query_name}\t{ref}\t{target}\t{target_hit_id}\t{forward['pident']:.2f}\t{forward['qcovs']:.1f}\t"
                    f"{forward['length']}\t{forward['evalue']:.2g}\t{forward['bitscore']:.1f}\t\t\t\t\t0\t{status}\n"
                )
                continue

            hit_header, hit_seq = target_record
            back_query = (
                out_dir
                / "target_best_hits"
                / f"{target}_{normalize_acc(target_hit_id)}.faa"
            )
            fasta_write_one(back_query, hit_header, hit_seq)

            reverse = run_blastp_table_top1(back_query, ref_db_prefix, evalue=evalue)

            if reverse is None:
                reverse_best = ""
                reverse_description = ""
                reciprocal_pass = False
                status = "NO_REVERSE_HIT"
                reverse_pident = ""
                reverse_qcovs = ""
                reverse_evalue = ""
            else:
                reverse_best = reverse["sseqid"]
                reverse_best_norm = normalize_acc(reverse_best)

                reverse_description = fasta_get_full_header_by_id(
                    ref_proteome, reverse_best
                )

                reciprocal_pass = (
                    reverse_best_norm in accepted_ref_ids
                    or reverse_best_norm.split(".")[0] in accepted_ref_ids
                )

                if reciprocal_pass:
                    status = "PASS"
                else:
                    status = f"REVERSE_BEST_IS_{reverse_best_norm}"

                reverse_pident = f"{reverse['pident']:.2f}"
                reverse_qcovs = f"{reverse['qcovs']:.1f}"
                reverse_evalue = f"{reverse['evalue']:.2g}"

            row = {
                "query": query_name,
                "reference": ref,
                "target": target,
                "forward_best_hit": target_hit_id,
                "forward_pident": forward["pident"],
                "forward_qcovs": forward["qcovs"],
                "forward_length": forward["length"],
                "forward_evalue": forward["evalue"],
                "forward_bitscore": forward["bitscore"],
                "reverse_best_hit": reverse_best,
                "reverse_description": reverse_description,
                "reverse_pident": None if reverse is None else reverse["pident"],
                "reverse_qcovs": None if reverse is None else reverse["qcovs"],
                "reverse_evalue": None if reverse is None else reverse["evalue"],
                "reciprocal_pass": reciprocal_pass,
                "status": status,
            }
            rows.append(row)

            tsv_lines.append(
                f"{query_name}\t{ref}\t{target}\t"
                f"{target_hit_id}\t{forward['pident']:.2f}\t{forward['qcovs']:.1f}\t"
                f"{forward['length']}\t{forward['evalue']:.2g}\t{forward['bitscore']:.1f}\t"
                f"{reverse_best}\t{reverse_description}\t{reverse_pident}\t{reverse_qcovs}\t{reverse_evalue}\t"
                f"{1 if reciprocal_pass else 0}\t{status}\n"
            )

        except Exception as e:
            status = f"FAILED: {e}"
            rows.append({"target": target, "status": status})
            tsv_lines.append(
                f"{query_name}\t{ref}\t{target}\t\t\t\t\t\t\t\t\t\t\t0\t{status}\n"
            )

    out_json = out_dir / "blastp_ortholog.json"
    out_tsv = out_dir / "blastp_ortholog.tsv"

    out_json.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    out_tsv.write_text("".join(tsv_lines), encoding="utf-8")

    console.print(
        Panel.fit(
            "[bold green]BLASTP ortholog screen finished[/bold green]\n\n"
            f"[bold]Query:[/bold] {query_name}\n"
            f"[bold]Reference:[/bold] {ref}\n"
            f"[bold]Targets:[/bold] {len(target_list)}\n"
            f"[bold]TSV:[/bold] {out_tsv}\n"
            f"[bold]JSON:[/bold] {out_json}"
        )
    )


def _safe_cache_name(x: str) -> str:
    return "".join(c if c.isalnum() or c in {"_", "-", "."} else "_" for c in str(x))


def _query_label_for_cache(source_label: str) -> str:
    return _safe_cache_name(source_label)


def get_canonical_focal_cache_dir(
    paths: ProjectPaths,
    target_assembly: str,
    gene_name: str,
    query_label: str,
    flank: int,
) -> Path:
    return (
        paths.genomes_dir
        / target_assembly
        / "canonical_focal_cache"
        / _safe_cache_name(gene_name)
        / query_label
        / f"flank_{flank}"
    )


def load_cached_focal_hit(cache_dir: Path) -> dict | None:
    meta = cache_dir / "focal_hit.json"
    if not meta.exists():
        return None
    try:
        return json.loads(meta.read_text(encoding="utf-8"))
    except Exception:
        return None


def save_cached_focal_hit(cache_dir: Path, payload: dict):
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = cache_dir / "focal_hit.json"
    meta.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def ensure_region_db(region_fna: Path, db_prefix: Path):
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not any(db_prefix.parent.glob(db_prefix.name + ".*")):
        make_nucl_db(region_fna, db_prefix)


def resolve_focal_query_candidates(
    *,
    paths: ProjectPaths,
    ref_files,
    gene_name: str,
    focal_accession: str | None,
    focal_faa: Path | None,
    focal_fallback_accessions: list[str],
) -> list[tuple[str, Path]]:

    # Returns a list of (source_label, fasta_path) candidates in priority order.

    candidates: list[tuple[str, Path]] = []
    seen_labels: set[str] = set()

    def _add_candidate(label: str, fa: Path):
        if label in seen_labels:
            return
        seen_labels.add(label)
        candidates.append((label, fa))

    if focal_faa is not None:
        _add_candidate(f"{gene_name}__faa__{focal_faa.stem}", focal_faa)

    elif focal_accession is not None:
        fa = fetch_protein_fasta(
            accession=focal_accession,
            out_fasta=paths.data_dir / "queries" / f"{focal_accession}.faa",
        ).fasta_path
        _add_candidate(f"{gene_name}__acc__{focal_accession}", fa)

    else:
        focal_prot = protein_id_for_gene_symbol_first(ref_files.gff3, gene_name)
        if focal_prot:
            fa = fetch_protein_fasta(
                accession=focal_prot,
                out_fasta=paths.data_dir / "queries" / f"{focal_prot}.faa",
            ).fasta_path
            _add_candidate(f"{gene_name}__acc__{focal_prot}", fa)

    for acc in focal_fallback_accessions:
        fa = fetch_protein_fasta(
            accession=acc,
            out_fasta=paths.data_dir / "queries" / f"{acc}.faa",
        ).fasta_path
        _add_candidate(f"{gene_name}__acc__{acc}", fa)

    if not candidates:
        raise RuntimeError(
            f"Could not resolve any focal protein query for {gene_name}. "
            f"Provide --focal-accession or --focal-faa."
        )

    return candidates


def resolve_reference_block(
    *,
    paths: ProjectPaths,
    ref_assembly: str,
    ref_files,
    gene_name: str,
    flank: int,
    n: int,
    focal_query_candidates: list[tuple[str, Path]],
    use_canonical_focal_cache: bool,
):
    """Resolve the reference focal block, using annotation first and sequence fallback if needed."""

    found = find_gene_by_name_fast(ref_files.gff3, gene_name=gene_name)

    if found:
        region_start = max(1, found.start_1based - flank)
        region_end = found.end_1based + flank

        genes = genes_in_region_fast(
            ref_files.gff3,
            found.contig,
            region_start,
            region_end,
        )

        focal = None
        for gene in genes:
            if found.gene_id and gene.gene_id == found.gene_id:
                focal = gene
                break

        if focal is None:
            for gene in genes:
                if gene.gene_name and gene.gene_name.lower() == gene_name.lower():
                    focal = gene
                    break

        if focal is None:
            raise RuntimeError(
                "Reference focal gene was found but not recovered inside its own region."
            )

        block = neighbors_around_focal(genes, focal, n_each_side=n)
        return block, focal, region_start, region_end, "name"

    if not ref_files.genome_fasta:
        raise RuntimeError(
            f"Reference {ref_assembly} has no genome FASTA, cannot resolve focal locus by sequence."
        )

    hit, region_start, region_end, _region_fa, _region_db, used_label, _used_cache = (
        resolve_species_focal_region(
            paths=paths,
            assembly=ref_assembly,
            species_files=ref_files,
            gene_name=gene_name,
            flank=flank,
            focal_query_candidates=focal_query_candidates,
            use_canonical_focal_cache=use_canonical_focal_cache,
        )
    )

    genes = genes_in_region_fast(
        ref_files.gff3,
        hit.sseqid,
        region_start,
        region_end,
    )

    if not genes:
        raise RuntimeError(
            f"No annotated genes found in fallback region for {ref_assembly}."
        )

    hit_mid = (hit.start_1based + hit.end_1based) // 2
    focal = min(
        genes,
        key=lambda gene: abs(((gene.start_1based + gene.end_1based) // 2) - hit_mid),
    )

    block = neighbors_around_focal(genes, focal, n_each_side=n)
    return block, focal, region_start, region_end, f"sequence:{used_label}"


@app.command()
def run():

    console.print("[yellow]run() not implemented yet.[/yellow] Use: synteny init")


if __name__ == "__main__":
    app()
