from pathlib import Path

GENES = ["BGN", "FMOD", "OGN", "EPYC", "PRELP"]


def create_phylogeny_run(
    root: Path,
    run_name: str = "slrp_phylogeny_2026_06",
) -> dict[str, Path]:
    """Create the directory structure for a phylogeny/domain run."""
    run_dir = root / "analyses" / run_name

    folders = {
        "run": run_dir,
        "blastp": run_dir / "01_blastp",
        "candidates": run_dir / "02_candidates",
        "domains": run_dir / "03_domains",
        "alignments": run_dir / "04_alignments",
        "trees": run_dir / "05_trees",
        "itol": run_dir / "06_itol",
        "logs": run_dir / "logs",
    }

    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)

    return folders
