"""Parse one explicit SignalP 6 prediction table against the canonical ID map."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
SIGNALP_DIR = BASE / "signal_peptides"
DEFAULT_MAPPING = SIGNALP_DIR / "input" / "signalp_canonical_id_mapping.tsv"
OUT_DIR = SIGNALP_DIR / "tables"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_cs(cs_text: str) -> tuple[str, str]:
    # Extract the optional cleavage position and probability from SignalP text.
    pos_match = re.search(r"CS pos:\s*([0-9]+-[0-9]+)", cs_text)
    prob_match = re.search(r"Pr:\s*([0-9.]+)", cs_text)
    return (
        pos_match.group(1) if pos_match else "",
        prob_match.group(1) if prob_match else "",
    )


def quality(
    prediction: str, sp_score: float, cleavage_probability: str
) -> tuple[str, str, str, str]:
    # Add check/watch labels without changing SignalP's SP/no-SP prediction.
    notes: list[str] = []
    strict_check = "no"
    watch = "no"

    if prediction != "SP":
        return (
            "no_signal_peptide",
            "yes",
            "yes",
            "No SignalP Sec/SPI signal peptide predicted",
        )

    label = "strong_signal_peptide"
    if sp_score < 0.90:
        label = "weak_signal_peptide"
        strict_check = "yes"
        watch = "yes"
        notes.append("Low SignalP SP probability")

    if cleavage_probability:
        cs_probability = float(cleavage_probability)
        if cs_probability < 0.80:
            if label == "strong_signal_peptide":
                label = "low_cleavage_probability"
            strict_check = "yes"
            watch = "yes"
            notes.append("Low cleavage-site probability")
        elif cs_probability < 0.85:
            if label == "strong_signal_peptide":
                label = "watch_cleavage_probability"
            watch = "yes"
            notes.append("Slightly lower cleavage-site probability")

    return label, strict_check, watch, "; ".join(notes)


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    # Require exact agreement between returned IDs and the submitted ID map.
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--job-id", required=True)
    parser.add_argument(
        "--model-mode", required=True, choices=("fast", "slow-sequential")
    )
    parser.add_argument("--run-date", required=True)
    parser.add_argument(
        "--output-prefix",
        default="signalp",
        help="Output filename prefix in the tables directory (default: signalp)",
    )
    args = parser.parse_args()

    mapping_rows = read_tsv(args.mapping)
    mapping = {row["signalp_id"]: row for row in mapping_rows}
    if not mapping_rows or len(mapping) != len(mapping_rows):
        raise ValueError("SignalP mapping must contain at least one row and unique IDs")

    # Ignore format comments but reject malformed or duplicate result rows.
    predictions: dict[str, list[str]] = {}
    with args.input.open(encoding="utf-8", errors="strict") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(f"Malformed SignalP line: {line}")
            if parts[0] in predictions:
                raise ValueError(f"Duplicate SignalP ID: {parts[0]}")
            predictions[parts[0]] = parts

    if set(predictions) != set(mapping):
        raise ValueError(
            "SignalP IDs differ from canonical mapping; "
            f"missing={sorted(set(mapping) - set(predictions))}; "
            f"extra={sorted(set(predictions) - set(mapping))}"
        )

    source = (
        f"SignalP-6.0 canonical {args.model_mode} ({args.run_date}; job {args.job_id})"
    )
    # Preserve raw scores beside the derived check/watch labels.
    rows: list[dict[str, str]] = []
    for mapped in mapping_rows:
        seq_id = mapped["signalp_id"]
        parts = predictions[seq_id]
        prediction = parts[1]
        other_score = float(parts[2])
        sp_score = float(parts[3])
        cs_text = parts[4] if len(parts) > 4 else ""
        cleavage_site, cleavage_probability = parse_cs(cs_text)
        quality_label, strict_check, watch, notes = quality(
            prediction, sp_score, cleavage_probability
        )
        rows.append(
            {
                "gene": mapped["gene"],
                "species": mapped["species"],
                "accession": mapped["accession"],
                "signalp_id": seq_id,
                "signal_peptide": "yes" if prediction == "SP" else "no",
                "signalp_prediction": prediction,
                "other_score": f"{other_score:.6f}",
                "sp_score": f"{sp_score:.6f}",
                "cleavage_site": cleavage_site,
                "cleavage_probability": cleavage_probability,
                "signalp_quality": quality_label,
                "strict_check": strict_check,
                "watch": watch,
                "notes": notes,
                "signalp_source": source,
                "signalp_job_id": args.job_id,
                "signalp_model_mode": args.model_mode,
                "signalp_rerun_required": "no",
            }
        )

    fields = list(rows[0])
    write_tsv(OUT_DIR / f"{args.output_prefix}_summary_clean.tsv", rows, fields)
    write_tsv(
        OUT_DIR / f"{args.output_prefix}_candidates_to_check.tsv",
        [row for row in rows if row["strict_check"] == "yes"],
        fields,
    )
    write_tsv(
        OUT_DIR / f"{args.output_prefix}_candidates_watch.tsv",
        [row for row in rows if row["watch"] == "yes"],
        fields,
    )

    print(f"Read: {args.input}")
    print(f"Validated predictions: {len(rows)}")
    print(f"Signal peptide yes: {sum(row['signal_peptide'] == 'yes' for row in rows)}")
    print(f"Signal peptide no: {sum(row['signal_peptide'] == 'no' for row in rows)}")
    print(f"Strict check: {sum(row['strict_check'] == 'yes' for row in rows)}")
    print(f"Watch: {sum(row['watch'] == 'yes' for row in rows)}")
    print(f"Source: {source}")


if __name__ == "__main__":
    main()
