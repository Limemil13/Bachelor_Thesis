#!/usr/bin/env python3
"""Curate the compact 16-species LUM panel from reciprocal BLASTP evidence."""

from __future__ import annotations

import csv
import textwrap
from pathlib import Path

from build_reciprocal_blastp_panel import read_fasta

BASE = Path(__file__).resolve().parents[1]
PROJECT_ROOT = BASE.parents[1]
DIAGNOSTIC_DIR = BASE / "diagnostics/LUM" / "reciprocal_blastp"
EVIDENCE = DIAGNOSTIC_DIR / "reciprocal_blastp.tsv"
ALL_CANDIDATES = DIAGNOSTIC_DIR / "LUM_candidates_all.faa"
OUT_FASTA = BASE / "candidates" / "LUM_domain_input.faa"
OUT_SELECTION = BASE / "candidates" / "LUM_candidate_selection.tsv"

SPECIES_ORDER = (
    "human",
    "mouse",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "bos_taurus",
    "monodelphis_domestica",
    "chicken",
    "anolis_carolinensis",
    "frog",
    "zebrafish",
    "lepisosteus_oculatus",
    "latimeria_chalumnae",
    "callorhinchus_milii",
    "scyliorhinus_canicula",
    "lamprey",
    "amphioxus",
    "whale_shark",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sequence_by_species(path: Path) -> dict[str, tuple[str, str]]:
    result: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        species = header.split()[0].split("|")[0]
        if species in result:
            raise ValueError(f"Duplicate candidate species in {path}: {species}")
        result[species] = (header, sequence)
    return result


def proteome_record(species: str, accession: str) -> tuple[str, str]:
    proteome = PROJECT_ROOT / "data" / "ncbi" / "proteomes" / species / "protein.faa"
    for header, sequence in read_fasta(proteome):
        if header.split()[0] == accession:
            return header, sequence
    raise KeyError(f"{accession} not found in {proteome}")


def main() -> None:
    evidence = {row["species"]: row for row in read_tsv(EVIDENCE)}
    candidates = sequence_by_species(ALL_CANDIDATES)
    if set(evidence) != set(SPECIES_ORDER) or set(candidates) != set(SPECIES_ORDER):
        raise ValueError(
            "LUM reciprocal evidence must contain the established 17-species screen"
        )

    selected: list[tuple[str, str]] = []
    selection_rows: list[dict[str, str]] = []
    for species in SPECIES_ORDER:
        row = evidence[species]
        raw_accession = row["forward_best_hit"]
        chosen_accession = raw_accession
        include = row["reciprocal_pass"] == "yes"
        rationale = "reciprocal_best_hit_returns_human_LUM"
        review_flag = ""

        if species == "bos_taurus":
            # The top BLAST hit is an X1 model with a 43-aa noncanonical
            # N-terminal extension.  The curated RefSeq precursor is identical
            # across the complete query-aligned region and is the cleaner
            # representative protein.
            chosen_accession = "NP_776359.1"
            rationale = "curated_RefSeq_precursor_preferred_over_43aa_extended_X1_model"
        elif species == "lamprey":
            review_flag = (
                "shared_with_previous_FMOD_panel;deep_lineage_class_II_ambiguity"
            )
        elif species == "amphioxus":
            include = False
            rationale = (
                "top_hit_is_654aa_non_LUM_and_reverse_best_hit_is_XP_005270514.1"
            )
            review_flag = "no_defensible_full_length_reciprocal_LUM_candidate"

        if include:
            if chosen_accession == raw_accession:
                _, sequence = candidates[species]
            else:
                _, sequence = proteome_record(species, chosen_accession)
            header = (
                f"{species}|LUM|{chosen_accession} "
                f"forward_pident={float(row['forward_pident']):.2f} "
                f"forward_qcovs={float(row['forward_qcovs']):.1f} "
                f"reverse_hit={row['reverse_best_hit']}"
            )
            selected.append((header, sequence))
            selected_length = str(len(sequence))
        else:
            selected_length = ""

        selection_rows.append(
            {
                "gene": "LUM",
                "species": species,
                "raw_best_hit": raw_accession,
                "raw_best_hit_length": row["target_protein_length"],
                "reverse_best_hit": row["reverse_best_hit"],
                "reciprocal_pass": row["reciprocal_pass"],
                "included_in_canonical_panel": "yes" if include else "no",
                "selected_accession": chosen_accession if include else "",
                "selected_length": selected_length,
                "selection_rationale": rationale,
                "manual_review_flag": review_flag,
            }
        )

    if len(selected) != 16:
        raise ValueError(f"Expected 16 selected LUM proteins, found {len(selected)}")
    accessions = [header.split()[0].split("|")[2] for header, _ in selected]
    if len(accessions) != len(set(accessions)):
        raise ValueError("Duplicate LUM accessions in canonical panel")

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with OUT_FASTA.open("w", encoding="utf-8") as handle:
        for header, sequence in selected:
            handle.write(f">{header}\n")
            handle.write("\n".join(textwrap.wrap(sequence, 80)) + "\n")
    with OUT_SELECTION.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(selection_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(selection_rows)

    print(f"Wrote {OUT_FASTA} ({len(selected)} proteins)")
    print(f"Wrote {OUT_SELECTION}")


if __name__ == "__main__":
    main()
