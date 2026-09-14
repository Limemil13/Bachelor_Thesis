#!/usr/bin/env python3

import argparse
import re
import statistics
from pathlib import Path


def read_fasta(path):
    records = []
    header = None
    seq = []

    with Path(path).open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

    if header is not None:
        records.append((header, "".join(seq)))

    return records


def write_fasta(records, path):
    with Path(path).open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            out.writelines(seq[i : i + 80] + "\n" for i in range(0, len(seq), 80))


def sanitize(text):
    text = text.strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def parse_ncbi_fields(header):
    """Return organism, GeneID and isoform from keyed NCBI header brackets."""
    fields = {}
    bare = []
    for content in re.findall(r"\[([^\]]+)\]", header):
        if "=" in content:
            key, value = content.split("=", 1)
            fields[key.strip().lower()] = value.strip()
        else:
            bare.append(content.strip())
    organism = fields.get("organism")
    if not organism and bare:
        organism = bare[0]
    return (
        organism or "unknown_species",
        fields.get("geneid", ""),
        fields.get("isoform", ""),
    )


def parse_accession(header):
    return header.split()[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--min-abs", type=int, default=180)
    parser.add_argument("--max-abs", type=int, default=600)
    parser.add_argument("--lower", type=float, default=0.70)
    parser.add_argument("--upper", type=float, default=1.30)
    parser.add_argument("--one-per-species", action="store_true")
    parser.add_argument("--one-per-locus", action="store_true")
    args = parser.parse_args()

    records = read_fasta(args.input)

    parsed = []
    for header, seq in records:
        acc = parse_accession(header)
        species, gene_id, isoform = parse_ncbi_fields(header)
        clean_species = sanitize(species)
        length = len(seq)
        clean_header = f"{clean_species}|{args.gene}|{acc}"
        parsed.append(
            {
                "old_header": header,
                "new_header": clean_header,
                "seq": seq,
                "acc": acc,
                "species": species,
                "clean_species": clean_species,
                "gene_id": gene_id,
                "isoform": isoform,
                "length": length,
                "status": "",
                "reason": "",
            }
        )

    abs_ok = [r for r in parsed if args.min_abs <= r["length"] <= args.max_abs]

    if not abs_ok:
        raise SystemExit("No sequences passed absolute length filter.")

    median_len = statistics.median([r["length"] for r in abs_ok])
    min_rel = args.lower * median_len
    max_rel = args.upper * median_len

    rel_ok = []
    for r in parsed:
        if not (args.min_abs <= r["length"] <= args.max_abs):
            r["status"] = "excluded"
            r["reason"] = "absolute_length"
        elif not (min_rel <= r["length"] <= max_rel):
            r["status"] = "excluded"
            r["reason"] = "relative_length"
        else:
            rel_ok.append(r)

    selected = rel_ok

    if args.one_per_species and args.one_per_locus:
        raise SystemExit("Choose only one of --one-per-species or --one-per-locus.")

    if args.one_per_species or args.one_per_locus:
        grouped = {}
        for r in rel_ok:
            if args.one_per_locus:
                # GeneID is the stable NCBI locus key.  Falling back to the
                # accession avoids collapsing unrelated records when GeneID is absent.
                group_key = r["gene_id"] or r["acc"]
            else:
                group_key = r["clean_species"]
            grouped.setdefault(group_key, []).append(r)

        selected_ids = set()
        for _, group in grouped.items():
            best = sorted(
                group, key=lambda r: (abs(r["length"] - median_len), r["length"])
            )[0]
            selected_ids.add(best["acc"])

        for r in rel_ok:
            if r["acc"] in selected_ids:
                r["status"] = "kept"
                r["reason"] = "representative"
            else:
                r["status"] = "excluded"
                r["reason"] = (
                    "extra_isoform_same_locus"
                    if args.one_per_locus
                    else "extra_candidate_same_species"
                )

        selected = [r for r in rel_ok if r["acc"] in selected_ids]
    else:
        for r in rel_ok:
            r["status"] = "kept"
            r["reason"] = "length_ok"

    write_fasta([(r["new_header"], r["seq"]) for r in selected], args.output)

    with Path(args.report).open("w") as out:
        out.write(
            "gene\taccession\tspecies\tgene_id\tisoform\tlength\tstatus\treason\told_header\tnew_header\n"
        )
        for r in parsed:
            out.write(
                f"{args.gene}\t{r['acc']}\t{r['species']}\t{r['gene_id']}\t{r['isoform']}\t{r['length']}\t"
                f"{r['status']}\t{r['reason']}\t{r['old_header']}\t{r['new_header']}\n"
            )

    print(f"Input sequences: {len(parsed)}")
    print(f"Median length after absolute filter: {median_len}")
    print(f"Kept sequences: {len(selected)}")
    print(f"Wrote: {args.output}")
    print(f"Wrote: {args.report}")


if __name__ == "__main__":
    main()
