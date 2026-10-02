#!/usr/bin/env python3
"""Compare gene neighborhoods around named loci in GFF3 annotations.

The comparison uses exact, case-insensitive gene-symbol matches. This is a
conservative first-pass synteny test: shared symbols support a corresponding
locus, while missing symbols can also reflect annotation gaps or renamed genes.

Example
-------
python compare_gff_neighborhoods.py \
  --locus chicken_BGN chicken.gff BGN \
  --locus human_BGN human.gff BGN \
  --locus human_ASPN human.gff ASPN@NC_000009.12 \
  --flank-count 10 \
  --output-dir neighborhood_results
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from urllib.parse import unquote


@dataclass(frozen=True)
class Gene:
    seqid: str
    start: int
    end: int
    strand: str
    symbol: str
    gene_id: str
    description: str
    biotype: str


@dataclass(frozen=True)
class Locus:
    label: str
    gff_path: Path
    selector: str
    assembly: str
    focal: Gene
    upstream: tuple[Gene, ...]
    downstream: tuple[Gene, ...]

    @property
    def genes_in_genomic_order(self) -> tuple[Gene, ...]:
        return self.upstream + (self.focal,) + self.downstream

    @property
    def flanking_symbols(self) -> tuple[str, ...]:
        genes = self.upstream + self.downstream
        return tuple(gene.symbol.upper() for gene in genes)

    @property
    def oriented_flanking_symbols(self) -> tuple[str, ...]:
        symbols = tuple(
            gene.symbol.upper()
            for gene in self.genes_in_genomic_order
            if gene != self.focal
        )
        return tuple(reversed(symbols)) if self.focal.strand == "-" else symbols


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def read_assembly_name(path: Path) -> str:
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("#!genome-build-accession"):
                return raw_line.split(maxsplit=1)[1].strip()
            if not raw_line.startswith("#"):
                break
    return "not reported"


def gene_symbol(attributes: dict[str, str]) -> str:
    symbol = (
        attributes.get("gene")
        or attributes.get("Name")
        or attributes.get("gene_name")
        or attributes.get("locus_tag")
        or attributes.get("ID", "unnamed")
    )
    return symbol.removeprefix("gene-")


def read_genes(path: Path, include_noncoding: bool) -> list[Gene]:
    genes: list[Gene] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            biotype = attributes.get("gene_biotype", "not reported")
            if not include_noncoding and biotype not in {
                "not reported",
                "protein_coding",
            }:
                continue
            genes.append(
                Gene(
                    seqid=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    symbol=gene_symbol(attributes),
                    gene_id=attributes.get("ID", ""),
                    description=attributes.get("description", ""),
                    biotype=biotype,
                )
            )
    return genes


def find_focal_gene(genes: list[Gene], selector: str) -> Gene:
    symbol_selector, separator, seqid_selector = selector.rpartition("@")
    if not separator:
        symbol_selector, seqid_selector = selector, ""
    wanted = symbol_selector.upper()
    matches = [
        gene
        for gene in genes
        if wanted in {gene.symbol.upper(), gene.gene_id.removeprefix("gene-").upper()}
        and (not seqid_selector or gene.seqid == seqid_selector)
    ]
    if not matches:
        raise ValueError(f"No exact gene match for {selector!r}")
    if len(matches) > 1:
        locations = ", ".join(
            f"{gene.seqid}:{gene.start}-{gene.end}" for gene in matches[:10]
        )
        raise ValueError(f"Multiple exact matches for {selector!r}: {locations}")
    return matches[0]


def gene_distance(focal: Gene, neighbor: Gene) -> int:
    if neighbor.end < focal.start:
        return focal.start - neighbor.end
    if neighbor.start > focal.end:
        return neighbor.start - focal.end
    return 0


def make_locus(
    label: str,
    path: Path,
    selector: str,
    genes: list[Gene],
    flank_count: int,
    max_distance_bp: int | None,
) -> Locus:
    focal = find_focal_gene(genes, selector)
    chromosome_genes = sorted(
        (gene for gene in genes if gene.seqid == focal.seqid),
        key=lambda gene: (gene.start, gene.end, gene.symbol),
    )
    focal_index = chromosome_genes.index(focal)
    upstream = chromosome_genes[:focal_index]
    downstream = chromosome_genes[focal_index + 1 :]
    if max_distance_bp is not None:
        upstream = [
            gene for gene in upstream if gene_distance(focal, gene) <= max_distance_bp
        ]
        downstream = [
            gene for gene in downstream if gene_distance(focal, gene) <= max_distance_bp
        ]
    return Locus(
        label=label,
        gff_path=path,
        selector=selector,
        assembly=read_assembly_name(path),
        focal=focal,
        upstream=tuple(upstream[-flank_count:]),
        downstream=tuple(downstream[:flank_count]),
    )


def longest_common_subsequence(left: tuple[str, ...], right: tuple[str, ...]) -> int:
    previous = [0] * (len(right) + 1)
    for left_symbol in left:
        current = [0]
        for index, right_symbol in enumerate(right, start=1):
            if left_symbol == right_symbol:
                current.append(previous[index - 1] + 1)
            else:
                current.append(max(current[-1], previous[index]))
        previous = current
    return previous[-1]


def compare_loci(left: Locus, right: Locus) -> dict[str, str | int]:
    left_set = set(left.flanking_symbols)
    right_set = set(right.flanking_symbols)
    shared = sorted(left_set & right_set)
    union = left_set | right_set
    left_order = tuple(
        symbol for symbol in left.oriented_flanking_symbols if symbol in shared
    )
    right_order = tuple(
        symbol for symbol in right.oriented_flanking_symbols if symbol in shared
    )
    forward = longest_common_subsequence(left_order, right_order)
    reversed_order = longest_common_subsequence(
        left_order, tuple(reversed(right_order))
    )
    if reversed_order > forward:
        order = "reversed"
        ordered_count = reversed_order
    else:
        order = "same"
        ordered_count = forward
    return {
        "locus_a": left.label,
        "focal_a": left.focal.symbol,
        "locus_b": right.label,
        "focal_b": right.focal.symbol,
        "flanking_genes_a": len(left.flanking_symbols),
        "flanking_genes_b": len(right.flanking_symbols),
        "shared_flanking_genes": len(shared),
        "shared_symbols": ",".join(shared),
        "jaccard": f"{len(shared) / len(union):.3f}" if union else "0.000",
        "best_order_orientation": order if shared else "not assessable",
        "shared_genes_in_order": ordered_count,
        "order_fraction": f"{ordered_count / len(shared):.3f}" if shared else "0.000",
    }


def neighborhood_rows(locus: Locus) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    genes = locus.genes_in_genomic_order
    focal_index = genes.index(locus.focal)
    for index, gene in enumerate(genes):
        relative_index = index - focal_index
        rows.append(
            {
                "locus": locus.label,
                "assembly": locus.assembly,
                "gff": locus.gff_path.name,
                "focal_selector": locus.selector,
                "seqid": gene.seqid,
                "relative_index": relative_index,
                "role": "focal" if relative_index == 0 else "flanking",
                "symbol": gene.symbol,
                "gene_id": gene.gene_id,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "distance_from_focal_bp": gene_distance(locus.focal, gene),
                "biotype": gene.biotype,
                "description": gene.description,
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str | int]]) -> None:
    if not rows:
        raise ValueError(f"No rows available for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--locus",
        action="append",
        nargs=3,
        required=True,
        metavar=("LABEL", "GFF3", "GENE"),
        help=(
            "Dataset label, GFF3 path, and exact focal gene symbol; repeat as "
            "needed. Use GENE@SEQID when alternate loci share a symbol."
        ),
    )
    parser.add_argument(
        "--flank-count",
        type=int,
        default=10,
        help="Maximum number of protein-coding genes to retain on each side (default: 10).",
    )
    parser.add_argument(
        "--max-distance-bp",
        type=int,
        help="Optional maximum distance from the focal gene boundary.",
    )
    parser.add_argument(
        "--include-noncoding",
        action="store_true",
        help="Include genes explicitly annotated with non-protein-coding biotypes.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for neighborhoods.tsv and pairwise_comparison.tsv.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.flank_count < 1:
        raise ValueError("--flank-count must be at least 1")
    if args.max_distance_bp is not None and args.max_distance_bp < 0:
        raise ValueError("--max-distance-bp cannot be negative")

    gene_cache: dict[Path, list[Gene]] = {}
    loci: list[Locus] = []
    for label, raw_path, selector in args.locus:
        path = Path(raw_path).resolve()
        if path not in gene_cache:
            gene_cache[path] = read_genes(
                path, include_noncoding=args.include_noncoding
            )
        genes = gene_cache[path]
        loci.append(
            make_locus(
                label,
                path,
                selector,
                genes,
                flank_count=args.flank_count,
                max_distance_bp=args.max_distance_bp,
            )
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    neighborhood_output = args.output_dir / "neighborhoods.tsv"
    comparison_output = args.output_dir / "pairwise_comparison.tsv"
    neighborhood_data = [row for locus in loci for row in neighborhood_rows(locus)]
    comparison_data = [
        compare_loci(left, right) for left, right in combinations(loci, 2)
    ]
    write_tsv(neighborhood_output, neighborhood_data)
    write_tsv(comparison_output, comparison_data)

    for locus in loci:
        symbols = " - ".join(gene.symbol for gene in locus.genes_in_genomic_order)
        print(
            f"{locus.label}: {locus.focal.seqid}:{locus.focal.start}-{locus.focal.end} "
            f"({locus.focal.strand})\n  {symbols}"
        )
    print(f"Wrote {neighborhood_output}")
    print(f"Wrote {comparison_output}")


if __name__ == "__main__":
    main()
