#!/usr/bin/env python3
"""Cluster TBLASTN HSPs into candidate genomic windows."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Cluster:
    # Group nearby HSPs on one contig and strand.
    seqid: str
    strand: str
    start: int
    end: int
    best_by_query: dict[str, float] = field(default_factory=dict)
    hsp_count: int = 0

    def add(self, query: str, start: int, end: int, bitscore: float) -> None:
        # Expand the span and retain the best score from each query.
        self.start = min(self.start, start)
        self.end = max(self.end, end)
        self.best_by_query[query] = max(bitscore, self.best_by_query.get(query, 0.0))
        self.hsp_count += 1

    @property
    def score(self) -> float:
        # Count only the best HSP per query when ranking a window.
        return sum(self.best_by_query.values())


def parse_args() -> argparse.Namespace:
    # Keep clustering and padding settings configurable.
    parser = argparse.ArgumentParser()
    parser.add_argument("hits", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--padding", type=int, default=100_000)
    parser.add_argument("--cluster-gap", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hits: dict[tuple[str, str], list[tuple[int, int, str, float]]] = defaultdict(list)

    with args.hits.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            # Normalize minus-strand subject coordinates before clustering.
            if len(row) < 12:
                continue
            query, seqid = row[0], row[1]
            sstart, send = int(row[8]), int(row[9])
            strand = "+" if sstart <= send else "-"
            start, end = sorted((sstart, send))
            hits[(seqid, strand)].append((start, end, query, float(row[11])))

    clusters: list[Cluster] = []
    for (seqid, strand), records in hits.items():
        current: Cluster | None = None
        for start, end, query, bitscore in sorted(records):
            # Start a new cluster after the configured genomic gap.
            if current is None or start > current.end + args.cluster_gap:
                current = Cluster(seqid, strand, start, end)
                clusters.append(current)
            current.add(query, start, end, bitscore)

    # Rank clusters and send only the top windows to Miniprot.
    ranked = sorted(
        clusters,
        key=lambda item: (item.score, len(item.best_by_query), item.hsp_count),
        reverse=True,
    )[: args.top]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["rank", "seqid", "start", "end", "strand", "score", "query_count", "hsp_count"]
        )
        for rank, cluster in enumerate(ranked, start=1):
            # Add flanking sequence for spliced alignment and clip at base 1.
            writer.writerow(
                [
                    rank,
                    cluster.seqid,
                    max(1, cluster.start - args.padding),
                    cluster.end + args.padding,
                    cluster.strand,
                    f"{cluster.score:.2f}",
                    len(cluster.best_by_query),
                    cluster.hsp_count,
                ]
            )


if __name__ == "__main__":
    main()
