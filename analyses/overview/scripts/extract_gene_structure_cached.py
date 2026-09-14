#!/usr/bin/env python3
"""Run SynVoy's gene-structure extractor while parsing each GFF only once.

The upstream extractor loops over GFFs and then genes but reparses the complete
GFF inside every gene call.  A one-entry cache preserves the implementation and
results while avoiding seven identical parses per reference annotation.
"""

from __future__ import annotations

import argparse
import functools
import importlib.util
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--upstream-script", required=True, type=Path)
    known, remaining = parser.parse_known_args()

    spec = importlib.util.spec_from_file_location(
        "synvoy_extract_gene_structure", known.upstream_script
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {known.upstream_script}")
    workflow = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(workflow)

    workflow.parse_gff = functools.lru_cache(maxsize=1)(workflow.parse_gff)
    sys.argv = [str(known.upstream_script), *remaining]
    workflow.main()


if __name__ == "__main__":
    main()
