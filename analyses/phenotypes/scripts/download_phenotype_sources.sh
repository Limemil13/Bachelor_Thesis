#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
RAW="$REPO/analyses/phenotypes/raw"
mkdir -p "$RAW"

curl -L --fail --retry 4 \
  https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt \
  -o "$RAW/MGI_GenePheno.rpt"
curl -L --fail --retry 4 \
  https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt \
  -o "$RAW/MGI_PhenoGenoMP.rpt"
curl -L --fail --retry 4 \
  https://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology \
  -o "$RAW/MPheno_OBO.ontology"

HPO_RELEASE=v2026-06-23
curl -L --fail --retry 4 \
  "https://github.com/obophenotype/human-phenotype-ontology/releases/download/$HPO_RELEASE/genes_to_phenotype.txt" \
  -o "$RAW/HPO_${HPO_RELEASE}_genes_to_phenotype.txt"
curl -L --fail --retry 4 \
  "https://github.com/obophenotype/human-phenotype-ontology/releases/download/$HPO_RELEASE/genes_to_disease.txt" \
  -o "$RAW/HPO_${HPO_RELEASE}_genes_to_disease.txt"

echo "Phenotype source downloads complete: $RAW"
