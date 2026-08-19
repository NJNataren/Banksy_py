#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

PROJECT="ptmt_pc55"
CONFIG_DIR="config/00_clustering/${PROJECT}"
CLUSTER_PREFIX="labels_scaled_gaussian_pc55_nc0.20_r"

shopt -s nullglob
CONFIGS=("${CONFIG_DIR}"/*.json)

if (( ${#CONFIGS[@]} == 0 )); then
  echo "No JSON configs found in ${CONFIG_DIR}" >&2
  exit 1
fi

for CONFIG in "${CONFIGS[@]}"; do
  SAMPLE="$(basename "${CONFIG}" .json)"
  SAMPLE_DIR="data/xenium/processed/${PROJECT}/${SAMPLE}"
  OUTPUT_DIR="data/xenium/output/${PROJECT}/${SAMPLE}/clustree_qc"

  mapfile -t CLUSTER_CSVS < <(compgen -G "${SAMPLE_DIR}/${SAMPLE}_cell_cluster_id_across_clustering_res_*.csv" || true)
  if (( ${#CLUSTER_CSVS[@]} == 0 )); then
    echo "Skipping ${SAMPLE}: no cluster CSV found in ${SAMPLE_DIR}" >&2
    continue
  fi
  if (( ${#CLUSTER_CSVS[@]} > 1 )); then
    echo "Skipping ${SAMPLE}: multiple cluster CSVs found; choose one explicitly" >&2
    printf '  %s
' "${CLUSTER_CSVS[@]}" >&2
    continue
  fi

  CLUSTER_CSV="${CLUSTER_CSVS[0]}"
  echo "Running clustree for ${SAMPLE} using ${CLUSTER_CSV}"
  Rscript 01a_clustree_cluster_resolution_qc.R     --cluster_csv "${CLUSTER_CSV}"     --dataset_name "${SAMPLE}"     --cluster_prefix "${CLUSTER_PREFIX}"     --output_dir "${OUTPUT_DIR}"
done
