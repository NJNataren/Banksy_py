#!/bin/bash

#SBATCH --job-name=00_ck_dirs_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/xenium_clustering_%x_%j.out
#SBATCH --error=logs/xenium_clustering_%x_%j.err
#SBATCH --partition=sacgf
#SBATCH --export=NONE

set -euo pipefail

REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"
CONFIG="config/00_clustering/vbct/small/CK_skin_res.json"
PROJECT="vbct"
DATASET_NAME="CK_skin_res"

PROCESSED_DIR="data/xenium/processed/${PROJECT}/${DATASET_NAME}"
OUTPUT_DIR="data/xenium/output/${PROJECT}/${DATASET_NAME}"
QC_DIR="data/xenium/output/${PROJECT}/QC_testing/${DATASET_NAME}"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"

echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Config file: ${CONFIG}"
echo "---------------------------------"
echo "Directory state before running script 00:"
for path in "${PROCESSED_DIR}" "${OUTPUT_DIR}" "${QC_DIR}"; do
    if [[ -d "${path}" ]]; then
        echo "EXISTS: ${path}"
    else
        echo "MISSING: ${path}"
    fi
done
echo "---------------------------------"
echo "Config contents:"
cat "${CONFIG}"
echo "================================="

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

python 00_xenium_clustering_clean_adata.py --config "${CONFIG}"

echo "================================="
echo "Directory state after running script 00:"
for path in "${PROCESSED_DIR}" "${OUTPUT_DIR}" "${QC_DIR}"; do
    if [[ -d "${path}" ]]; then
        echo "EXISTS: ${path}"
    else
        echo "MISSING: ${path}"
    fi
done
echo "Job finished: $(date '+%Y-%m-%d %H:%M:%S')"
