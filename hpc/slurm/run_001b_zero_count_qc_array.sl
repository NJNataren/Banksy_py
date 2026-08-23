#!/bin/bash

#SBATCH --job-name=001b_zero_count_qc
#SBATCH --array=0-23%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf

set -euo pipefail

REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"
CONFIG_DIR="${CONFIG_DIR:-config/00_clustering/ptmt}"

# Example submissions:
# sbatch --array=0-23%4 --export=ALL,CONFIG_DIR=config/00_clustering/ptmt hpc/slurm/run_001b_zero_count_qc_array.sl
# sbatch --array=0-23%4 --export=ALL,CONFIG_DIR=config/00_clustering/ptmt_pc55 hpc/slurm/run_001b_zero_count_qc_array.sl


echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
mkdir -p logs

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIGS=("${CONFIG_DIR}"/*.json)
N_CONFIGS="${#CONFIGS[@]}"

if [[ "${N_CONFIGS}" -eq 0 || ! -e "${CONFIGS[0]}" ]]; then
    echo "No JSON configs found in ${CONFIG_DIR}"
    exit 1
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${N_CONFIGS}" ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} is outside config count ${N_CONFIGS}; exiting cleanly."
    exit 0
fi

CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Config directory: ${CONFIG_DIR}"
echo "Config count: ${N_CONFIGS}"
echo "Config file: ${CONFIG}"
echo "---------------------------------"
echo "Config contents:"
cat "${CONFIG}"
echo "================================="

python 001b_plot_zero_count_from_raw_adata.py --config "${CONFIG}"

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"
