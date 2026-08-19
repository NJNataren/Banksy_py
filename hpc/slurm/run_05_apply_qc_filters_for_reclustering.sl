#!/bin/bash

#SBATCH --job-name=05_xenium_qc_filter
#SBATCH --array=0-7%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf

set -euo pipefail

REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
mkdir -p logs

###############################
#   Conda environment setup   #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIG_DIR="${CONFIG_DIR:-config/05_apply_qc_filters/vbct}"
# Example override:
# sbatch --array=0-7%2 --export=ALL,CONFIG_DIR=config/05_apply_qc_filters/vbct hpc/slurm/run_05_apply_qc_filters_for_reclustering.sl
CONFIGS=("${CONFIG_DIR}"/*.json)

if (( ${#CONFIGS[@]} == 0 )); then
    echo "No JSON configs found in ${CONFIG_DIR}."
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID >= ${#CONFIGS[@]} )); then
    echo "Array task ID ${SLURM_ARRAY_TASK_ID} is outside the config count (${#CONFIGS[@]}) for ${CONFIG_DIR}."
    exit 0
fi

CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

echo "Task ${SLURM_ARRAY_TASK_ID} using config: ${CONFIG}"
echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Config directory: ${CONFIG_DIR}"
echo "Config file: ${CONFIG}"
echo "---------------------------------"
echo "Config contents:"
cat "${CONFIG}"
echo "================================="

python 05_apply_qc_filters_for_reclustering.py --config "${CONFIG}"

echo "QC filtering finished: $(date '+%Y-%m-%d %H:%M:%S')"
