#!/bin/bash

#SBATCH --job-name=01_ck_skin_res_qc_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=sacgf

set -euo pipefail

REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"
CONFIG="config/01_QC/vbct/small/CK_skin_res.json"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
mkdir -p logs

echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Config file: ${CONFIG}"
echo "---------------------------------"
echo "Config contents:"
cat "${CONFIG}"
echo "================================="

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

python 01_QC_xenium_spatial_clean_clustered.py --config "${CONFIG}"

echo "================================="
echo "CK_skin_res QC test finished: $(date '+%Y-%m-%d %H:%M:%S')"
