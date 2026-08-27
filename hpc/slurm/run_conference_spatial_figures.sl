#!/bin/bash

#SBATCH --job-name=conference_spatial_figures
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/conference_spatial_figures_%x_%j.out
#SBATCH --error=logs/conference_spatial_figures_%x_%j.err
#SBATCH --partition=sacgf

#SBATCH --export=None

set -euo pipefail

REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"
CONFIG="${CONFIG:-config/conference_figures/ck_skin_res_unfiltered_r1p10.json}"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
echo "Config: ${CONFIG}"

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
    --config "${CONFIG}"

echo "Conference figure job finished: $(date '+%Y-%m-%d %H:%M:%S')"
