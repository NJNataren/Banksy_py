#!/bin/bash

#SBATCH --job-name=03_xenium_dotplot_export
#SBATCH --array=0-15%1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf

set -euo pipefail

###############################
#   Conda environment setup   #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

# Override at submit time when needed, for example:
# sbatch --export=CONFIG_DIR=config/dotplot/export_summary/vbct_small run_03_xenium_dotplot_export_from_config.sl
CONFIG_DIR="${CONFIG_DIR:-config/dotplot/export_summary/vbct_small}"
shopt -s nullglob
CONFIGS=("$CONFIG_DIR"/*.json)

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "No JSON configs found in $CONFIG_DIR" >&2
    exit 1
fi

if [[ $SLURM_ARRAY_TASK_ID -ge ${#CONFIGS[@]} ]]; then
    echo "Array task $SLURM_ARRAY_TASK_ID is out of range for ${#CONFIGS[@]} configs in $CONFIG_DIR" >&2
    exit 1
fi

CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Config directory: $CONFIG_DIR"
echo "Config file: $CONFIG"
echo "---------------------------------"
echo "Config contents:"
cat "$CONFIG"
echo "================================="

python 03_export_dotplot_data_from_config.py --config "$CONFIG"

echo "Dotplot export finished: $(date '+%Y-%m-%d %H:%M:%S')"
