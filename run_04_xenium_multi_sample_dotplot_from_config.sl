#!/bin/bash

#SBATCH --job-name=04_xenium_multi_dotplot
#SBATCH --array=0-0%1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
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
# sbatch --export=CONFIG_DIR=config/04_plot_dotplot/archive/vbct_hpc_or_old_layout/vbct_small run_xenium_multi_sample_dotplot_from_config.sl
#CONFIG_DIR="${CONFIG_DIR:-config/04_plot_dotplot/archive/vbct_hpc_or_old_layout/vbct_small}"
CONFIG_DIR="${CONFIG_DIR:-config/04_plot_dotplot/archive/vbct_hpc_or_old_layout/vbct_large}"
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

python 04_plot_multi_sample_dotplot_from_config.py --config "$CONFIG"

echo "Multi-sample dotplot finished: $(date '+%Y-%m-%d %H:%M:%S')"
