#!/bin/bash

#SBATCH --job-name=01_xenium_color_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/xenium_clustering_color_test_%x_%j.out
#SBATCH --error=logs/xenium_clustering_color_test_%x_%j.err
#SBATCH --partition=sacgf

#SBATCH --export=None

set -euo pipefail

###############################
#   Conda environment setup   #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIG="config/01_clustering/testing/10850_run_3_1818_AMACR_neg_r1p10_color_test.json"

mkdir -p logs

## Print a timestamp for the job and the config contents
echo "Using config: $CONFIG"
echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Config file: $CONFIG"
echo "---------------------------------"
echo " Config contents:"
cat "$CONFIG"
echo "================================="

#############################
#   Run clustering script    #
#############################

python 01_xenium_clustering_clean_adata_test.py --config "$CONFIG"

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"
