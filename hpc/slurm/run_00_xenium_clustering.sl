#!/bin/bash

#SBATCH --job-name=00_xenium_s_clustering
#SBATCH --array=0-4 #limit to 2 so you're not using all the nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/xenium_clustering_%x_%A_%a.out
#SBATCH --error=logs/xenium_clustering_%x_%A_%a.err
#SBATCH --partition=sacgf


#SBATCH --export=None

set -euo pipefail

REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
echo "Checking required project directories:"
ls -ld data data/xenium data/xenium/processed

## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

#CONFIG=/hpcfs/users/a1210419/Banksy_py/config/vbct/CK_skin_res.json #local testing
CONFIG_DIR="${CONFIG_DIR:-config/00_clustering/vbct/small}"
# Example override:
# sbatch --export=ALL,CONFIG_DIR=config/00_clustering/ptmt_pc55 hpc/slurm/run_00_xenium_clustering.sl
CONFIGS=("${CONFIG_DIR}"/*.json)

if (( SLURM_ARRAY_TASK_ID >= ${#CONFIGS[@]} )); then
    echo "Array task ID ${SLURM_ARRAY_TASK_ID} is outside the config count (${#CONFIGS[@]}) for ${CONFIG_DIR}."
    exit 0
fi

CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

## Print a timestamp for each job and the config contents
echo "Task $SLURM_ARRAY_TASK_ID using config: $CONFIG"
echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Config file: $CONFIG"
echo "---------------------------------"
echo " Config contents:"
cat $CONFIG
echo "================================="

#########################
#	Run clustering	#
#########################

#python 00_xenium_clustering_clean_adata.py --config $CONFIG
python 00_xenium_clustering_clean_adata.py --config $CONFIG

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"
