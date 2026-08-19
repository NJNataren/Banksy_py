#!/bin/bash

#SBATCH --job-name=00_xenium_clustering
#SBATCH --array=0-2%2 # runs at most 2 tasks at a time out of 3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # increase from 8 - parallelisation of Leiden clustering will benefit
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=logs/xenium_clustering_%x_%A_%a.out
#SBATCH --error=logs/xenium_clustering_%x_%A_%a.err
#SBATCH --partition=sacgf
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=nathalie.nataren@adelaide.edu.au

#SBATCH --export=NONE

set -euo pipefail

REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
mkdir -p logs
## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

#CONFIG=/hpcfs/users/a1210419/Banksy_py/config/vbct/CK_skin_res.json
CONFIG_DIR="config/00_clustering/vbct/large"
CONFIGS=($CONFIG_DIR/*.json)
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

python 00_xenium_clustering_clean_adata.py --config $CONFIG

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"
