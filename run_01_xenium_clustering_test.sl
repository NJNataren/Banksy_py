#!/bin/bash

#SBATCH --job-name=01_xenium_clustering
#SBATCH --array=0-4%2 #limit to 2 so you're not using all the nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/xenium_clustering_%x_%A_%a.out
#SBATCH --error=logs/xenium_clustering_%x_%A_%a.err
#SBATCH --partition=sacgf


#SBATCH --export=None

set -euo pipefail
## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

#CONFIG=/hpcfs/users/a1210419/Banksy_py/config/vbct/CK_skin_res.json #local testing
#CONFIG_DIR="config/clustering/vbct/small"
CONFIG_DIR="config/clustering/ptmt"
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
#	Run the QC script	#
#########################

#python 01_xenium_clustering.py --config $CONFIG
python 01_xenium_clustering_clean_adata_test.py --config $CONFIG

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"
