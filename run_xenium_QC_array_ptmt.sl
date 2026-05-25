#!/bin/bash


#SBATCH --job-name=00_QC_testing_xenium_spatial
#SBATCH --array=0-7%2
#number of samples to process
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf

set -euo pipefail

## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIG_DIR="config/QC/ptmt"
CONFIGS=($CONFIG_DIR/*.json)
CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

echo "Task $SLURM_ARRAY_TASK_ID using config: $CONFIG"


#########################
#	Run the QC script	#
#########################

python 00_QC_testing_xenium_spatial_PTMT_v2.py --config $CONFIG
