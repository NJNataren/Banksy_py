#!/bin/bash

#SBATCH --job-name=00_QC_testing_xenium_spatial
#SBATCH --array=0-7%2
#number of samples to process
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf

#SBATCH --export=None


## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIG_DIR="config/vbct"
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

python 00_QC_testing_xenium_spatial.py --config $CONFIG
