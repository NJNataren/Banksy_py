#!/bin/bash

#SBATCH --job-name=01_xenium_clustering
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/xenium_clustering_%x_%j.out
#SBATCH --error=logs/xenium_clustering_%x_%j.err
#SBATCH --partition=sacgf


#SBATCH --export=None


## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

CONFIG=/hpcfs/users/a1210419/Banksy_py/config/vbct/clustering/small/BE_brain_non_res.json #local testing
## CONFIG_DIR="config/clustering/vbct/small"
## CONFIGS=($CONFIG_DIR/*.json)
## CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

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

python 01_xenium_clustering.py --config $CONFIG

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')"