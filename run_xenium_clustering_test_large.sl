#!/bin/bash

#SBATCH --job-name=01_xenium_clustering
#SBATCH --array=0-2 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/xenium_clustering_%x_%A_%a.out
#SBATCH --error=logs/xenium_clustering_%x_%A_%a.err
#SBATCH --partition=sacgf

#SBATCH --export=None


## Load conda environment

###############################
#	Conda environment setup	  #
###############################

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

#CONFIG=/hpcfs/users/a1210419/Banksy_py/config/vbct/CK_skin_res.json
CONFIG_DIR="config/vbct/large"
CONFIGS=($CONFIG_DIR/*.json)
CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

## Print a timestamp for each job and the config contents
echo "Task $SLURM_ARRAY_TASK_ID using config: $CONFIG"
echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
#echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Config file: $CONFIG"
echo "---------------------------------"
echo " Config contents:"
cat $CONFIG
echo "================================="

#########################
#	Run the QC script	#
#########################

python 01_xenium_clustering.py --config $CONFIG

echo "Sample finished: $(date '+%Y-%m-%d %H:%M:%S')
