#!/bin/bash
#SBATCH --job-name=10_cell_type_spatial
#SBATCH --array=0-7%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf
set -euo pipefail
REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"
cd "$REPO_DIR"
source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy
CONFIG=config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json
# Derive task order from the samples this config uses (cohort order by default).
SAMPLE=$(python - "$CONFIG" "${SLURM_ARRAY_TASK_ID:?Run through sbatch}" <<'PYCODE'
import sys
from workflow_config import load_cohort, load_json, select_samples
cfg = load_json(sys.argv[1])
samples = select_samples(load_cohort(cfg['cohort']), cfg.get('samples'))
print(samples[int(sys.argv[2])]['sample'])
PYCODE
)
echo "Sample: $SAMPLE"
python -u 10_cell_type_spatial_analysis.py --config "$CONFIG" --sample "$SAMPLE" --validate-only
python -u 10_cell_type_spatial_analysis.py --config "$CONFIG" --sample "$SAMPLE"
