#!/bin/bash
#SBATCH --job-name=10_cell_type_spatial
#SBATCH --array=0-7%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=sacgf
set -euo pipefail
REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"
cd "$REPO_DIR"
source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy
CONFIG=config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json
# Derive task order from the annotation config.
SAMPLE=$(python - "$CONFIG" "${SLURM_ARRAY_TASK_ID:?Run through sbatch}" <<'PYCODE'
import json
import sys
from pathlib import Path
cfg = json.loads(Path(sys.argv[1]).read_text())
annotation = json.loads(Path(cfg['annotation_config']).read_text())
print(annotation['samples'][int(sys.argv[2])]['sample'])
PYCODE
)
echo "Sample: $SAMPLE"
python -u 10_squidpy_recluster_spatial_analysis.py --config "$CONFIG" --sample "$SAMPLE" --validate-only
python -u 10_squidpy_recluster_spatial_analysis.py --config "$CONFIG" --sample "$SAMPLE"
