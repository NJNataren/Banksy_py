#!/bin/bash
#SBATCH --job-name=08_cell_type_annotation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=sacgf

set -euo pipefail

REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"
CONFIG="${CONFIG:-config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions.json}"

cd "$REPO_DIR"
mkdir -p logs

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

EXTRA_ARGS=()
if [[ "${OVERWRITE:-0}" == "1" ]]; then
    EXTRA_ARGS+=(--overwrite)
fi

echo "Repository: $REPO_DIR"
echo "Config: $CONFIG"
echo "Overwrite existing outputs: ${OVERWRITE:-0}"
echo "Started: $(date)"

# Script 08 writes shared count tables and a combined eight-sample UMAP, so all
# configured samples run in one job rather than competing array tasks.
python -u 08_annotate_cell_types_and_plot_umaps.py \
    --config "$CONFIG" \
    --input-root "$REPO_DIR" \
    --output-root "$REPO_DIR" \
    "${EXTRA_ARGS[@]}"

echo "Finished: $(date)"
