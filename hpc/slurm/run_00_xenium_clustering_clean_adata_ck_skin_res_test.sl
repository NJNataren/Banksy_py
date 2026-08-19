#!/bin/bash

#SBATCH --job-name=00_ck_skin_res_clean_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.bootstrap.out
#SBATCH --error=%x_%j.bootstrap.err
#SBATCH --partition=sacgf
#SBATCH --export=NONE

set -euo pipefail

REPO_DIR="${REPO_DIR:-/scratchdata1/users/a1210419/Banksy_py}"

echo "Submitting directory: ${SLURM_SUBMIT_DIR:-unknown}"
echo "Changing to repository directory: ${REPO_DIR}"
cd "${REPO_DIR}"
echo "Running from: $(pwd)"
mkdir -p logs

CONFIG="config/00_clustering/vbct/small/CK_skin_res.json"
PROJECT="vbct"
DATASET_NAME="CK_skin_res"
RES_STR="0.70_0.80_0.90_1.00"
OUTPUT_DIR="data/xenium/output/${PROJECT}/${DATASET_NAME}"
RUN_TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUTPUT_DIR}/slurm_test_runs/${RUN_TIMESTAMP}_${SLURM_JOB_ID}"

mkdir -p "${RUN_DIR}"

exec > >(tee "${RUN_DIR}/00_xenium_clustering_clean_adata.out") \
     2> >(tee "${RUN_DIR}/00_xenium_clustering_clean_adata.err" >&2)

echo "================================="
echo "Job started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Config file: ${CONFIG}"
echo "Run directory: ${RUN_DIR}"
echo "---------------------------------"
echo "Config contents:"
cat "${CONFIG}"
echo "================================="

source /hpcfs/users/a1210419/miniforge3/etc/profile.d/conda.sh
conda activate banksy

python 00_xenium_clustering_clean_adata.py --config "${CONFIG}"

echo "================================="
echo "Inspecting clean AnnData metadata transfer"
python - <<PY
import anndata as ad

path = "data/xenium/processed/${PROJECT}/${DATASET_NAME}/adata_expression_clean_${DATASET_NAME}_with_banksy_clusters_${RES_STR}.h5ad"
adata = ad.read_h5ad(path)
label_cols = [
    col for col in adata.obs.columns
    if col.startswith("labels_scaled_gaussian_pc30_nc0.20")
]

print(f"Read: {path}")
print(adata)
print("BANKSY label columns:", label_cols)
print("obsm keys:", list(adata.obsm.keys()))

required_obsm = [
    "xy",
    "spatial",
    "X_umap_scaled_gaussian_pc30_nc0.20",
]
missing_obsm = [key for key in required_obsm if key not in adata.obsm]
if missing_obsm:
    raise SystemExit(f"Missing expected obsm keys: {missing_obsm}")
if len(label_cols) != 4:
    raise SystemExit(f"Expected 4 BANKSY label columns, found {len(label_cols)}")
PY

echo "Job finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Timestamped test output written to: ${RUN_DIR}"
