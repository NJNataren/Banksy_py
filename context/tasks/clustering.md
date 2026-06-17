# Task: Xenium BANKSY Clustering

## Goal

Maintain and extend the config-driven Xenium clustering workflow for BANKSY-based spatial clustering of metastatic melanoma Xenium samples, while preserving clean gene-expression values for downstream marker analysis.

## Current Status

- Active clustering script: `01_xenium_clustering_clean_adata_test.py`.
- The older `01_xenium_clustering.py` and `01_xenium_clustering.ipynb` are no longer the active implementation in this working tree.
- The active script reads one JSON config via `--config`.
- Clustering configs live under `config/01_clustering/`, split by study/group and sample.
- PTMT configs currently live under `config/01_clustering/ptmt/` and include `raw_subdir` so samples can read raw files from run-specific folders such as `ptmt/Run_1`, `ptmt/Run_2`, and `ptmt/Run_3`.
- Main PTMT Slurm entrypoint: `run_01_xenium_clustering_test.sl`.
- That Slurm script currently uses `CONFIG_DIR="config/01_clustering/ptmt"` and runs `python 01_xenium_clustering_clean_adata_test.py --config $CONFIG`.
- If running all current PTMT configs, update/check the Slurm array range against the number of JSON files in `config/01_clustering/ptmt/`.

## Key Files

- `01_xenium_clustering_clean_adata_test.py`
- `config/01_clustering/`
- `config/01_clustering/ptmt/`
- `run_01_xenium_clustering_test.sl`
- `banksy/`
- `banksy_utils/`

## Config Contract

Each clustering config should define:

- `dataset_name`: sample name used for input and output filenames.
- `pc_label`: number of principal components as a filename/config label, for example `"35"`.
- `lambda_label`: BANKSY lambda label, for example `"0.20"`.
- `res_label`: list of Leiden/BANKSY resolutions as strings.
- `nbr_weight_decay`: neighbour weighting scheme, commonly `"scaled_gaussian"`.
- `coord_keys`: coordinate columns/obsm key, usually `["x", "y", "xy"]`.
- `raw_subdir`: optional subdirectory under `data/xenium/raw_data/`, for example `"ptmt/Run_2"`.

The script reads raw input from:

```text
data/xenium/raw_data/<raw_subdir>/<dataset_name>_raw.h5ad
```

If `raw_subdir` is omitted, the script falls back to:

```text
data/xenium/raw_data/<dataset_name>_raw.h5ad
```

## Workflow Outputs

For each sample, the active script:

- filters zero-count cells using `nCount_Xenium > 0`;
- downcasts the AnnData object to reduce memory load;
- saves raw counts in `adata.layers["counts"]`;
- normalizes and log-transforms expression;
- saves a clean expression AnnData before BANKSY feature expansion;
- runs BANKSY graph construction, PCA/UMAP, and Leiden clustering across configured resolutions;
- saves `adata_spatial_<dataset_name>_<resolution>.h5ad` objects;
- saves the BANKSY dictionary and `results_df`;
- exports a cell-cluster identity table across resolutions;
- copies BANKSY cluster labels back onto the clean expression AnnData by cell ID;
- runs `sc.tl.rank_genes_groups()` on clean expression values, not on the BANKSY-expanded spatial object;
- saves marker CSVs, marker plots, total-count violin plots, and the clean AnnData with BANKSY labels/marker rankings.

## Constraints

- Keep clustering config-driven via `--config path/to/sample.json`.
- Use `raw_subdir` or another config field for run-specific raw data locations; do not hard-code run folders in the script.
- Avoid hard-coding sample names or paths when configs already provide them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as machine-specific or potentially large.
- Do not run full clustering jobs locally unless explicitly requested.
- Prefer targeted syntax checks or argument/config smoke tests when feasible.
- Preserve clean-expression marker analysis: biological marker ranking/plots should use `clean_adata`, while BANKSY spatial objects are for clustering structure.

## Next Steps

- Rename or promote `01_xenium_clustering_clean_adata_test.py` once the clean-AnnData workflow is considered production-ready.
- Before running a Slurm array, check that `#SBATCH --array` matches the number of intended configs.
- Consider reducing leftover notebook-export inspection cells/comments in the active script when stabilizing it.
