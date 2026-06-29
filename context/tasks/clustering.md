# Task: Xenium BANKSY Clustering

## Goal

Maintain and extend the config-driven Xenium clustering workflow for BANKSY-based spatial clustering of metastatic melanoma Xenium samples, while preserving clean gene-expression values for downstream marker analysis.

## Current Status

- Active clustering script: `01_xenium_clustering_clean_adata_test.py`.
- The older `01_xenium_clustering.py` and `01_xenium_clustering.ipynb` are no longer the active implementation in this working tree.
- The active script reads one JSON config via `--config`.
- Clustering configs live under `config/01_clustering/`, split by study/group and sample.
- Processed and output paths are now explicitly project-scoped in configs and scripts, for example `data/xenium/processed/vbct/<dataset_name>/` and `data/xenium/output/ptmt/<dataset_name>/`.
- The previous helper-function approach for inferring/scoping project paths was removed; configs should state the project/path intent directly.
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
- `project`: project/study folder used for processed and output paths, usually `"vbct"` or `"ptmt"`.
- `raw_subdir`: optional subdirectory under `data/xenium/raw_data/`, for example `"ptmt/Run_2"`.

The script reads raw input from:

```text
data/xenium/raw_data/<raw_subdir>/<dataset_name>_raw.h5ad
```

If `raw_subdir` is omitted, the script falls back to:

```text
data/xenium/raw_data/<dataset_name>_raw.h5ad
```

The script writes project-scoped outputs to:

```text
data/xenium/processed/<project>/<dataset_name>/
data/xenium/output/<project>/<dataset_name>/
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

Important saved clean AnnData outputs:

- clean pre-BANKSY expression object: `data/xenium/processed/<project>/<dataset_name>/adata_clean_<dataset_name>.h5ad`;
- clean expression object with all transferred BANKSY labels and marker-ranking results: `data/xenium/processed/<project>/<dataset_name>/adata_clean_with_banksy_labels_<dataset_name>.h5ad`.

Plotting details fixed recently:

- BANKSY labels are treated as ordered categorical labels so cluster ordering is numeric rather than lexicographic, preventing `0, 1, 10, ...` ordering in count and marker plots.
- Full-figure spatial, UMAP, PCA, and connection plots use a shared categorical colour mapping so cluster colours stay consistent across panels.
- The colour mapping supports more than 20 clusters, tested with `10850_run_3_1818_AMACR_neg` at resolution `1.10`.

## Clean Object Transfer Recommendations

When updating `01_xenium_clustering_clean_adata_test.py`, use it as the place where BANKSY-derived metadata is copied onto the clean expression AnnData, because this script has access to the clean object, `banksy_dict`, `results_df`, and per-resolution spatial objects at the same time.

Recommended metadata to copy into the clean object:

- BANKSY cluster labels for each configured resolution, already stored as columns such as `labels_scaled_gaussian_pc30_nc0.20_r1.00`.
- BANKSY UMAP embeddings, stored in `.obsm` with explicit names such as `X_umap_scaled_gaussian_pc30_nc0.20`. Resolution usually should not be part of the UMAP key unless the embedding itself is resolution-specific.
- BANKSY PCA embeddings, if available and useful for debugging/reproducibility, with explicit `.obsm` keys such as `X_pca_scaled_gaussian_pc30_nc0.20`.
- Spatial coordinates standardized as `clean_adata.obsm["spatial"] = clean_adata.obs[["x", "y"]].to_numpy()` while preserving the existing `xy` convention.
- Human-readable cluster annotation columns when `new_labels` is available, for example `labels_scaled_gaussian_pc30_nc0.20_r1.00_ann`.
- QC-relevant `.obs` columns such as `nCount_Xenium`, `nFeature_Xenium`, `transcript_counts`, `cell_area`, `nucleus_area`, and aggregate control/codeword count columns.

Do not copy BANKSY-expanded expression values into the clean expression matrix. In particular, do not replace `clean_adata.X` with `adata_spatial.X`, and do not add neighbour-derived `_nbr_0`/`_nbr_1` features as biological expression values. Avoid storing large BANKSY graphs or full `banksy_dict` contents in the clean object unless a downstream workflow explicitly needs them.

## Constraints

- Keep clustering config-driven via `--config path/to/sample.json`.
- Use `raw_subdir` for run-specific raw data locations; do not hard-code run folders in the script.
- Keep processed/output project routing explicit in JSON configs or simple script path construction using `project`; do not reintroduce hidden helper functions for path rewriting.
- Avoid hard-coding sample names or paths when configs already provide them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as machine-specific or potentially large.
- Do not run full clustering jobs locally unless explicitly requested.
- Prefer targeted syntax checks or argument/config smoke tests when feasible.
- Preserve clean-expression marker analysis: biological marker ranking/plots should use `clean_adata`, while BANKSY spatial objects are for clustering structure.

## Next Steps

- Rename or promote `01_xenium_clustering_clean_adata_test.py` once the clean-AnnData workflow is considered production-ready.
- Before running a Slurm array, check that `#SBATCH --array` matches the number of intended configs.
- Consider reducing leftover notebook-export inspection cells/comments in the active script when stabilizing it.
- Re-run clustering for any samples whose downstream workflows need the new project-scoped processed paths or the latest colour/order fixes.
