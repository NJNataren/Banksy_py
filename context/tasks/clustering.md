# Task: Xenium BANKSY Clustering

## Goal

Maintain and extend the config-driven Xenium clustering workflow for BANKSY-based spatial clustering of metastatic melanoma Xenium samples, while preserving clean gene-expression values for downstream marker analysis.

## Current Status

- Active clustering script: `00_xenium_clustering_clean_adata.py`.
- Active clean-clustered QC script: `01_QC_xenium_spatial_clean_clustered.py`.
- Workflow order is now:
  - `00_xenium_clustering_clean_adata.py`: clustering plus clean AnnData generation/metadata transfer.
  - `01_QC_xenium_spatial_clean_clustered.py`: QC/inspection of the clean clustered AnnData object.
- The older `01_xenium_clustering.py` and `01_xenium_clustering.ipynb` are no longer the active implementation in this working tree.
- Both active scripts can read one JSON config via `--config`; recent notebook-friendly updates also allow no-config local CK testing fallbacks.
- Clustering configs live under `config/01_clustering/`, split by study/group and sample.
- Processed and output paths are now explicitly project-scoped in configs and scripts, for example `data/xenium/processed/vbct/<dataset_name>/` and `data/xenium/output/ptmt/<dataset_name>/`.
- The previous helper-function approach for inferring/scoping project paths was removed; configs should state the project/path intent directly.
- PTMT configs currently live under `config/01_clustering/ptmt/` and include `raw_subdir` so samples can read raw files from run-specific folders such as `ptmt/Run_1`, `ptmt/Run_2`, and `ptmt/Run_3`.
- Main clustering Slurm entrypoints include `run_00_xenium_clustering.sl` and `run_00_xenium_clustering_large.sl`.
- `run_00_xenium_clustering.sl` now anchors execution to the HPC repository path before running Python, then prints working-directory and required-directory diagnostics.
- If running all current PTMT configs, update/check the Slurm array range against the number of JSON files in `config/01_clustering/ptmt/`.

## Session Handoff

Recent work focused on making the clean AnnData object the central handoff between clustering and QC, then making both scripts easier to run as notebooks for local analysis.

- Script order/name convention was changed so clustering is step 00 and clean-clustered QC is step 01.
- Slurm wrappers and README mentions were updated earlier to point at the renamed scripts; some wrapper names may still contain legacy words, but commands should point to the active scripts.
- `00_xenium_clustering_clean_adata.py` now transfers BANKSY-derived metadata back onto the clean expression AnnData without replacing clean expression values.
- `01_QC_xenium_spatial_clean_clustered.py` was locally tested on `CK_skin_res` using the clean clustered object and completed successfully.
- Notebook versions were generated for local interactive analysis: `00_xenium_clustering_clean_adata.ipynb` and `01_QC_xenium_spatial_clean_clustered.ipynb`.
- Both active scripts were syntax-checked with `python -m py_compile`.
- Current git status at the end of the session included modified `00_xenium_clustering_clean_adata.py`, modified `run_00_xenium_clustering.sl`, deleted `config/01_clustering/vbct/large/MF_skin_non_res_roi.json`, and untracked `config/01_clustering/vbct/large/MF_skin_non_res.json`.
- The `MF_skin_non_res_roi.json` to `MF_skin_non_res.json` change appears to be an intentional user rename; do not revert it.

## Key Files

- `00_xenium_clustering_clean_adata.py`
- `00_xenium_clustering_clean_adata.ipynb`
- `01_QC_xenium_spatial_clean_clustered.py`
- `01_QC_xenium_spatial_clean_clustered.ipynb`
- `config/01_clustering/`
- `config/01_clustering/ptmt/`
- `run_00_xenium_clustering.sl`
- `run_00_xenium_clustering_large.sl`
- `run_01_xenium_QC_array_ptmt.sl`
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
- `new_labels`: used by the 01 QC script for human-readable cluster annotations. Do not use `new_labels` in the 00 clustering script.

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
- clean expression object with transferred BANKSY labels/embeddings: `data/xenium/processed/<project>/<dataset_name>/adata_expression_clean_<dataset_name>_with_banksy_clusters_<resolutions>.h5ad`.

For `CK_skin_res`, the clean clustered object was confirmed locally at:

```text
data/xenium/processed/vbct/CK_skin_res/adata_expression_clean_CK_skin_res_with_banksy_clusters_0.70_0.80_0.90_1.00.h5ad
```

It contained label columns such as `labels_scaled_gaussian_pc30_nc0.20_r0.70`, `labels_scaled_gaussian_pc30_nc0.20_r0.80`, `labels_scaled_gaussian_pc30_nc0.20_r0.90`, and `labels_scaled_gaussian_pc30_nc0.20_r1.00`, plus `.obsm` keys `xy`, `spatial`, `X_umap_scaled_gaussian_pc30_nc0.20`, and `X_pca_scaled_gaussian_pc30_nc0.20`.

Plotting details fixed recently:

- BANKSY labels are treated as ordered categorical labels so cluster ordering is numeric rather than lexicographic, preventing `0, 1, 10, ...` ordering in count and marker plots.
- Full-figure spatial, UMAP, PCA, and connection plots use a shared categorical colour mapping so cluster colours stay consistent across panels.
- The colour mapping supports more than 20 clusters, tested with `10850_run_3_1818_AMACR_neg` at resolution `1.10`.

## Clean Object Transfer Recommendations

When updating `00_xenium_clustering_clean_adata.py`, use it as the place where BANKSY-derived metadata is copied onto the clean expression AnnData, because this script has access to the clean object, `banksy_dict`, `results_df`, and per-resolution spatial objects at the same time.

Current implemented transfer behavior:

- `make_banksy_label_col()` builds dynamic BANKSY label names.
- `copy_obsm_aligned()` copies embeddings by matching observation names.
- BANKSY labels are copied to clean `.obs` for each configured resolution.
- UMAP/PCA are copied from source keys such as `reduced_pc_30_umap` and `reduced_pc_30` into clean object keys such as `.obsm["X_umap_scaled_gaussian_pc30_nc0.20"]` and `.obsm["X_pca_scaled_gaussian_pc30_nc0.20"]`.
- Spatial coordinates are standardized as `.obsm["xy"]` and `.obsm["spatial"]`.
- The script prints availability of key QC/control `.obs` columns.
- Shape and `var_names` guards confirm the clean expression object is not changed during transfer.

Do not copy BANKSY-expanded expression values into the clean expression matrix. In particular, do not replace `clean_adata.X` with `adata_spatial.X`, and do not add neighbour-derived `_nbr_0`/`_nbr_1` features as biological expression values. Avoid storing large BANKSY graphs or full `banksy_dict` contents in the clean object unless a downstream workflow explicitly needs them.

Do not add `new_labels` handling to the 00 clustering script. `new_labels` belongs in the 01 QC config/script layer where human-readable labels are used for interpretation and plotting.

## Notebook-Friendly Config Handling

Both active scripts were adjusted to work in notebooks and command-line runs.

For `00_xenium_clustering_clean_adata.py`:

- `--config` is optional.
- `parse_known_args()` avoids Jupyter kernel argument failures.
- No-config fallback targets local `CK_skin_res` testing with `project="vbct"`, `raw_subdir="vbct"`, `pc_label="30"`, `lambda_label="0.20"`, resolutions `["0.70", "0.80", "0.90", "1.00"]`, `nbr_weight_decay="scaled_gaussian"`, `coord_keys=["x", "y", "xy"]`, and `max_workers=8`.
- `max_workers` is resolved with `SLURM_CPUS_PER_TASK` first, then config `max_workers`, then default `4`.

For `01_QC_xenium_spatial_clean_clustered.py`:

- `--config` is optional.
- `parse_known_args()` avoids Jupyter kernel argument failures.
- No-config fallback also targets local `CK_skin_res` testing.
- The local fallback includes `new_labels` for the CK test clusters.
- Generic Scanpy aliases are created after loading the clean object: `X_umap_scaled_gaussian_pc30_nc0.20` to `X_umap` and `X_pca_scaled_gaussian_pc30_nc0.20` to `X_pca`.

The QC script also includes fixes from the successful local CK run:

- sparse/matrix-like `.X` values are flattened safely before percentile heatmap calculations;
- marker annotation master file fallback paths include the current May 2026 Tier 1 marker file;
- stale `threshold_passed_cat` references were replaced with existing `min_trans_passed_cat`.

## Local Validation

Useful commands that passed during this session:

```bash
python -m py_compile 00_xenium_clustering_clean_adata.py
python -m py_compile 01_QC_xenium_spatial_clean_clustered.py
MPLBACKEND=Agg python 01_QC_xenium_spatial_clean_clustered.py --config config/01_clustering/vbct/small/CK_skin_res.json
bash -n run_00_xenium_clustering.sl
```

The full local 01 QC run generated outputs under:

```text
data/xenium/output/vbct/QC_testing/CK_skin_res/
data/xenium/output/vbct/CK_skin_res/
data/xenium/processed/vbct/CK_skin_res/
```

## HPC Large-Sample Write Issue

A larger sample failed during `write_h5ad()` in `00_xenium_clustering_clean_adata.py` with h5py reporting:

```text
error message = 'No such file or directory', flags = 13, o_flags = 242
```

The failure occurred while creating a spatial AnnData `.h5ad`, not while reading the config. `config/01_clustering/vbct/large/CK_bowel_res.json` appears complete and includes the expected fields: `project`, `raw_subdir`, `dataset_name`, `pc_label`, `lambda_label`, `res_label`, `nbr_weight_decay`, `coord_keys`, and `new_labels`.

Most likely explanation: the Slurm job was running from the wrong working directory or hitting an HPC filesystem/mount/path issue. The script uses relative paths under `data/xenium`, so the working directory matters.

Implemented diagnostics:

- `run_00_xenium_clustering.sl` sets `REPO_DIR="/scratchdata1/users/a1210419/Banksy_py"` and `cd`s there before running Python.
- The wrapper prints submit directory, working directory, and `ls -ld data data/xenium data/xenium/processed`.
- The 00 script prints the spatial output path and whether the parent directory exists before writing.
- The 00 script raises a clearer `FileNotFoundError` if the spatial output parent directory is missing.

Next log lines to check after rerunning on HPC:

```text
Running from: /scratchdata1/users/a1210419/Banksy_py
Writing spatial AnnData to: data/xenium/processed/vbct/CK_bowel_res/adata_spatial_CK_bowel_res_1p1.h5ad
Spatial output parent directory exists: True
```

If the parent exists and h5py still reports no such file or directory, suspect a compute-node filesystem/HDF5 issue rather than a missing config key.

## Constraints

- Keep clustering config-driven via `--config path/to/sample.json`.
- Use `raw_subdir` for run-specific raw data locations; do not hard-code run folders in the script.
- Keep processed/output project routing explicit in JSON configs or simple script path construction using `project`; do not reintroduce hidden helper functions for path rewriting.
- Avoid hard-coding sample names or paths when configs already provide them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as machine-specific or potentially large.
- Do not run full clustering jobs locally unless explicitly requested.
- Prefer targeted syntax checks or argument/config smoke tests when feasible.
- Preserve clean-expression marker analysis: biological marker ranking/plots should use `clean_adata`, while BANKSY spatial objects are for clustering structure.

## Next Session: Clustree Resolution QC

Add an R-based clustree step as `01a_clustree_cluster_resolution_qc.R` or similar. It should sit after script 00 and alongside script 01, before marker/dotplot review with scripts 03/04. The purpose is to assess cluster stability and splitting across BANKSY Leiden resolutions before choosing a resolution for annotation and downstream biological interpretation.

Minimum input should be the script 00 cell-by-resolution cluster assignment CSV:

```text
data/xenium/processed/<project>/<dataset_name>/<dataset_name>_cell_cluster_id_across_clustering_res_<resolutions>.csv
```

That table should have one row per cell and one cluster-label column per resolution, with a consistent prefix such as `labels_scaled_gaussian_pc30_nc0.20_r`. This is enough for a first `clustree::clustree()` plot. Optional later input can be a QC `.obs` export from script 01 joined by cell ID to colour or annotate clustree nodes by metrics such as median `nCount_Xenium`, median `cell_area`, percent high-area cells, or negative-control burden.

Recommended first implementation:

- Use `01a_merge_cluster_resolution_csvs.R` when large-sample cluster assignment CSVs are split by resolution; it joins one-resolution CSVs by `index` cell ID and writes a merged clustree-ready table.
- Create a small config/CLI-driven R script that accepts `--cluster_csv`, `--dataset_name`, `--cluster_prefix`, and `--output_dir`.
- Save PNG and PDF outputs under the sample output/QC area, for example `data/xenium/output/<project>/<dataset_name>/clustree_qc/`.
- Set up local clustree dependencies with `setup_clustree_r_env.R`; this is intended for local analysis rather than HPC execution.
- Optionally pass `--qc_config config/01_QC/<project>/<sample>.json` to add annotation labels from `new_labels` onto only the configured `cluster_col` resolution and export a label table.
- Keep it separate from `01_QC_xenium_spatial_clean_clustered.py` because clustree is R-based and has a distinct resolution-stability purpose.
- Start with cluster-resolution topology only; add QC-aware node overlays after the basic plot works.

Use clustree together with script 01 QC and script 03/04 marker dotplots when deciding whether to accept a clustering solution, apply QC masks and rerun script 00, tune BANKSY parameters, or consider re-segmentation.

## Next Steps

- Before running a Slurm array, check that `#SBATCH --array` matches the number of intended configs.
- Consider reducing leftover notebook-export inspection cells/comments in the active script when stabilizing it.
- Keep `new_labels` config handling in 01 QC only unless the workflow explicitly changes.
-Generate dot plots on a per sample basis across all available resolutions, this way we can assess maker expression across resoltions and comepare them back to clustree results.
  -The dot plot should annotate cell types where available using cell type annotations from qc 01 script configs.