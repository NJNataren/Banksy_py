# Xenium BANKSY 00-04 Workflow: Concise Guide

This is the short working guide for the Xenium BANKSY 00-04 workflow. For exact historical commands, legacy variants, and detailed path inventories, see `README_00_04_XENIUM_ANALYSIS_WORKFLOW_V1.md`.

For agent context, prefer `context/XENIUM_00_04_AGENT_CONTEXT.md` to avoid flooding context windows.

## Core Principle

Do not use `adata_spatial.X` from BANKSY spatial objects for biological marker-expression dotplots. Those objects contain BANKSY-expanded/scaled features, including neighbour-derived values. Use clean expression AnnData objects with BANKSY cluster labels copied into `.obs`.

## Workflow

0. `00_xenium_clustering_clean_adata.py`
   - BANKSY clustering from raw Xenium AnnData.
   - Saves clean expression AnnData with BANKSY labels and embeddings.
   - Writes PCA scree QC outputs under each sample output directory.

0a. `00a_plot_pca_scree_from_existing_adata.py`
   - Backfills PCA scree plots from existing processed AnnData without rerunning BANKSY.
   - Use `run_00a_pca_scree_array.sl` for queued HPC array runs.

1a-prep. `01a_merge_cluster_resolution_csvs.R`
   - Local helper for large samples whose script 00 cluster assignment CSVs were generated one resolution at a time.
   - Joins split one-resolution CSVs by cell ID into a single clustree-ready table.

1a. `01a_clustree_cluster_resolution_qc.R`
   - Local R-based clustree resolution-stability QC from script 00 cluster assignment CSVs.
   - Optional but recommended before choosing annotation resolutions.
   - Set up local R dependencies with `Rscript setup_clustree_r_env.R`.

1. `01_QC_xenium_spatial_clean_clustered.py`
   - QC and spatial inspection from the clean clustered object made by script 00.

2. `02_create_expression_adata_with_banksy_clusters.py`
   - Legacy/archived-label path.
   - Use only when regenerating clean expression objects from existing BANKSY outputs or attaching archived labels.

3. `03_export_dotplot_data_from_clean_script00_config.py`
   - Preferred current-cluster dotplot export path.
   - Reads script 00 clean clustered objects directly.

4. `03_export_dotplot_data_from_config.py`
   - Script 02 / archived-label export path.

5. `04_plot_multi_sample_dotplot_from_config_local.py`
   - Local plotting from script 03 summary CSVs.

## Execution Split

- HPC: scripts 00-03, plus optional script 00a for scree QC backfill.
- Local: script 01a for clustree resolution QC and script 04 after copying summary CSVs from the HPC.
- Slurm array ranges should match the number of JSON configs in the selected leaf config directory; `run_00a_pca_scree_array.sl` exits cleanly for array IDs beyond the config count.

## Important Config Roots

```text
config/00_clustering/
config/02_create_expression/
config/03_export_summary/
config/04_plot_dotplot/local/
```

Preferred clean script 00 export configs:

```text
config/03_export_summary/vbct_clean_script00/
config/03_export_summary/ptmt_clean_script00/
```

## Key Paths

Raw input:

```text
data/xenium/raw_data/<raw_subdir>/<dataset_name>_raw.h5ad
data/xenium/raw_data/<dataset_name>_raw.h5ad
```

Preferred clean clustered object:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_with_banksy_clusters_<resolutions>.h5ad
```

PCA scree QC outputs from script 00 or 00a:

```text
data/xenium/output/<project>/<sample>/pca_qc/<sample>_pca_scree_plot.png
data/xenium/output/<project>/<sample>/pca_qc/<sample>_pca_scree_variance.csv
```

Merged split-resolution CSV outputs from script 01a-prep:

```text
data/xenium/processed/<project>/<sample>/<sample>_cell_cluster_id_across_clustering_res_<resolutions>.csv
data/xenium/processed/<project>/<sample>/<sample>_cell_cluster_id_across_clustering_res_<resolutions>_merge_summary.csv
```

Clustree resolution QC outputs from script 01a:

```text
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc.png
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc.pdf
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc_sc3_stability.png
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc_sc3_stability.pdf
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc_annotated.png
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_resolution_qc_annotated.pdf
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_input_columns.csv
data/xenium/output/<project>/<sample>/clustree_qc/<sample>_clustree_annotation_labels.csv
```

Clean script 00 dotplot exports:

```text
data/xenium/processed/vbct/cross_sample_dotplot_exports/*_dotplot_summary.csv
data/xenium/processed/ptmt/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv
```

Local plot outputs:

```text
figures/dotplots/
```

## Local R Setup for Clustree

Use a project-local `renv` environment for the R-based clustree QC step:

```bash
Rscript setup_clustree_r_env.R
```

If `renv.lock` already exists, this restores the recorded package versions. Otherwise it initializes `renv`, installs the required clustree packages, and writes `renv.lock`.

Example merge for split large-sample resolution CSVs:

```bash
Rscript 01a_merge_cluster_resolution_csvs.R \
  --input_dir data/xenium/processed/vbct/CK_bowel_res \
  --dataset_name CK_bowel_res \
  --cluster_prefix labels_scaled_gaussian_pc30_nc0.20_r
```

The merge helper consumes one-resolution CSVs, skips existing merged or summary CSVs, checks that cell IDs match across files, sorts resolution columns numerically, and writes a merged CSV named from the available resolutions. Rerun it after new 0.90 or 1.00 files are generated.

Example CK clustree local run:

```bash
Rscript 01a_clustree_cluster_resolution_qc.R \
  --cluster_csv data/xenium/processed/vbct/CK_skin_res/CK_skin_res_cell_cluster_id_across_clustering_res_0.50_0.60_0.70_0.80_0.90_1.00_1.10.csv \
  --dataset_name CK_skin_res \
  --cluster_prefix labels_scaled_gaussian_pc30_nc0.20_r \
  --output_dir data/xenium/output/vbct/CK_skin_res/clustree_qc \
  --qc_config config/01_QC/vbct/CK_skin_res.json
```

## Archived Labels

Archived-label workflows require:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

Script 02 can attach `archive_cluster_id`, `archive_cell_type_label`, and `has_archive_label`. Cells without archived matches are retained with explicit missing-label values.

## Validation

Use targeted checks rather than full local workflow runs:

```bash
python -m py_compile <script.py>
python -m json.tool <config.json>
```

## Cautions

- Avoid hard-coded sample names or paths when configs define them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as large or machine-specific.
- Prefer `03_export_dotplot_data_from_clean_script00_config.py` for current BANKSY cluster dotplots when script 00 clean objects already exist.
- Use `00a_plot_pca_scree_from_existing_adata.py` instead of rerunning script 00 when only scree plots are missing.
- Use the full workflow README only when its extra detail is actually needed.

