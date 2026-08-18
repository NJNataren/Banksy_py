# Xenium BANKSY Workflow: Concise Agent Guide

Use this file as the quick orientation map for the Xenium BANKSY pipeline. Use `README_00_04_XENIUM_ANALYSIS_WORKFLOW.md` only when implementation details, commands, or path inventories are needed. Use `re-clustering.md` for the planned QC-filtered reclustering branch.

## Core Rule

Do not use `adata_spatial.X` from BANKSY spatial objects for biological marker-expression dotplots. Those matrices contain BANKSY-expanded/scaled neighbour features. Use clean expression AnnData objects with BANKSY cluster labels copied into `.obs`.

## Pipeline Map

0. `00_xenium_clustering_clean_adata.py`
   - Raw Xenium AnnData -> BANKSY clustering.
   - Saves clean log-normalized expression AnnData with BANKSY labels/embeddings.
   - Also writes PCA scree QC and pre-filter zero-count QC.

0a. `00a_plot_pca_scree_from_existing_adata.py`
   - Backfills PCA scree QC from existing AnnData without rerunning BANKSY.

1. `01_QC_xenium_spatial_clean_clustered.py`
   - QC/spatial review on script 00 clean clustered objects.
   - Produces filtering-decision masks/diagnostics for later reclustering.

1a. `01a_clustree_cluster_resolution_qc.R`
   - Local clustree resolution-stability QC from script 00 cluster CSVs.
   - `01a_merge_cluster_resolution_csvs.R` merges split per-resolution CSVs first when needed.

2. `02_create_expression_adata_with_banksy_clusters.py`
   - Legacy/helper path for rebuilding clean expression objects from existing BANKSY outputs or archived labels.

3. `03_export_dotplot_data_from_clean_script00_config.py`
   - Preferred export path for current BANKSY cluster dotplots from script 00 clean objects.

3-legacy. `03_export_dotplot_data_from_config.py`
   - Export path for script 02 or archived-label objects.

4. `04_plot_multi_sample_dotplot_from_config_local.py`
   - Local plotting from script 03 summary CSVs.

5. `05_apply_qc_filters_for_reclustering.py`
   - Applies reviewed script 01 QC filters.
   - Writes QC-annotated and QC-filtered clean objects for reclustering.

6. Planned: `06_recluster_qc_filtered_with_banksy.py`
   - Rerun the BANKSY portion of script 00 on script 05 filtered objects.
   - Include planned neighbour analysis. See `re-clustering.md`.

## Where Things Run

- HPC: scripts 00, 00a, 01, 02, 03, 05, and planned 06.
- Local: script 01a clustree QC and script 04 plotting after outputs are copied back.
- Slurm array ranges must match the number of JSON configs in the selected leaf config directory.
- Some wrappers accept `CONFIG_DIR` overrides, especially `run_00_xenium_clustering.sl` and `run_01_xenium_QC_array_ptmt.sl`.

## Config Roots

```text
config/00_clustering/
config/01_QC/
config/02_create_expression/
config/03_export_summary/
config/04_plot_dotplot/active/
config/05_apply_qc_filters/
config/06_recluster_qc_filtered/        planned
```

Useful project split:

```text
vbct          standard VBCT runs
ptmt          standard PTMT/PC35 runs
ptmt_pc55     PTMT PC55 comparison runs; keep separate from ptmt
```

## Key Objects

Raw input:

```text
data/xenium/raw_data/<raw_subdir>/<sample>_raw.h5ad
```

Preferred clean clustered object from script 00:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_with_banksy_clusters_<resolutions>.h5ad
```

QC-filtered reclustering input from script 05:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_qc_filtered_qc_v1.h5ad
```

Cluster assignment table used by clustree:

```text
data/xenium/processed/<project>/<sample>/<sample>_cell_cluster_id_across_clustering_res_<resolutions>.csv
```

Dotplot summary exports:

```text
data/xenium/processed/<project>/cross_sample_dotplot_exports/*_dotplot_summary.csv
```

Local dotplot figures:

```text
figures/dotplots/
```

## Current Preferences

- Prefer clean script 00 objects for marker analysis and dotplot exports.
- Prefer `03_export_dotplot_data_from_clean_script00_config.py` for current cluster dotplots.
- Use script 02 only for legacy/archived-label workflows or object reconstruction.
- Use script 00a when only PCA scree QC is missing.
- Keep PC55 outputs under `ptmt_pc55`; do not mix with standard `ptmt` outputs.
- Keep QC-filter labels in reclustering output filenames to avoid overwriting original script 00 outputs.

## Validation

Use targeted checks rather than full local workflow runs:

```bash
python -m py_compile <script.py>
python -m json.tool <config.json>
```

## Cautions

- Avoid hard-coded sample names or paths when configs define them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as large or machine-specific.
- Preserve clean expression values for biological summaries; BANKSY-expanded matrices are for clustering, not expression interpretation.
