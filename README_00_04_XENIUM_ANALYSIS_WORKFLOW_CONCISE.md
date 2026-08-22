# Xenium BANKSY Workflow: Concise Agent Guide

Use this file as the quick orientation map for the Xenium BANKSY pipeline. Use `README_00_04_XENIUM_ANALYSIS_WORKFLOW.md` when implementation details, commands, or path inventories are needed. Use `re-clustering.md` only as historical notes for the QC-filtered reclustering branch.

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
   - Local clustree resolution-stability QC from script 00 cluster CSVs or script 06 recluster cluster CSVs.
   - `01a_merge_cluster_resolution_csvs.R` merges split per-resolution CSVs first when needed.
   - `01a_run_clustree_from_config.py` runs clustree batches from JSON configs, including post-07 VBCT small recluster plots.

2. `02_create_expression_adata_with_banksy_clusters.py`
   - Legacy/helper path for rebuilding clean expression objects from existing BANKSY outputs or archived labels.

3. `03_export_dotplot_data_from_clean_script00_config.py`
   - Preferred export path for current BANKSY cluster dotplots from script 00 clean objects.
   - Also works for script 06 clean recluster objects when configs explicitly list the recluster `.obs` label columns.

3-legacy. `03_export_dotplot_data_from_config.py`
   - Export path for script 02 or archived-label objects.

4. `04_plot_multi_sample_dotplot_from_config_local.py`
   - Local plotting from script 03 summary CSVs.

5. `05_apply_qc_filters_for_reclustering.py`
   - Applies reviewed script 01 QC filters as masks only.
   - Retains all cells and writes a QC-annotated clean object for provenance.
   - Adds `qc_keep_for_reclustering`, `qc_filter_status`, per-filter fail masks, and fail-reason columns.

6. `06_recluster_qc_annotated_with_banksy.py`
   - Reruns the BANKSY portion of script 00 from a script 05 QC-annotated object.
   - Default `recluster_inclusion` is `qc_pass_only`; optional `all_cells` sensitivity mode is supported.
   - Saves `results_df`, full-cell cluster tables, clean expression objects with recluster labels/embeddings, marker tables/plots, BANKSY plots, and cluster-count QC plots. Heavy per-resolution BANKSY spatial objects and `banksy_dict` are opt-in only.
   - Writes clean-object UMAP cluster plots under `umap_cluster_plot/` and cluster-plus-QC UMAPs under `umap_qc/`.

7. `07_squidpy_recluster_spatial_analysis.py`
   - Runs downstream Squidpy spatial interpretation on the clean expression object used by script 06.
   - Includes neighbourhood enrichment, centrality, co-occurrence, Moran's I, and top Moran gene spatial plots.

## Where Things Run

- HPC: scripts 00, 00a, 01, 02, 03, 05, 06, and 07.
- Local: script 01a clustree QC and script 04 plotting after outputs are copied back.
- Slurm array ranges must match the number of JSON configs in the selected leaf config directory.
- Some wrappers accept `CONFIG_DIR` overrides, especially `hpc/slurm/run_00_xenium_clustering.sl` and `hpc/slurm/run_01_xenium_QC_array_ptmt.sl`.

## Config Roots

```text
config/00_clustering/
config/01_QC/
config/01a_clustree/
config/02_create_expression/
config/03_export_summary/
config/03_export_summary/active/vbct_recluster_qc_pass_only/
config/04_plot_dotplot/active/
config/04_plot_dotplot/active/vbct_recluster_qc_pass_only/
config/05_apply_qc_filters/
config/06_recluster_qc_annotated/
config/07_squidpy_recluster_analysis/
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

QC-annotated reclustering input from script 05:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_qc_annotated_<input_label>.h5ad
```

Script 06 full provenance object:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_qc_annotated_<input_label>_with_banksy_reclusters_<run_label>_<resolutions>.h5ad
```

Script 06 clean object used for BANKSY and script 07:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_recluster_<run_label>_cells_used_for_banksy_with_clusters_<resolutions>.h5ad
```

Original script 00 cluster assignment table used by clustree:

```text
data/xenium/processed/<project>/<sample>/<sample>_cell_cluster_id_across_clustering_res_<resolutions>.csv
```

Script 06 recluster cluster assignment table:

```text
data/xenium/processed/<project>/<sample>/<sample>_recluster_<run_label>_cell_cluster_id_across_clustering_res_<resolutions>.csv
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

- Prefer clean script 00 objects for original marker analysis and dotplot exports.
- For marker-expression dotplots after QC reclustering, use the script 06 clean recluster object containing only cells used by BANKSY: `adata_expression_clean_<sample>_recluster_<run_label>_cells_used_for_banksy_with_clusters_<resolutions>.h5ad`.
- Prefer `03_export_dotplot_data_from_clean_script00_config.py` for current cluster dotplots.
- Use script 02 only for legacy/archived-label workflows or object reconstruction.
- Use script 00a when only PCA scree QC is missing.
- Keep PC55 outputs under `ptmt_pc55`; do not mix with standard `ptmt` outputs.
- Keep `input_label` and `run_label` in reclustering output filenames to avoid overwriting original script 00 outputs.
- Default script 06 to `recluster_inclusion = "qc_pass_only"`; use `all_cells` only as a sensitivity analysis.
- Keep Squidpy analyses in script 07 so BANKSY does not need to be rerun when spatial interpretation settings change.
- For post-07 recluster clustree plots, use `./helper_scripts/local_runs/run_01a_clustree_from_config_local.sh` or `python3 01a_run_clustree_from_config.py --config config/01a_clustree/vbct/small_recluster_qc_pass_only.json`.

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

## Current CK Smoke-Test State

- CK smoke testing used project `vbct`, sample `CK_skin_res`, and `run_label = filtered_qc_v1_qc_pass_only_smoke`.
- Script 05 kept 6,387 / 6,760 cells for reclustering; 373 cells were excluded by QC masks.
- Script 06 completed qc-pass-only reclustering for resolutions `0.50` through `1.50`.
- Script 07 completed on the clean recluster object for analysis resolution `1.50`.
- Script 03/04 marker dotplots were generated using configs under `vbct_recluster_qc_pass_only`.
- The all-gene script 03 export contains all 304 genes. Script 04 now has configs for canonical markers, reviewed/manual-filtered all-gene plotting, and all-304-gene plotting.

## Current Clustree State

- `01a_clustree_cluster_resolution_qc.R` supports script 06 recluster CSVs via `--cluster_suffix`, for example `_recluster_filtered_qc_v1_qc_pass_only` or `_recluster_filtered_qc_v1_qc_pass_only_smoke`.
- Config-driven VBCT small post-07 clustree config:

```text
config/01a_clustree/vbct/small_recluster_qc_pass_only.json
```

- Local helper command:

```bash
./helper_scripts/local_runs/run_01a_clustree_from_config_local.sh
```

- Dry-run preview command:

```bash
./helper_scripts/local_runs/run_01a_clustree_from_config_local.sh --dry-run
```

- The runner reuses script 06 configs from `config/06_recluster_qc_annotated/vbct/small` and resolves `<sample>_recluster_<run_label>_cell_cluster_id_across_clustering_res_<resolutions>.csv` automatically.
- Skip `--qc_config` for post-recluster clustree unless annotation configs explicitly refer to recluster cluster columns.

