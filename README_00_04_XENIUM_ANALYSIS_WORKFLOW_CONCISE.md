# Xenium BANKSY 00-04 Workflow: Concise Guide

This is the short working guide for the Xenium BANKSY 00-04 workflow. For exact historical commands, legacy variants, and detailed path inventories, see `README_00_04_XENIUM_ANALYSIS_WORKFLOW_V1.md`.

For agent context, prefer `context/XENIUM_00_04_AGENT_CONTEXT.md` to avoid flooding context windows.

## Core Principle

Do not use `adata_spatial.X` from BANKSY spatial objects for biological marker-expression dotplots. Those objects contain BANKSY-expanded/scaled features, including neighbour-derived values. Use clean expression AnnData objects with BANKSY cluster labels copied into `.obs`.

## Workflow

0. `00_xenium_clustering_clean_adata.py`
   - BANKSY clustering from raw Xenium AnnData.
   - Saves clean expression AnnData with BANKSY labels and embeddings.

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

- HPC: scripts 00-03.
- Local: script 04 after copying summary CSVs from the HPC.
- Slurm array ranges must match the number of JSON configs in the selected leaf config directory.

## Important Config Roots

```text
config/01_clustering/
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

Clean script 00 dotplot exports:

```text
data/xenium/processed/vbct/cross_sample_dotplot_exports/*_dotplot_summary.csv
data/xenium/processed/ptmt/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv
```

Local plot outputs:

```text
figures/dotplots/
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
- Use the full V1 README only when its extra detail is actually needed.

