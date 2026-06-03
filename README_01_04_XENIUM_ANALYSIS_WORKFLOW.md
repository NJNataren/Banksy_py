# Xenium BANKSY Clustering and Dotplot Workflow: Scripts 01-04

This README documents the local/HPC workflow built around the `01`-`04` scripts in this repository. The workflow starts with Xenium expression data, performs BANKSY spatial clustering, creates clean expression objects with BANKSY labels, exports dotplot-ready summaries, and finally renders multi-sample dotplots.

The main biological use case is metastatic melanoma Xenium spatial transcriptomics analysis for immune checkpoint inhibitor response/resistance studies.

## Important Concept

Do not use `adata_spatial.X` from saved BANKSY spatial objects as marker-expression values for biological dotplots.

The BANKSY spatial objects contain BANKSY-expanded/scaled features in `.X`, including neighbour-derived features. Those values are appropriate for clustering, but not for marker-expression dotplots. The dotplot workflow therefore creates a clean expression AnnData object from the original Xenium expression matrix and copies BANKSY cluster labels into `.obs`.

## Workflow Overview

1. `01_xenium_clustering.py`
   Run BANKSY clustering from raw Xenium AnnData files.

2. `02_create_expression_adata_with_banksy_clusters.py`
   Create clean expression AnnData files and copy BANKSY cluster labels into `.obs`.

3. `03_export_dotplot_data_from_config.py`
   Export long-format dotplot summary CSVs from the clean expression objects.

4. `04_plot_multi_sample_dotplot_from_config.py` or `04_plot_multi_sample_dotplot_from_config_local.py`
   Render multi-sample dotplots from the exported CSVs.

## Environment

The recent local/HPC workflow uses the conda environment:

```bash
conda activate banksy
```

or for one-off local runs:

```bash
conda run -n banksy python <script.py> --config <config.json>
```

Core dependencies include `scanpy`, `anndata`, `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `igraph`, `leidenalg`, and `umap-learn`.

## Script 01: BANKSY Clustering

### Script

```text
01_xenium_clustering.py
```

Notebook/source companion:

```text
01_xenium_clustering.ipynb
```

A related test/clean-adata notebook/script also exists:

```text
01_xenium_clustering_clean_adata_test.py
01_xenium_clustering_clean_adata_test.ipynb
```

### Purpose

Runs BANKSY spatial clustering for a Xenium sample. It reads a sample-specific JSON config, loads the raw AnnData object, filters zero-count cells, downcasts to float32, normalizes/log-transforms expression for clustering, builds the BANKSY spatial graph/matrix, runs PCA/UMAP, performs Leiden clustering across configured resolutions, plots clustering results, and saves clustering outputs.

### Main Config Inputs

Clustering configs live under:

```text
config/clustering/
```

Current examples include:

```text
config/clustering/vbct/small/*.json
config/clustering/vbct/large/*.json
config/clustering/ptmt/*.json
```

Important config keys used by `01_xenium_clustering.py` include:

```text
dataset_name
pc_label
lambda_label
res_label
nbr_weight_decay
coord_keys
```

### Main Data Inputs

For a sample called `<dataset_name>`, the script expects raw Xenium AnnData here:

```text
data/xenium/raw_data/<dataset_name>_raw.h5ad
```

The AnnData object should include coordinate columns in `.obs`, usually:

```text
x
y
```

and the count/QC column used for filtering:

```text
nCount_Xenium
```

### Main Outputs

Per-sample outputs are written under:

```text
data/xenium/processed/<dataset_name>/
data/xenium/output/<dataset_name>/
```

Important outputs include:

```text
data/xenium/processed/<dataset_name>/<dataset_name>_float_32.h5ad
data/xenium/processed/<dataset_name>/adata_spatial_<dataset_name>_<resolution>.h5ad
data/xenium/processed/<dataset_name>/<dataset_name>_pc<pc_label>_nc<lambda_label>_r<resolution>_banksy_dict.pkl.gz
data/xenium/processed/<dataset_name>/<dataset_name>_cell_cluster_id_across_clustering_res_<resolution>.csv
```

The `adata_spatial_*.h5ad` files are needed by script 02 for BANKSY labels, but their `.X` values should not be used directly for marker-expression dotplots.

### Example Run

```bash
conda run -n banksy python 01_xenium_clustering.py \
  --config config/clustering/vbct/small/CK_skin_res.json
```

## Script 02: Create Clean Expression AnnData With BANKSY Labels

### Script

```text
02_create_expression_adata_with_banksy_clusters.py
```

### Purpose

Creates a clean expression AnnData object for dotplot export. It starts from the original expression matrix, filters cells consistently with clustering, saves raw count-like values in a layer, normalizes/log-transforms expression, stores log-normalized expression in `.X` and `.raw`, and copies BANKSY cluster labels from one or more `adata_spatial_*.h5ad` files into `.obs`.

### Main Config Inputs

Configs live under:

```text
config/dotplot/create_expression/
```

Current examples include:

```text
config/dotplot/create_expression/vbct_small/*.json
config/dotplot/testing/create_expression/*.json
```

Important config keys:

```text
expression_adata_path
output_h5ad
objects
filter_obs_key
filter_min
counts_layer
```

Each `objects` entry defines a BANKSY object and label column to copy, for example:

```json
{
  "adata_path": "data/xenium/processed/<sample>/adata_spatial_<sample>_0.70.h5ad",
  "groupby": "banksy_cluster_pc20_nc0.20_r0.70"
}
```

### Main Data Inputs

Usually:

```text
data/xenium/processed/<sample>/<sample>_float_32.h5ad
```

plus one or more BANKSY spatial objects from script 01:

```text
data/xenium/processed/<sample>/adata_spatial_<sample>_<resolution>.h5ad
```

### Main Outputs

Clean expression-plus-cluster AnnData files, usually under a processed dotplot directory such as:

```text
data/xenium/processed/cross_sample_dotplot_exports/
```

The exact output path is controlled by `output_h5ad` in each config.

### Example Run

```bash
conda run -n banksy python 02_create_expression_adata_with_banksy_clusters.py \
  --config config/dotplot/create_expression/vbct_small/CK_skin_res.json
```

## Script 03: Export Dotplot Summary CSVs

### Script

```text
03_export_dotplot_data_from_config.py
```

### Purpose

Exports long-format dotplot summaries from clean expression AnnData files created by script 02. Each output row corresponds to one sample/resolution/cluster/gene dot.

Dot color is based on:

```text
mean_expression
```

Dot size is based on:

```text
percent_expressing
```

### Main Config Inputs

Configs live under:

```text
config/dotplot/export_summary/
```

Current production-style examples include:

```text
config/dotplot/export_summary/vbct_small/all_genes/*.json
config/dotplot/export_summary/vbct_small/canonical_markers/*.json
```

Testing examples include:

```text
config/dotplot/testing/export_summary/*.json
config/dotplot/testing/legacy/*.json
```

Important config keys:

```text
objects
output_csv
expression_source
use_all_genes
marker_file
gene_column
marker_group_column
split_by
output_csv_template
```

For all-gene exports:

```json
"use_all_genes": true
```

For marker-list exports:

```json
"marker_file": "data/xenium/raw_data/gene_markers/xenium_melanoma_canonical_genes_May_2026.csv",
"gene_column": "Gene",
"marker_group_column": "Tier_1_annotation"
```

### Main Data Inputs

Clean expression AnnData files from script 02, with BANKSY labels in `.obs`.

Marker files when not using all-gene mode, for example:

```text
data/xenium/raw_data/gene_markers/xenium_melanoma_canonical_genes_May_2026.csv
data/xenium/raw_data/gene_markers/xenium_gene_list_annotation_master_v1_May_2026_Tier1_v3.csv
```

### Main Outputs

Long-format CSVs with columns including:

```text
sample
resolution
cluster_id
sample_cluster
groupby
gene
marker_group
mean_expression
percent_expressing
n_cells
expression_source
adata_path
```

Current local plotting inputs include:

```text
data/xenium/processed/cross_sample_dotplot_exports/all_genes_all_res/*_all_genes_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/canonical_all_res/*_canonical_markers_dotplot_summary.csv
```

### Example Runs

All genes:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/dotplot/export_summary/vbct_small/all_genes/CK_skin_res_all_genes.json
```

Canonical markers:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/dotplot/export_summary/vbct_small/canonical_markers/CK_skin_res_canonical_markers.json
```

### HPC Slurm Entrypoint

```text
run_03_xenium_dotplot_export_from_config.sl
```

Example:

```bash
sbatch --export=CONFIG_DIR=config/dotplot/export_summary/vbct_small/all_genes \
  run_03_xenium_dotplot_export_from_config.sl
```

Make sure the Slurm array range matches the number of JSON configs in the chosen directory.

## Script 04: Render Multi-Sample Dotplots

### Scripts

HPC-oriented script:

```text
04_plot_multi_sample_dotplot_from_config.py
```

Local working/tuning script:

```text
04_plot_multi_sample_dotplot_from_config_local.py
```

The local script currently has the most recent layout and review features.

### Purpose

Reads one or more script 03 dotplot summary CSVs and renders a combined multi-sample dotplot PNG. It can select one resolution per sample, group/order genes, highlight selected genes, exclude reviewed genes, and optionally z-score expression values for dot color.

### Main Config Inputs

Plot configs live under:

```text
config/dotplot/plot_dotplot/
```

Current local configs:

```text
config/dotplot/plot_dotplot/local/all_genes_multi_sample_local.json
config/dotplot/plot_dotplot/local/canonical_markers_multi_sample_local.json
```

HPC examples:

```text
config/dotplot/plot_dotplot/vbct_small/canonical_markers_multi_sample.json
config/dotplot/plot_dotplot/vbct_large/canonical_markers_multi_sample.json
```

Important config keys:

```text
output_png
title
dotplot_summary_csvs
sample_order
sample_filters
gene_file
gene_column
gene_group_column
highlight_gene_file
highlight_gene_column
gene_review_file
gene_review_gene_column
gene_review_keep_column
z_score_expression
z_score_expression_column
z_score_clip
color_value_column
colorbar_label
figure
```

### Main Data Inputs

Script 03 summary CSVs, for example:

```text
data/xenium/processed/cross_sample_dotplot_exports/all_genes_all_res/*_all_genes_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/canonical_all_res/*_canonical_markers_dotplot_summary.csv
```

Gene annotation/grouping files:

```text
data/xenium/raw_data/gene_markers/xenium_gene_list_annotation_master_v1_May_2026_Tier1_v3.csv
data/xenium/raw_data/gene_markers/xenium_melanoma_canonical_genes_May_2026.csv
```

Review/exclusion CSVs:

```text
data/xenium/raw_data/gene_markers/xenium_gene_review_for_dotplot.csv
data/xenium/raw_data/gene_markers/xenium_canonical_gene_review_for_dotplot.csv
```

### Gene Review Files

The review CSVs allow collaborators to remove genes without editing the plotting script.

Set `keep_for_dotplot` to one of these values to exclude a gene:

```text
FALSE
false
no
0
drop
exclude
```

Optional fields:

```text
drop_reason
reviewer_notes
```

The local script applies gene exclusions before z-scoring, so z-scores are calculated only on retained genes.

### Z-Scoring

When enabled:

```json
"z_score_expression": true,
"z_score_expression_column": "mean_expression",
"z_score_clip": 2.5,
"color_value_column": "mean_expression_zscore",
"colorbar_label": "Mean expression z-score"
```

The script z-scores expression within each gene across the plotted sample/cluster rows. Dot color then represents relative expression for each gene, while dot size remains percent expressing.

Genes with zero variance after filtering receive a z-score of 0 for plotting.

### Current Local Outputs

All-gene z-score dotplot:

```text
figures/dotplots/local_all_genes_tier1_canonical_highlight_zscore_dotplot.png
```

Canonical marker z-score dotplot:

```text
figures/dotplots/local_canonical_markers_multi_sample_zscore_dotplot.png
```

Both local configs currently use 150 DPI.

### Example Local Runs

All genes:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/dotplot/plot_dotplot/local/all_genes_multi_sample_local.json
```

Canonical markers:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/dotplot/plot_dotplot/local/canonical_markers_multi_sample_local.json
```

### HPC Slurm Entrypoint

```text
run_04_xenium_multi_sample_dotplot_from_config.sl
```

Example:

```bash
sbatch --export=CONFIG_DIR=config/dotplot/plot_dotplot/vbct_small \
  run_04_xenium_multi_sample_dotplot_from_config.sl
```

## Current Local Gene Review Notes

Low-expression review decisions from the canonical review file have been copied into the all-gene review file for matching canonical genes.

Currently copied excluded genes include:

```text
CLEC9A
FCER1A
LGALS2
CLDN5
SELE
SOX17
AQP3
IFNG
```

Strict low-expression candidates identified from the all-gene dotplot summaries were:

```text
CD1B
PLAC9
```

Using a looser max-percent-expressing cutoff also flagged:

```text
LINC00636
```

These should be reviewed biologically before final exclusion.

## File Dependency Summary

### Raw/Input Data

```text
data/xenium/raw_data/<sample>_raw.h5ad
data/xenium/raw_data/gene_markers/*.csv
```

### Script 01 Outputs Used Later

```text
data/xenium/processed/<sample>/<sample>_float_32.h5ad
data/xenium/processed/<sample>/adata_spatial_<sample>_<resolution>.h5ad
```

### Script 02 Outputs Used Later

Clean expression-plus-cluster `.h5ad` files defined by each `output_h5ad` in:

```text
config/dotplot/create_expression/**/*.json
```

### Script 03 Outputs Used Later

```text
data/xenium/processed/cross_sample_dotplot_exports/all_genes_all_res/*_all_genes_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/canonical_all_res/*_canonical_markers_dotplot_summary.csv
```

### Script 04 Outputs

```text
figures/dotplots/*.png
```

## Validation Commands

Syntax check:

```bash
python -m py_compile 02_create_expression_adata_with_banksy_clusters.py \
  03_export_dotplot_data_from_config.py \
  04_plot_multi_sample_dotplot_from_config_local.py
```

Validate one JSON config:

```bash
python -m json.tool config/dotplot/plot_dotplot/local/all_genes_multi_sample_local.json
```

Render current local plots:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/dotplot/plot_dotplot/local/all_genes_multi_sample_local.json

conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/dotplot/plot_dotplot/local/canonical_markers_multi_sample_local.json
```

## Cautions

- Do not use BANKSY-expanded `adata_spatial.X` values for biological marker dotplots.
- Avoid hard-coding sample names or paths when configs already define them.
- Full script 01 and multi-sample script 02/03 runs can be compute-heavy; run locally only with a clear sample/test config.
- Files under `data/xenium/`, `figures/`, `logs/`, and `hpc/` may be large or machine-specific.
- Review CSVs currently live under `data/xenium/raw_data/gene_markers/`; check whether they are git-ignored before relying on git to track them.
- Slurm array ranges must match the number of configs in the selected config directory.
