# Xenium BANKSY QC, Clustering, and Dotplot Workflow: Scripts 00-04

This README documents the local/HPC workflow built around the `00`-`04` scripts in this repository. The workflow starts with Xenium expression data, runs sample QC, performs BANKSY spatial clustering, creates clean expression objects with BANKSY labels and optional archived labels, exports dotplot-ready summaries, and finally renders multi-sample dotplots.

The main biological use case is metastatic melanoma Xenium spatial transcriptomics analysis for immune checkpoint inhibitor response/resistance studies.

## Important Concept

Do not use `adata_spatial.X` from saved BANKSY spatial objects as marker-expression values for biological dotplots.

The BANKSY spatial objects contain BANKSY-expanded/scaled features in `.X`, including neighbour-derived features. Those values are appropriate for clustering, but not for marker-expression dotplots. The dotplot workflow therefore creates a clean expression AnnData object from the original Xenium expression matrix and copies grouping metadata into `.obs`.

Current supported grouping metadata include:

```text
BANKSY cluster labels from current spatial objects
Archived cell-type labels from older annotated AnnData objects
```

## Workflow Overview

0. `00_QC_testing_xenium_spatial.py` or `00_QC_testing_xenium_spatial_PTMT_v2.py`
   Run per-sample Xenium QC, inspection plots, and threshold summaries from raw AnnData files.

1. `01_xenium_clustering.py`
   Run BANKSY clustering from raw or processed Xenium AnnData files.

2. `02_create_expression_adata_with_banksy_clusters.py`
   Create clean expression AnnData files and copy BANKSY cluster labels plus optional archived cell-type labels into `.obs`.

3. `03_export_dotplot_data_from_config.py`
   Export long-format dotplot summary CSVs from the clean expression objects.

4. `04_plot_multi_sample_dotplot_from_config_local.py`
   Render multi-sample dotplots from the exported CSVs. The local script is the current plotting/tuning entrypoint.

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


## Where To Run Each Step

The current workflow is split between HPC computation and local plotting:

```text
HPC:   00_QC_testing_xenium_spatial.py / 00_QC_testing_xenium_spatial_PTMT_v2.py
HPC:   01_xenium_clustering.py
HPC:   02_create_expression_adata_with_banksy_clusters.py
HPC:   03_export_dotplot_data_from_config.py
Local: 04_plot_multi_sample_dotplot_from_config_local.py
```

Scripts 00-03 read or write full AnnData objects or large expression summaries, so they should normally be run on the HPC. Script 04 reads only script 03 CSV summaries and writes PNGs, so it is intended for local plotting and layout tuning after the relevant summary CSVs have been copied back from the HPC.

Typical full order:

```text
1. Run script 00 on the HPC for per-sample QC and inspection outputs.
2. Run script 01 on the HPC for BANKSY clustering.
3. Ensure archive_spatial_cell_type_labels.csv exists on the HPC when archived labels are needed.
4. Run script 02 on the HPC to regenerate clean AnnData objects with current and/or archived labels in .obs.
5. Run script 03 on the HPC to export dotplot summaries.
6. Copy the script 03 summary CSVs locally.
7. Run script 04 locally with the desired plot config.
```

For Slurm arrays, the array range must match the number of JSON config files in the selected leaf config directory. For eight sample configs, use:

```bash
--array=0-7
```


## Script 00: Xenium QC And Inspection

### Scripts

VBCT workflow:

```text
00_QC_testing_xenium_spatial.py
```

PTMT workflow:

```text
00_QC_testing_xenium_spatial_PTMT_v2.py
```

Notebook companions also exist:

```text
00_QC_testing_xenium_spatial.ipynb
00_QC_testing_xenium_spatial_PTMT_v2.ipynb
```

### Purpose

Runs per-sample Xenium QC from a JSON config. The script reads the raw Xenium AnnData object, creates sample-specific output directories, inspects count and gene distributions, generates spatial QC plots, applies transcript-count threshold annotations, writes QC tables, and creates exploratory plots that help decide whether the sample is suitable for downstream BANKSY clustering.

Script 00 is mainly an inspection/QC stage. It should be run before script 01 when starting from new raw Xenium samples or when QC thresholds and sample-level QC outputs need to be refreshed.

### Main Config Inputs

QC configs live under:

```text
config/00_QC/
```

Current config directories include:

```text
config/00_QC/vbct/*.json
config/00_QC/ptmt/*.json
```

Important config keys include:

```text
dataset_name
pc_label
lambda_label
res_label
nbr_weight_decay
coord_keys
new_labels
```

### Main Data Inputs

For a sample called `<dataset_name>`, script 00 expects raw Xenium AnnData here:

```text
data/xenium/raw_data/<dataset_name>_raw.h5ad
```

The object should include spatial coordinates and count/QC metadata in `.obs`, including commonly used columns such as:

```text
x
y
nCount_Xenium
nFeature_Xenium
```

### Main Outputs

Script 00 creates per-sample processed and output directories:

```text
data/xenium/processed/<dataset_name>/
data/xenium/output/<dataset_name>/
data/xenium/output/QC_testing/<dataset_name>/
```

Common QC outputs include:

```text
data/xenium/output/QC_testing/<dataset_name>/banksy_counts_and_genes_plot_<dataset_name>.png
data/xenium/output/QC_testing/<dataset_name>/transcript_knee_plot_<dataset_name>.png
data/xenium/output/QC_testing/<dataset_name>/tissue_spatial_scatter_min_counts_threshold_passed_cells_<dataset_name>.png
data/xenium/output/QC_testing/<dataset_name>/tissue_spatial_scatter_transcripts_qc_dynamic_range_<dataset_name>.png
data/xenium/output/QC_testing/<dataset_name>/obs_<dataset_name>.csv
data/xenium/output/QC_testing/<dataset_name>/gene_counts_across_all_cells_<dataset_name>_test.csv
```

The script also writes a float32 AnnData copy used by downstream workflows. Check the script output messages and paths when running on the HPC, because historical script variants have used slightly different raw/processed destinations for this file.

### How To Run On The HPC

For VBCT samples, use:

```text
run_00_xenium_QC_array_vbct.sl
```

This wrapper uses:

```text
CONFIG_DIR=config/00_QC/vbct
python 00_QC_testing_xenium_spatial.py --config <config>
```

Submit the full VBCT QC array with:

```bash
sbatch --array=0-7 run_00_xenium_QC_array_vbct.sl
```

For PTMT samples, use:

```text
run_00_xenium_QC_array_ptmt.sl
```

This wrapper uses:

```text
CONFIG_DIR=config/00_QC/ptmt
python 00_QC_testing_xenium_spatial_PTMT_v2.py --config <config>
```

Submit the full PTMT QC array with:

```bash
sbatch --array=0-7 run_00_xenium_QC_array_ptmt.sl
```

For a single config, run directly inside the `banksy` environment:

```bash
conda run -n banksy python 00_QC_testing_xenium_spatial.py \
  --config config/00_QC/vbct/CK_skin_res.json
```

As with the other Slurm wrappers, make sure the `#SBATCH --array` range matches the number of JSON files in the selected config directory.

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
config/01_clustering/
```

Current examples include:

```text
config/01_clustering/vbct/small/*.json
config/01_clustering/vbct/large/*.json
config/01_clustering/ptmt/*.json
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

### How To Run On The HPC

For a single sample, run script 01 directly inside the `banksy` environment:

```bash
conda run -n banksy python 01_xenium_clustering.py \
  --config config/01_clustering/vbct/small/CK_skin_res.json
```

For array jobs, use or adapt one of the Slurm wrappers:

```text
run_01_xenium_clustering_test.sl
run_01_xenium_clustering_test_large.sl
run_xenium_clustering_test_single.sl
```

Before submitting, check the wrapper's `CONFIG_DIR`, `#SBATCH --array`, and final `python ...` command. In this repository, some script 01 wrappers are testing-oriented and may call `01_xenium_clustering_clean_adata_test.py` rather than `01_xenium_clustering.py`.

Example Slurm submission pattern:

```bash
sbatch run_01_xenium_clustering_test.sl
```

If you change the config directory, make sure the array range matches the number of JSON configs in that directory.

## Script 02: Create Clean Expression AnnData With BANKSY Labels

### Script

```text
02_create_expression_adata_with_banksy_clusters.py
```

### Purpose

Creates a clean expression AnnData object for dotplot export. It starts from the original expression matrix, filters cells consistently with clustering, saves raw count-like values in a layer, normalizes/log-transforms expression, stores log-normalized expression in `.X` and `.raw`, and copies BANKSY cluster labels from one or more `adata_spatial_*.h5ad` files into `.obs`.

Script 02 can also attach archived cell-level labels exported from older annotated AnnData objects. In the current `vbct_small` configs, archived metadata are added as:

```text
archive_cluster_id
archive_cell_type_label
has_archive_label
```

Cells without an archived match are not dropped. They are assigned explicit missing-label values such as:

```text
archive_cluster_id = missing_archive_cluster
archive_cell_type_label = missing_archive_label
```

### Main Config Inputs

Configs live under:

```text
config/02_create_expression/
```

Current examples include:

```text
config/02_create_expression/vbct_small/*.json
config/02_create_expression/testing/*.json
```

Important config keys:

```text
expression_adata_path
output_h5ad
objects
filter_obs_key
filter_min
counts_layer
archive_labels
```

Each `objects` entry defines a BANKSY object and label column to copy, for example:

```json
{
  "adata_path": "data/xenium/processed/<sample>/adata_spatial_<sample>_0.70.h5ad",
  "groupby": "labels_scaled_gaussian_pc20_nc0.20_r0.70"
}
```

The optional `archive_labels` block defines a cell-level CSV to join by cell ID:

```json
"archive_labels": {
  "label_csv": "data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv",
  "sample": "BE_brain_non_res",
  "columns": {
    "cluster_id": "archive_cluster_id",
    "cell_type_label": "archive_cell_type_label"
  },
  "missing_values": {
    "archive_cluster_id": "missing_archive_cluster",
    "archive_cell_type_label": "missing_archive_label"
  },
  "has_label_column": "has_archive_label"
}
```

For samples whose archived labels use a different sample name, set `archive_labels.sample` explicitly. For example, the current `MF_skin_non_res` config reads archived labels from `MF_skin_non_res_roi`.

### Main Data Inputs

Usually:

```text
data/xenium/processed/<sample>/<sample>_float_32.h5ad
```

plus one or more BANKSY spatial objects from script 01:

```text
data/xenium/processed/<sample>/adata_spatial_<sample>_<resolution>.h5ad
```

When using archived labels, script 02 also needs the cell-level archived label table:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

This table can be created from archived annotated AnnData objects with:

```bash
conda run -n banksy python helper_scripts/archive_labels/export_archive_cell_type_labels.py \
  --config helper_scripts/archive_labels/config/archive_spatial_cluster_labels.json \
  --output-csv data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

### Main Outputs

Clean expression-plus-metadata AnnData files. The exact output path is controlled by `output_h5ad` in each config, usually under the sample processed directory, for example:

```text
data/xenium/processed/<sample>/<sample>_normalised_log1p_with_banksy_clusters.h5ad
```

These objects contain clean log-normalized expression values plus `.obs` grouping columns for downstream dotplots.

### How To Run On The HPC

Before running archived-label configs, make sure this file exists on the HPC:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

For one sample/config:

```bash
conda run -n banksy python 02_create_expression_adata_with_banksy_clusters.py \
  --config config/02_create_expression/vbct_small/CK_skin_res.json
```

For the full `vbct_small` config set, use the Slurm wrapper:

```bash
sbatch --array=0-7 \
  --export=CONFIG_DIR=config/02_create_expression/vbct_small \
  run_02_create_expression_adata_with_banksy_clusters.sl
```

Script 02 must be rerun after changing metadata transfer logic, because script 03 can only group by columns that already exist in the clean `.h5ad` objects.

## Script 03: Export Dotplot Summary CSVs

### Script

```text
03_export_dotplot_data_from_config.py
```

### Purpose

Exports long-format dotplot summaries from clean expression AnnData files created by script 02. Each output row corresponds to one sample/group/gene dot.

The grouping column is chosen by each object config's `groupby` key. Current common groupings are:

```text
labels_scaled_gaussian_pc*_nc*_r*    current BANKSY clusters
archive_cell_type_label              archived/manual cell-type labels
```

When grouping by archived labels, cells with no archived match are retained as the explicit `missing_archive_label` group.

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
config/03_export_summary/
```

Current production-style examples include:

```text
config/03_export_summary/vbct_small/all_genes/*.json
config/03_export_summary/vbct_small/canonical_markers/*.json
config/03_export_summary/vbct_small/all_genes_archive_labels/*.json
config/03_export_summary/vbct_small/canonical_markers_archive_labels/*.json
```

Testing examples include:

```text
config/03_export_summary/testing/export_summary/*.json
config/03_export_summary/testing/legacy/*.json
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
missing_group_label
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

For archived-label exports, each object uses the archived label column as its grouping key:

```json
{
  "sample": "BE_brain_non_res",
  "resolution": "archive",
  "adata_path": "data/xenium/processed/BE_brain_non_res/BE_brain_non_res_normalised_log1p_with_banksy_clusters.h5ad",
  "groupby": "archive_cell_type_label",
  "groupby_label": "archive_cell_type_label",
  "missing_group_label": "missing_archive_label"
}
```

### Main Data Inputs

Clean expression AnnData files from script 02, with BANKSY labels and optional archived labels in `.obs`.

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
group_id
group_label
sample_cluster
sample_group
groupby
groupby_label
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
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels/*_canonical_markers_archive_cell_type_labels_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels_all_genes/*_all_genes_archive_cell_type_labels_dotplot_summary.csv
```

### Example Runs

All genes:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/03_export_summary/vbct_small/all_genes/CK_skin_res_all_genes.json
```

Canonical markers:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/03_export_summary/vbct_small/canonical_markers/CK_skin_res_canonical_markers.json
```

Canonical markers grouped by archived labels:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/03_export_summary/vbct_small/canonical_markers_archive_labels/CK_skin_res_canonical_markers_archive_cell_type_labels.json
```

All genes grouped by archived labels:

```bash
conda run -n banksy python 03_export_dotplot_data_from_config.py \
  --config config/03_export_summary/vbct_small/all_genes_archive_labels/CK_skin_res_all_genes_archive_cell_type_labels.json
```

### How To Run On The HPC

Slurm wrapper:

```text
run_03_xenium_dotplot_export_from_config.sl
```

Run one of the leaf config directories below. The wrapper expects JSON files directly inside `CONFIG_DIR`, not one level above multiple config subdirectories.

Cluster-grouped all genes:

```bash
sbatch --array=0-7 \
  --export=CONFIG_DIR=config/03_export_summary/vbct_small/all_genes \
  run_03_xenium_dotplot_export_from_config.sl
```

Cluster-grouped canonical markers:

```bash
sbatch --array=0-7 \
  --export=CONFIG_DIR=config/03_export_summary/vbct_small/canonical_markers \
  run_03_xenium_dotplot_export_from_config.sl
```

Archived-label canonical markers:

```bash
sbatch --array=0-7 \
  --export=CONFIG_DIR=config/03_export_summary/vbct_small/canonical_markers_archive_labels \
  run_03_xenium_dotplot_export_from_config.sl
```

Archived-label all genes:

```bash
sbatch --array=0-7 \
  --export=CONFIG_DIR=config/03_export_summary/vbct_small/all_genes_archive_labels \
  run_03_xenium_dotplot_export_from_config.sl
```

Make sure the Slurm array range matches the number of JSON configs in the chosen directory. If there are eight configs, valid task IDs are `0-7`.

## Script 04: Render Multi-Sample Dotplots

### Scripts

Local working/tuning script:

```text
04_plot_multi_sample_dotplot_from_config_local.py
```

The local script currently has the most recent layout, review, z-scoring, and generic grouping features. It can plot either BANKSY-cluster summaries or archived-label summaries, depending on the config passed to `--config`.

### Purpose

Reads one or more script 03 dotplot summary CSVs and renders a combined multi-sample dotplot PNG. It can select one resolution per sample for cluster-based outputs, group/order genes, highlight selected genes, exclude reviewed genes, optionally z-score expression values for dot color, and label rows using generic sample/group metadata from script 03. One run produces one plot, determined by the JSON config.

### Main Config Inputs

Plot configs live under:

```text
config/04_plot_dotplot/
```

Current local configs:

```text
config/04_plot_dotplot/local/all_genes_multi_sample_local.json
config/04_plot_dotplot/local/canonical_markers_multi_sample_local.json
config/04_plot_dotplot/local/all_genes_archive_cell_type_labels_multi_sample_local.json
config/04_plot_dotplot/local/canonical_markers_archive_cell_type_labels_multi_sample_local.json
```

HPC examples:

```text
config/04_plot_dotplot/vbct_small/canonical_markers_multi_sample.json
config/04_plot_dotplot/vbct_large/canonical_markers_multi_sample.json
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
y_label_template
y_axis_label
figure
```

### Main Data Inputs

Script 03 summary CSVs, for example:

```text
data/xenium/processed/cross_sample_dotplot_exports/all_genes_all_res/*_all_genes_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/canonical_all_res/*_canonical_markers_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels/*_canonical_markers_archive_cell_type_labels_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels_all_genes/*_all_genes_archive_cell_type_labels_dotplot_summary.csv
```

For local plotting after HPC script 03 runs, copy the summary CSVs locally to the same relative paths expected by the plot config. You do not need to copy regenerated `.h5ad` files locally for script 04.

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

The script z-scores expression within each gene across the plotted sample/group rows. Dot color then represents relative expression for each gene, while dot size remains percent expressing.

Genes with zero variance after filtering receive a z-score of 0 for plotting.

### Current Local Outputs

All-gene cluster-grouped z-score dotplot:

```text
figures/dotplots/local_all_genes_tier1_canonical_highlight_zscore_dotplot.png
```

Canonical marker cluster-grouped z-score dotplot:

```text
figures/dotplots/local_canonical_markers_multi_sample_zscore_dotplot.png
```

All-gene archived-label z-score dotplot:

```text
figures/dotplots/local_all_genes_archive_cell_type_labels_tier1_canonical_highlight_zscore_dotplot.png
```

Canonical marker archived-label z-score dotplot:

```text
figures/dotplots/local_canonical_markers_archive_cell_type_labels_multi_sample_zscore_dotplot.png
```

These local configs currently use 150 DPI.

### How To Run Locally

Before running script 04 locally, copy the relevant script 03 summary CSVs from the HPC into the local path expected by the plot config. For archived-label all-gene plots, for example, copy the HPC files from:

```text
/scratchdata1/users/a1210419/Banksy_py/data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels_all_genes/
```

to the matching local repo path:

```text
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels_all_genes/
```

Then run exactly one plot config per command. One script 04 run produces one PNG.

All genes, cluster grouped:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/all_genes_multi_sample_local.json
```

Canonical markers, cluster grouped:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/canonical_markers_multi_sample_local.json
```

All genes grouped by archived labels:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/all_genes_archive_cell_type_labels_multi_sample_local.json
```

Canonical markers grouped by archived labels:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/canonical_markers_archive_cell_type_labels_multi_sample_local.json
```

### HPC Slurm Entrypoint

A Slurm wrapper exists:

```text
run_04_xenium_multi_sample_dotplot_from_config.sl
```

The current plotting/tuning workflow uses `04_plot_multi_sample_dotplot_from_config_local.py` locally after script 03 CSVs have been copied from the HPC. If using the Slurm wrapper, first confirm that the script it calls is present and aligned with the local plotting script.

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

### Script 00 Outputs Used Later

QC/inspection outputs are primarily reviewed by humans, but downstream scripts may use per-sample AnnData files and QC metadata generated during this stage:

```text
data/xenium/processed/<sample>/
data/xenium/output/QC_testing/<sample>/
```

### Script 01 Outputs Used Later

```text
data/xenium/processed/<sample>/<sample>_float_32.h5ad
data/xenium/processed/<sample>/adata_spatial_<sample>_<resolution>.h5ad
```

### Script 02 Outputs Used Later

Clean expression-plus-metadata `.h5ad` files defined by each `output_h5ad` in:

```text
config/02_create_expression/**/*.json
```

For archived-label workflows, these objects must contain:

```text
archive_cluster_id
archive_cell_type_label
has_archive_label
```

### Script 03 Outputs Used Later

```text
data/xenium/processed/cross_sample_dotplot_exports/all_genes_all_res/*_all_genes_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/canonical_all_res/*_canonical_markers_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels/*_canonical_markers_archive_cell_type_labels_dotplot_summary.csv
data/xenium/processed/cross_sample_dotplot_exports/archive_cell_type_labels_all_genes/*_all_genes_archive_cell_type_labels_dotplot_summary.csv
```

### Script 04 Outputs

```text
figures/dotplots/*.png
```

## Validation Commands

Syntax check:

```bash
python -m py_compile 00_QC_testing_xenium_spatial.py \
  00_QC_testing_xenium_spatial_PTMT_v2.py \
  01_xenium_clustering.py \
  02_create_expression_adata_with_banksy_clusters.py \
  03_export_dotplot_data_from_config.py \
  04_plot_multi_sample_dotplot_from_config_local.py
```

Validate one JSON config:

```bash
python -m json.tool config/04_plot_dotplot/local/all_genes_multi_sample_local.json
```

Render current local plots:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/all_genes_multi_sample_local.json

conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/canonical_markers_multi_sample_local.json

conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/all_genes_archive_cell_type_labels_multi_sample_local.json

conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/local/canonical_markers_archive_cell_type_labels_multi_sample_local.json
```

## Cautions

- Do not use BANKSY-expanded `adata_spatial.X` values for biological marker dotplots.
- Avoid hard-coding sample names or paths when configs already define them.
- Full script 00/01 and multi-sample script 02/03 runs can be compute-heavy; run locally only with a clear sample/test config.
- Archived-label script 02 runs require `data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv` to exist on the machine running script 02.
- Files under `data/xenium/`, `figures/`, `logs/`, and `hpc/` may be large or machine-specific.
- Review CSVs currently live under `data/xenium/raw_data/gene_markers/`; check whether they are git-ignored before relying on git to track them.
- Slurm array ranges must match the number of configs in the selected config directory.
