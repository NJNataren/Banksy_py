# Sankey And Archived Label Mapping Helpers

This directory contains helper scripts for exporting old archived cell labels and
mapping those labels onto the current clean dotplot workflow. The scripts are
not part of the main `00`-`04` analysis sequence, but they support the archived
label workflow used by scripts `02`, `03`, and `04`.

## Purpose

The archived AnnData objects contain older manual or biological cluster labels.
These labels are useful for two related tasks:

```text
1. Transfer old cell-type labels onto current clean AnnData objects.
2. Compare old annotated clusters with newer BANKSY clusters by cell overlap.
```

The main outputs are cell-level CSVs under:

```text
data/xenium/processed/cluster_label_exports/
```

## Files

### `05_export_cluster_labels_from_adata_spatial.py`

Exports detailed cell-level cluster labels from archived or current
`adata_spatial*.h5ad` objects.

Config:

```text
helper_scripts/archive_labels/config/archive_spatial_cluster_labels.json
```

Inputs:

```text
data/xenium/output/archive/<sample>/<sample>_clustered_spatial_pc*_nc*_r*.h5ad
```

Outputs:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cluster_labels.csv
data/xenium/processed/cluster_label_exports/archive_spatial_cluster_label_summary.csv
```

The full cell-level output includes columns such as:

```text
cell_id
sample
resolution
pc_label
lambda_label
cluster_id
cluster_column
annotation_column
cluster_annotation
adata_path
x
y
cell type
nCount_Xenium
nFeature_Xenium
cell_area
```

`cluster_annotation` is the readable archived label when an annotation column is
available. The literal `.obs["cell type"]` column may be numeric in some older
objects, so use `cluster_annotation` when you want the old readable biological
label.

Run:

```bash
conda run -n banksy python sankey/05_export_cluster_labels_from_adata_spatial.py \
  --config helper_scripts/archive_labels/config/archive_spatial_cluster_labels.json
```

### `export_archive_cell_type_labels.py`

Exports a slimmer cell-level table for transferring archived labels into clean
expression AnnData objects. This is the helper used by script `02` via the
`archive_labels` config block.

Config:

```text
helper_scripts/archive_labels/config/archive_spatial_cluster_labels.json
```

Recommended output:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

Output columns include:

```text
cell_id
sample
resolution
pc_label
lambda_label
cluster_id
cell_type_label
cluster_column
cell_type_column
adata_path
```

`cell_type_label` is the old readable label copied from the archived AnnData
annotation column. Script `02_create_expression_adata_with_banksy_clusters.py`
joins this table to the clean expression AnnData by `cell_id` and writes:

```text
archive_cluster_id
archive_cell_type_label
has_archive_label
```

Run:

```bash
conda run -n banksy python helper_scripts/archive_labels/export_archive_cell_type_labels.py \
  --config helper_scripts/archive_labels/config/archive_spatial_cluster_labels.json \
  --output-csv data/xenium/processed/cluster_label_exports/archive_spatial_cell_type_labels.csv
```

### `06_map_old_cluster_labels_to_current_clusters.py`

Maps old archived labels onto current BANKSY clusters by cell overlap. This is
useful when you want to compare older manual annotations with current clustering
results, or create a guide that says which old label best matches each current
cluster.

Config:

```text
sankey/config/cluster_label_mapping/archive_to_current_dotplot_clusters.json
```

Inputs:

```text
data/xenium/processed/cluster_label_exports/archive_spatial_cluster_labels.csv
data/xenium/processed/<sample>/<sample>_normalised_log1p_with_banksy_clusters.h5ad
```

Main outputs:

```text
data/xenium/processed/cluster_label_exports/archive_to_current_dotplot_cluster_label_mapping.csv
data/xenium/processed/cluster_label_exports/archive_to_current_dotplot_cluster_label_mapping_summary.csv
data/xenium/processed/cluster_label_exports/archive_to_current_dotplot_cluster_label_guide.csv
data/xenium/processed/cluster_label_exports/archive_to_current_dotplot_cell_label_mapping.csv
```

The cluster-level mapping output reports, for each current cluster:

```text
current_cluster_id
old_cluster_id_majority
old_cluster_annotation_majority
new_cluster_n_cells
matched_old_n_cells
majority_old_fraction_of_new_cluster
majority_old_fraction_of_matched_cells
mapping_status
old_cluster_mix
```

Important config fields:

```text
old_cluster_labels_csv
current_adata_path
current_cluster_column
old_sample
old_labels_are_roi_subset
mapping_denominator
min_overlap_fraction
```

For ROI/subset archived labels, such as `MF_skin_non_res_roi`, the config can
set:

```json
"old_labels_are_roi_subset": true,
"mapping_denominator": "matched_cells"
```

This avoids penalizing a current cluster for cells that were never present in
the old ROI object.

Run:

```bash
conda run -n banksy python sankey/06_map_old_cluster_labels_to_current_clusters.py \
  --config sankey/config/cluster_label_mapping/archive_to_current_dotplot_clusters.json
```

## Recommended Run Order

For archived-label dotplots:

```text
1. Run export_archive_cell_type_labels.py to create archive_spatial_cell_type_labels.csv.
2. Run script 02 with configs in config/02_create_expression/ so clean AnnData
   objects receive archive_cell_type_label metadata.
3. Run script 03 with archive-label configs in config/03_export_summary/.
4. Copy script 03 summary CSVs locally and run script 04.
```

For old-to-current cluster comparison:

```text
1. Run 05_export_cluster_labels_from_adata_spatial.py.
2. Run script 02 if the current clean AnnData objects need to be regenerated.
3. Run 06_map_old_cluster_labels_to_current_clusters.py.
4. Review the cluster mapping and dotplot label guide CSVs.
```

## Relationship To The Main Workflow

The main numbered workflow uses configs under:

```text
config/00_QC/
config/01_clustering/
config/02_create_expression/
config/03_export_summary/
config/04_plot_dotplot/
```

The `sankey/` configs are kept here because they are support utilities for old
label export and mapping rather than one of the core `00`-`04` pipeline steps.

