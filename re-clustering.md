# QC-Filtered BANKSY Reclustering Notes

## Current State

- Script `05_apply_qc_filters_for_reclustering.py` applies reviewed QC filters to the script 01 QC-annotated clean expression object.
- Default filtered output path: `data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_qc_filtered_qc_v1.h5ad`.
- The duplicated `filtered_filtered` filename issue was fixed by using `adata_expression_clean_<sample>_qc_<output_label>.h5ad`.
- Existing CK skin config: `config/05_apply_qc_filters/vbct/CK_skin_res.json`.
- Existing Slurm wrapper: `run_05_apply_qc_filters_for_reclustering.sl`.

## Proposed Next Script

Create:

```text
06_recluster_qc_filtered_with_banksy.py
```

Create configs under:

```text
config/06_recluster_qc_filtered/
```

The script should rerun the BANKSY clustering portion of `00_xenium_clustering_clean_adata.py` on the QC-filtered clean expression object produced by script 05.

## Recommended Design

Keep these script 00 behaviours:

- config-driven `project`, `dataset_name`, `pc_label`, `lambda_label`, `res_label`, `nbr_weight_decay`, `coord_keys`, `max_workers`, and optional `scree_n_pcs`;
- directory creation for `data/xenium/processed/<project>/<sample>/` and `data/xenium/output/<project>/<sample>/`;
- coordinate setup from `obs["x"]` and `obs["y"]` into `obsm["xy"]` and `obsm["spatial"]`;
- raw/count layer preservation before normalization if needed;
- PCA scree QC before BANKSY matrix expansion;
- `initialize_banksy`, `generate_banksy_matrix`, nonspatial comparison, `pca_umap`, and parallel Leiden clustering;
- BANKSY result plots via `plot_results`;
- export of `adata_spatial_*`, `banksy_dict.pkl.gz`, `results_df_*`, and cell-by-resolution cluster CSV;
- transfer of BANKSY labels and embeddings back onto a clean expression AnnData;
- marker ranking on the clean expression matrix, not BANKSY-expanded `.X`.

Skip these script 00 behaviours:

- raw Xenium AnnData loading from `data/xenium/raw_data`;
- zero-count plotting/filtering;
- general script 01 QC diagnostics.

## Output Naming

Keep the QC-filter label in every reclustering output to avoid overwriting original script 00 outputs.

Example input:

```text
adata_expression_clean_CK_skin_res_qc_filtered_qc_v1.h5ad
```

Example outputs:

```text
adata_spatial_CK_skin_res_qc_filtered_qc_v1_0p5.h5ad
results_df_CK_skin_res_qc_filtered_qc_v1_pc30_nc0.20_r0.50_0.60_0.70.csv
CK_skin_res_qc_filtered_qc_v1_cell_cluster_id_across_clustering_res_0.50_0.60_0.70.csv
adata_expression_clean_CK_skin_res_qc_filtered_qc_v1_with_banksy_clusters_0.50_0.60_0.70.h5ad
```

## Draft Config Shape

```json
{
  "project": "vbct",
  "dataset_name": "CK_skin_res",
  "input_label": "filtered_qc_v1",
  "pc_label": "30",
  "lambda_label": "0.20",
  "res_label": ["0.50", "0.60", "0.70", "0.80", "0.90", "1.00", "1.10"],
  "nbr_weight_decay": "scaled_gaussian",
  "coord_keys": ["x", "y", "xy"],
  "max_workers": 8,
  "scree_n_pcs": 75,
  "run_neighbour_analysis": true
}
```

## Neighbour Analysis To Add

The user wants the neighbour analysis that is currently commented out in script 00 to be included in the reclustering workflow.

Suggested structure:

```text
06_recluster_qc_filtered_with_banksy.py
├── parse config and resolve paths
├── load QC-filtered clean expression AnnData
├── prepare coordinates and expression layers
├── optional PCA scree QC
├── run BANKSY graph/matrix construction
├── run PCA/UMAP and Leiden clustering
├── save BANKSY spatial objects and clean clustered object
├── run marker ranking on clean expression object
└── run neighbour analysis
    ├── choose cluster label columns / resolution(s)
    ├── build spatial nearest-neighbour graph from filtered cells
    ├── summarize neighbour composition by cluster or annotation
    ├── export neighbour count/proportion tables
    └── save neighbour composition plots
```

Open question for the next session: identify the exact commented neighbour-analysis block in `00_xenium_clustering_clean_adata.py` or related notebooks and decide whether it should run for all resolutions or only a selected annotation resolution.
