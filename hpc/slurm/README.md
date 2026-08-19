# Slurm Wrappers

Run these wrappers from the repository root on the HPC, for example:

```bash
sbatch hpc/slurm/run_00_xenium_clustering.sl
```

Most wrappers change to `REPO_DIR` before running the Python workflow scripts, so moving them under `hpc/slurm/` does not change script-relative config, data, output, or log paths. `REPO_DIR` defaults to `/scratchdata1/users/a1210419/Banksy_py` where the wrapper supports overrides.

## Workflow Wrappers

- `run_00_xenium_clustering.sl`: script 00 clustering array, configurable with `CONFIG_DIR`.
- `run_00_xenium_clustering_large.sl`: script 00 large VBCT clustering wrapper.
- `run_00a_pca_scree_array.sl`: script 00a PCA scree backfill array.
- `run_01_xenium_QC_array_vbct.sl`: script 01 VBCT QC array, configurable with `CONFIG_DIR`.
- `run_01_xenium_QC_array_ptmt.sl`: script 01 PTMT/PC55 QC array, configurable with `CONFIG_DIR`.
- `run_02_create_expression_adata_with_banksy_clusters.sl`: legacy/helper script 02 expression object creation.
- `run_03_xenium_dotplot_export_from_clean_script00_config.sl`: preferred script 03 clean-object dotplot export.
- `run_03_xenium_dotplot_export_from_config.sl`: legacy script 03 export path.
- `run_04_xenium_multi_sample_dotplot_from_config.sl`: script 04 multi-sample dotplot plotting.
- `run_05_apply_qc_filters_for_reclustering.sl`: script 05 QC filter annotation.
- `run_06_recluster_qc_annotated_with_banksy.sl`: script 06 BANKSY reclustering from QC-annotated objects; defaults to `config/06_recluster_qc_annotated/vbct`.
- `run_07_squidpy_recluster_spatial_analysis.sl`: script 07 Squidpy spatial analysis on reclustered objects; defaults to `config/07_squidpy_recluster_analysis/vbct`.

## Smoke-Test Wrappers

- `run_00_xenium_clustering_ck_skin_res_missing_dirs_test.sl`
- `run_00_xenium_clustering_clean_adata_ck_skin_res_test.sl`
- `run_01_xenium_QC_ck_skin_res_test.sl`
- `run_xenium_clustering_test_single.sl`
