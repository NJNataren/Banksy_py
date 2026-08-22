# Local Run Helpers

These scripts are intended for local or login-node helper runs rather than Slurm submission. They resolve the repository root from their own location before running, so they can be launched from the repo root or directly from this directory.

- `run_01a_clustree_ptmt_pc55_local.sh`: local clustree QC loop for PTMT PC55 script 00 cluster CSVs.
- `run_01a_clustree_from_config_local.sh`: config-driven local clustree runner; defaults to VBCT small post-07 recluster plots.
- `test_local_dotplot_export_CK_skin_res_0p5.sh`: local legacy dotplot export smoke test for CK skin resolution 0.5.
