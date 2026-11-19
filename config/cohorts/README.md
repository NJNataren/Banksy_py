# Cohort configs

One file per cohort and QC branch holds the facts that script 08 and every later
step need: sample names, response, tissue, branch, selected resolution, the
script 06 (`recluster_h5ad`) and script 08 (`annotated_h5ad`) object paths,
cluster column (`groupby`), UMAP key, coordinate keys, the reviewed mapping
CSV, and the shared palette. Added 2026-11-19 (post-08 cleanup, phase 1).

Step configs point here with `"cohort": "config/cohorts/<project>/<file>.json"`.
They may list `samples` to choose and order samples for that step and to add
plot-only fields (`title`, `facet`, `display_label`, `point_size`); they may
not redefine a cohort fact. Loader and rules: `workflow_config.py` at the
repository root.

Paths here are relative. Config files resolve against the repository root;
data paths resolve against each step's `--input-root` or `--output-root`
(both default to the repository root).

Readers: script 08 (annotate and replot configs), script 09, script 10
(cell-type mode and `run_10_cell_type_spatial_analysis.sl`), both
`helper_scripts/spatial_analysis/` helpers, and both `helper_scripts/poster/`
helpers. The conference-figure helpers still keep their own sample list.
