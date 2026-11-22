# Cohort configs

One file per cohort and QC branch holds the clustering facts that script 08 and
every later step need: sample names, response, tissue, branch, selected
resolution, the script 06 object (`recluster_h5ad`), cluster column
(`groupby`), UMAP key, and coordinate keys. Added 2026-11-19 (post-08 cleanup,
phase 1). Since 2026-11-22 the labels (mapping, palette) and the script 08
object paths belong to versioned annotation sets in
`config/annotations/<project>/<cohort file stem>/vN/` (README there), so one
cohort can carry several labellings.

Annotation sets point here with `"cohort"`; step configs name an annotation
set. Step configs may list `samples` to choose and order samples for that step
and to add plot-only fields (`title`, `facet`, `display_label`, `point_size`);
they may not redefine a cohort fact. Loader and rules: `xenium_workflow/cohorts.py` at the
repository root.

Paths here are relative. Config files resolve against the repository root;
data paths resolve against each step's `--input-root` or `--output-root`
(both default to the repository root).

Readers: script 08 (annotate and replot configs), script 09, script 10
(`10_cell_type_spatial_analysis.py` and `run_10_cell_type_spatial_analysis.sl`),
script 10a, both
`helper_scripts/spatial_analysis/` helpers, and both `helper_scripts/poster/`
helpers. The conference-figure helpers still keep their own sample list.
