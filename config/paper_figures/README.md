# Paper figure configs

Configs for `helper_scripts/paper_figures/`. Like the other post-08 steps, each
config names an annotation set (`"annotation"`) and a tier (`"label_level"`);
`--annotation`/`--label-level` override them, `--allow-draft` permits a draft
set, and `--input-root`/`--output-root` default to the repository root.
Figures go to `data/xenium/output/<project>/paper_figures/<cohort>/<version>/<tier>/<output_subdir>/`,
tables to the same path under `data/xenium/processed/`. Scripts are named by
content, not figure number.

## Cell-type marker dot plots (issue #64, added 2026-11-24)

Two configs, same markers and inputs:

- `vbct/cell_type_marker_dotplot.json`: one row per sample x cell type, rows
  grouped by sample in the config's `samples` order (`tissue_group` draws the
  thicker separators). Colour: z-score per gene within each sample.
- `vbct/cell_type_marker_dotplot_cohort.json`: one row per cell type; each
  sample's values averaged with equal sample weight. Colour: z-score per gene
  across cell types.
- `vbct/cell_type_marker_dotplot_by_sample_global_zscore.json`: supplementary
  detection panel; by-sample rows with each gene z-scored across all 50 rows,
  so colour is comparable between samples. It shows sample-level detection
  differences (MF near white throughout) and spillover (melanoma markers in
  HW immune rows); not for justifying labels.

Options: `"rows"`: `sample_cell_type` or `cell_type`; `"scaling"`:
`{"method": "zscore" | "minmax", "within_sample": bool, "clip": 2.5}`. The
z-score follows Seurat's DotPlot (per gene, sample SD, clipped at +/-2.5); the
mean being scaled is the mean of log-normalised expression (Seurat averages
expm1 values). With few rows the reachable z-score is capped (about +1.8 for 5
cell types, +2.3 for 7).

- Input: the script 03 per-cluster exports listed in `exports` (only
  `source_rows` and the cohort's groupby and resolution are used), labelled
  with the set's mapping. Clusters of one cell type are merged by cell-count
  weighting, which is exact because script 03 averages over all cells
  including zeros. No H5ADs or mount are needed. Merged counts are checked
  against script 09's counts for the same set when those exist.
- Dot size is percent expressing. Colour is not comparable between samples in
  the per-sample plot; min-max across all rows (`within_sample: false` with
  sample rows) was tried first and made low-detection samples (e.g. MF) pale.
- Tables (shared by both configs): `cell_type_expression_by_sample.csv` (all
  panel genes), `marker_specificity_per_sample.csv` and
  `marker_specificity_summary.csv` (each cell type against all other cells of
  the same sample; descriptive, no thresholds), and `<output_stem>_values.csv`
  (the plotted values).

```bash
for c in cell_type_marker_dotplot cell_type_marker_dotplot_cohort; do
  python helper_scripts/paper_figures/plot_cell_type_marker_dotplot_from_config.py \
    --config config/paper_figures/vbct/$c.json
done
```
