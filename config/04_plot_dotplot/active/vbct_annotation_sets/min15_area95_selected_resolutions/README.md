# Per-sample cluster marker dot plots labelled by an annotation set

One config per cohort sample. Each plots only that sample's clusters at the
cohort's selected resolution (min15_area95 branch), labelled
`cluster <id> · <tier label> (n=<cells>)` and grouped by label in the set's
vocabulary order. Added 2026-11-22.

- Input: the existing script 03 filter-sensitivity exports in
  `data/xenium/processed/vbct/cross_sample_dotplot_exports/filter_sensitivity/`
  (rows with `sample == min15_area95`). Their clusters and cell counts were
  checked against v1's `cluster_annotation_counts.csv` for all eight samples
  (127 clusters, identical), so no new export is needed.
- Genes: the 89 genes kept by
  `03_xenium_gene_review_for_dotplot_manual_annotation_2026-06-11.csv`,
  grouped by `Tier_1_annotation`; canonical genes in red.
- Colour: mean expression z-scored per gene across the plotted sample's
  clusters, clipped at +/-2.5; dot size is percent expressing.
- Output: `data/xenium/output/vbct/cell_type_annotation/<cohort>/<version>/<tier>/cluster_marker_dotplots/`.

```bash
for f in config/04_plot_dotplot/active/vbct_annotation_sets/min15_area95_selected_resolutions/*.json; do
  python 04_plot_multi_sample_dotplot_from_config_local.py --config "$f"
done
# same plots with another set or tier (drafts allowed; the title shows the status):
#   ... --config "$f" --annotation config/annotations/vbct/min15_area95_selected_resolutions/v2 [--label-level tier_2]
```

The `annotation` block is read by `apply_annotation_labels` in
`04_plot_multi_sample_dotplot_from_config_local.py`.
