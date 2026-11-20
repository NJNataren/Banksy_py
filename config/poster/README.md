# Poster helpers

Configs for `helper_scripts/poster/` (HCA 2026 poster figures), moved here from
`config/08_cell_type_annotation/vbct/` on 2026-11-20 (post-08 cleanup,
phase 3). The helpers read only `cell_type`, UMAP embeddings, and marker
expression from the script 08 annotated H5ADs (h5py, read-only).

- `vbct/poster_cohort_umap_grid.json`: UMAP grid; lists samples in poster order
  with poster titles and optional point sizes.
- `vbct/poster_cohort_marker_dotplot.json`,
  `vbct/poster_cohort_marker_dotplot_rdbu_large_font.json`: cohort marker dot
  plot (equal sample weight, per-gene min-max scaling); all cohort samples.

Samples, object paths, UMAP keys, and palette come from the cohort config
(`config/cohorts/vbct/min15_area95_selected_resolutions.json`). Figures go to
`data/xenium/output/vbct/poster/min15_area95_selected_resolutions/`, dot plot
tables to `data/xenium/processed/vbct/poster/min15_area95_selected_resolutions/cohort_marker_dotplot/`.
Figures made for the poster before 2026-11-20 remain in
`figures/conference_min15_area95_selected_resolutions_poster/`.

```bash
P=config/poster/vbct
python helper_scripts/poster/plot_cohort_umap_grid_from_config.py --config $P/poster_cohort_umap_grid.json
python helper_scripts/poster/plot_cohort_marker_dotplot_from_config.py --config $P/poster_cohort_marker_dotplot_rdbu_large_font.json
```

Poster decisions and status: `context/tasks/hca_poster.md`.
