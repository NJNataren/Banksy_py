# Poster helpers

Configs for `helper_scripts/poster/` (HCA 2026 poster figures), moved here from
`config/08_cell_type_annotation/vbct/` on 2026-11-20 (post-08 cleanup,
phase 3). The helpers read only the chosen tier's label column (e.g.
`cell_type_tier_1`), UMAP embeddings, and marker expression from the script 08
annotated H5ADs (h5py, read-only), after checking each object's annotation stamp.

- `vbct/poster_cohort_umap_grid.json`: UMAP grid; lists samples in poster order
  with poster titles and optional point sizes.
- `vbct/poster_cohort_marker_dotplot.json`,
  `vbct/poster_cohort_marker_dotplot_rdbu_large_font.json`: cohort marker dot
  plot (equal sample weight, per-gene min-max scaling); all cohort samples.

Since 2026-11-22 each config names an annotation set (`"annotation"`), a tier
(`"label_level"`), and an `output_subdir`; samples, object paths, UMAP keys,
and palette come from the set and its cohort (`config/annotations/README.md`).
`--annotation`/`--label-level` override the config; `--allow-draft` permits a
draft set. Marker keys in the dot plot config must be labels of the chosen tier.
Figures go to `data/xenium/output/<project>/poster/<cohort>/<version>/<tier>/<output_subdir>/`,
dot plot tables to the same path under `data/xenium/processed/`.
Figures made for the poster before 2026-11-20 remain in
`figures/conference_min15_area95_selected_resolutions_poster/`.

```bash
P=config/poster/vbct
python helper_scripts/poster/plot_cohort_umap_grid_from_config.py --config $P/poster_cohort_umap_grid.json
python helper_scripts/poster/plot_cohort_marker_dotplot_from_config.py --config $P/poster_cohort_marker_dotplot_rdbu_large_font.json
```

Poster decisions and status: `context/tasks/hca_poster.md`.
