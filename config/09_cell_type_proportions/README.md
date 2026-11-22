# Cell-type proportions

Run locally from Banksy_py/ in an environment with NumPy, pandas, and Matplotlib:

```bash
python 09_plot_cell_type_proportions.py --config config/09_cell_type_proportions/vbct/min15_area95_selected_resolutions.json
```

Use --validate-only to check inputs without writing files, or --overwrite to
replace existing outputs. Since 2026-11-22 the config names an annotation set
(`"annotation"`) and tier (`"label_level"`, default `tier_1`); `--annotation`
and `--label-level` override them, and `--allow-draft` permits a draft set. The
input is script 08's `cell_type_counts_<tier>.csv` for that set (under
`--input-root`); branch and resolution come from the set's cohort and the
palette from the set (see `config/annotations/README.md`).
`expected_total_cells` is an optional extra check.

The script checks sample/branch/resolution provenance, unique sample/type rows,
nonnegative integer counts, expected total cells, and the shared palette.
Missing types are filled with zero; each denominator is the total annotated
post-filtered population of that sample. Sample order follows the config;
bottom-to-top stack and legend order follow the palette. Do not pool QC branches.

Outputs are PNG, SVG, and PDF in
`data/xenium/output/<project>/cell_type_proportions/<cohort>/<version>/<tier>/` and the
proportion CSV in the same path under `data/xenium/processed/`. Unversioned
outputs in `.../cell_type_proportions/min15_area95_selected_resolutions/` predate
2026-11-22 (v1 labels); earlier ones are in
`figures/_pre_post08_cleanup_2026-11-20/cell_type_proportions/`.
The supplied config uses 300 DPI and editable vector text. Source count tables,
annotated objects, and UMAPs are preserved. Plot typography is configurable;
the supplied configs use 14-point legend text, 15-point legend titles and axis
labels, 13-point tick labels, and 16-point tissue-facet titles. They omit the
overall title and caption while retaining the shared cell-type legend and the
tissue headings needed to identify facets.

## Tissue facets

Use `vbct/min15_area95_selected_resolutions_by_tissue.json` with the same script
for Skin, Brain, Lung, and Gastrointestinal facets. Per-sample `facet` and
`display_label` fields in this config define grouping and labels (facets are
plot groupings, separate from the cohort's `tissue`); `plot.facet_order` defines
panel order. Panel widths follow sample counts, maintaining equal-width bars.
Outputs use the separate stem `cell_type_proportions_by_tissue`.
