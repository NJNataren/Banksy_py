# Cell-type proportions

Run locally from Banksy_py/ in an environment with NumPy, pandas, and Matplotlib:

```bash
python 09_plot_cell_type_proportions.py --config config/09_cell_type_proportions/vbct/min15_area95_selected_resolutions.json
```

Use --validate-only to check inputs without writing files, or --overwrite to
replace existing outputs. Paths are relative to the project root inferred from
the config project_root_parent_depth (3 for the supplied config).

The script checks sample/branch/resolution provenance, unique sample/type rows,
nonnegative integer counts, expected total cells, and the shared palette.
Missing types are filled with zero; each denominator is the total annotated
post-filtered population of that sample. Sample order follows the config;
bottom-to-top stack and legend order follow the palette. Do not pool QC branches.

Outputs are PNG, SVG, PDF, and a proportion CSV in figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/.
The supplied config uses 300 DPI and editable vector text. Source count tables,
annotated objects, and UMAPs are preserved. Plot typography is configurable;
the supplied configs use 14-point legend text, 15-point legend titles and axis
labels, 13-point tick labels, and 16-point tissue-facet titles. They omit the
overall title and caption while retaining the shared cell-type legend and the
tissue headings needed to identify facets.

## Tissue facets

Use `vbct/min15_area95_selected_resolutions_by_tissue.json` with the same script
for Skin, Brain, Lung, and Gastrointestinal facets. Per-sample `tissue` and
`display_label` fields define grouping and labels; `plot.facet_order` defines
panel order. Panel widths follow sample counts, maintaining equal-width bars.
Outputs use the separate stem `cell_type_proportions_by_tissue`.
