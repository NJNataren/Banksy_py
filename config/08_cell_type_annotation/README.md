# Script 08: reviewed cell-type annotation and UMAPs

`08_annotate_cell_types_and_plot_umaps.py` applies a reviewed cluster mapping to
clean script 06 recluster objects. It preserves source objects and original
cluster columns, writes separate annotated H5ADs and count audits, and uses the
same saved BANKSY coordinates for Matplotlib and Scanpy UMAP plots.

The active VBCT config uses Scanpy and writes PNG, SVG, and PDF. Dense point
clouds are rasterized inside SVG/PDF so Inkscape can edit vector text, labels,
legend, and layout without creating files containing millions of vector marks.
Set `rasterized_points` to `false` only when individual editable points are
required and the resulting file size is acceptable.

Validate without writing outputs:

```bash
conda run -n banksy python 08_annotate_cell_types_and_plot_umaps.py \
  --config config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions.json \
  --validate-only
```

Run all samples by removing `--validate-only`. Use `--sample <sample>` for a
subset, `--no-write-annotated` for figures and tables only, and `--overwrite`
to replace existing outputs. On HPC, pass
`--input-root /scratchdata1/users/a1210419/Banksy_py`.
