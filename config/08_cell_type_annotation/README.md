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
to replace existing outputs.

## HPC execution

Submit script 08 from the HPC repository root instead of running it on a login
node:

```bash
cd /scratchdata1/users/a1210419/Banksy_py
mkdir -p logs
sbatch hpc/slurm/run_08_cell_type_annotation.sl
```

The wrapper requests 8 CPUs, 180 GB RAM, and 8 hours. It overrides both input
and output roots with the HPC repository, avoiding the local mounted path stored
in the active config.

This is deliberately one job rather than a Slurm array. Script 08 writes two
shared count tables and a combined eight-sample UMAP after processing all
samples; independent per-sample tasks would race on those shared outputs and
produce incomplete summaries.

The default is fail-safe: existing outputs cause the script to stop. To replace
partial or previously completed script 08 outputs intentionally, submit:

```bash
sbatch --export=ALL,OVERWRITE=1 hpc/slurm/run_08_cell_type_annotation.sl
```

After script 08 finishes successfully, validate script 10 inputs and then
submit `hpc/slurm/run_10_cell_type_spatial_analysis.sl`.
