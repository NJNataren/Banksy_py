# Script 08: cell-type annotation and UMAPs

`08_annotate_cell_types_and_plot_umaps.py` applies an annotation set (see
`config/annotations/README.md`) to clean script 06 recluster objects. It writes
one label column per tier (`cell_type_tier_1`, `cell_type_tier_2`, ...), a
provenance stamp in `uns["cell_type_annotation"]`, per-tier Scanpy UMAPs,
count tables, a change table against the parent set, and a run manifest.
Source objects, original cluster columns, expression and embeddings are
preserved. Since 2026-11-22 the mapping and palette live in the annotation set,
not here, and every output path is derived from the set
(`xenium_workflow/annotation_sets.py`); the old single `cell_type` column is retired.

Configs here hold plot settings only:

- `vbct/min15_area95_selected_resolutions.json`: annotate and plot. Names the
  set with `"annotation"`; optional `plot_levels` limits which tiers are plotted
  (default: all).
- `vbct/min15_area95_selected_resolutions_clean_umaps.json`: replot clean
  presentation UMAPs (no titles, axis labels, or external legend; 12-point
  on-data labels) from the set's existing annotated objects
  (`"input_key": "annotated_h5ad"`, `output_subdir: clean_scanpy`). Writes no
  H5ADs or tables.

`--annotation <set dir>` runs either config against another version. Figures
are PNG, SVG and PDF; dense points are rasterized so Inkscape can edit text and
layout (`rasterized_points: false` only when editable points are required).

Locally, script 06 inputs are on the Phoenix mount:

```bash
A=config/annotations/vbct/min15_area95_selected_resolutions
conda run -n banksy python 08_annotate_cell_types_and_plot_umaps.py \
  --config config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions.json \
  --annotation $A/v2 --input-root ~/mnt/phoenix_scratch/Banksy_py --validate-only
```

Drop `--validate-only` to write outputs. `--sample <sample>` restricts samples
(tables then cover only those samples), `--no-write-annotated` skips H5ADs, and
`--overwrite` replaces existing outputs (normal while a set is a draft).

Outputs (under `--output-root`):

```text
data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/
  <sample>/adata_expression_clean_<sample>_<branch>_cell_type_annotated_<version>_r<res>.h5ad
  cluster_annotation_counts.csv     one row per source cluster, all tiers, n_cells
  cell_type_counts_<tier>.csv       per sample and label (read by script 09)
  changes_vs_<parent>.csv           clusters whose labels changed since the parent set
  run manifest (JSON)               inputs, outputs, set hash, cell totals
data/xenium/output/<project>/cell_type_annotation/<cohort>/<version>/<tier>/[<output_subdir>/]
  per_sample_scanpy/<sample>_cell_type_<tier>_umap_scanpy.{png,svg,pdf}
  combined_scanpy_cell_type_<tier>_umaps.{png,svg,pdf}
```

Outputs from before 2026-11-22 (unversioned `.../min15_area95_selected_resolutions/`
folders and `figures/_pre_post08_cleanup_2026-11-20/cell_type_annotation/`)
are from the same labels as v1 but carry no stamp.

## HPC execution

Submit from the HPC repository root (one job, not an array: script 08 writes
shared tables and combined UMAPs after all samples):

```bash
cd /scratchdata1/users/a1210419/Banksy_py
mkdir -p logs
sbatch hpc/slurm/run_08_cell_type_annotation.sl
# another version, or a rerun of a draft:
sbatch --export=ALL,ANNOTATION=config/annotations/vbct/min15_area95_selected_resolutions/<version>,OVERWRITE=1 \
  hpc/slurm/run_08_cell_type_annotation.sl
```

The wrapper requests 8 CPUs, 180 GB RAM, and 8 hours and passes the HPC
repository as input and output root. Existing outputs stop the run unless
`OVERWRITE=1`. After script 08 finishes (and the set is released), validate
and submit script 10 (`config/10_cell_type_spatial_analysis/vbct/README.md`).
