# Task: PTMT Dot Plots

PTMT (PC35) and PTMT PC55 marker dot plots made with script 03 (export, HPC)
and script 04 (render, local). Split from `dotplot.md` on 2026-11-24; the full
earlier text is in `context/tasks/dotplot/archive/dotplot_to_2026-11-24.md`. Cluster
annotation and shared-marker selection for PC55 run 1 are in
`marker_selection.md`. Script 03/04 options: `config/04_plot_dotplot/README.md`.

Keep `ptmt` and `ptmt_pc55` outputs separate.

## Current State

**Marker files** (`data/xenium/raw_data/gene_markers/ptmt/`):
- `ptmt_combined_xenium_panel_markers.csv`: one row per panel gene, built from
  `hImmuno_submitted_to_designer_tier1.csv` and
  `hProstate_100g_with_Tier_1_annotations.csv`; genes in both panels have
  `panel_source`/`Tier_1_annotation` joined by `;`. Script 03 uses
  `"gene_column": "Gene"`, `"marker_group_column": "Tier_1_annotation"`.
- `ptmt_combined_xenium_panel_gene_review_for_dotplot.csv`: review filter; its
  `colour_canonical` column drives red label highlighting.
- `2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv`:
  harmonised scRNA-seq markers (`plot_order`, `gene`, `major_class`,
  `minor_class`, `cell_types`).

**Script 03 exports** (each contains all resolutions for the panel genes):

| Branch | Config dir | Output |
|---|---|---|
| PTMT PC35 | `config/03_export_summary/active/ptmt_clean_script00/ptmt_panel_pc35/` | `data/xenium/processed/ptmt/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv` |
| PC55 panel | `config/03_export_summary/active/ptmt_pc55_clean_script00/ptmt_panel/` | `data/xenium/processed/ptmt_pc55/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv` |
| PC55 all genes | `config/03_export_summary/active/ptmt_pc55_clean_script00/all_genes/` | only needed for literally every gene in the object, not the panel |

```bash
# HPC (Nathalie runs), PC35 example:
sbatch --array=0-23%3 \
  --export=CONFIG_DIR=config/03_export_summary/active/ptmt_clean_script00/ptmt_panel_pc35 \
  hpc/slurm/run_03_xenium_dotplot_export_from_clean_script00_config.sl
```

**Script 04 configs** (local):
- PC35, all 24 samples at r0.7:
  `config/04_plot_dotplot/active/ptmt_clean_script00/multi_sample/ptmt_panel_multi_sample_local.json`
  → `figures/dotplots/local_ptmt_panel_multi_sample_r0p70_zscore_dotplot.png`.
- PC55, in `config/04_plot_dotplot/active/ptmt_pc55_clean_script00/`:
  - `multi_sample/ptmt_pc55_panel_multi_sample_local.json`: all samples, panel genes.
  - `multi_sample/ptmt_pc55_minor_class_run_1_r_0p70_local.json`,
    `..._r1p00_gene_cell_types_labels_local.json`, `..._r1p50_local.json`,
    `..._r1p50_gene_cell_types_labels_local.json`: eight run 1 samples, Kat's
    scRNA-seq markers grouped by `minor_class`, ordered by `plot_order`, no gene
    clustering; label variants show `gene - cell_types`.
  - `multi_sample/ptmt_pc55_run1_r1p5_{broad,minor_class}_marker_{review,priority_review}_local.json`:
    marker-selection review plots (see `marker_selection.md`).
  - `per_sample_panel_all_res/run_1/*.json`, `per_sample_all_genes_all_res/*.json`.
- Superseded PC55 trials: `config/04_plot_dotplot/archive/ptmt_pc55_config_trials/`.
- All render at 150 DPI.

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r1p00_gene_cell_types_labels_local.json
```

## Decisions

- 2026-08: PC55 run 1 plots use Kat's harmonised scRNA-seq markers, grouped by
  `minor_class`, ordered by `plot_order`; gene clustering off.

## Next Step

Review the `cell_types` column of Kat's scRNA-seq marker file and choose one
display label per gene (some rows list several, `;`-separated). Store it in a
new display-label column (do not overwrite the source columns) and point
`gene_label_columns` in the r1.00 label config at it. Not started; lower
priority than the VBCT work in `dotplot.md`.

## History

- To 2026-11-24: `context/tasks/dotplot/archive/dotplot_to_2026-11-24.md` (config
  reorganisation, PC55 run 1 plots, superseded dot-spacing task).
