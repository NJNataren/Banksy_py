# Script 04 Dotplot Configs

`active/` contains current runnable plot configs for the script 00 clean-object workflow.

`archive/` contains older dated runs, archive-label configs, and old HPC/legacy layouts retained for provenance.

## Active Layout

```text
active/
├── vbct_annotation_sets/min15_area95_selected_resolutions/   per-sample cluster plots labelled by an annotation set (README there)
├── vbct_filter_sensitivity/                                  QC-branch pseudo-sample plots
│   └── selected_resolutions/
├── vbct_recluster_qc_pass_only/
├── vbct_recluster_grouped_area95_trial/
├── vbct_clean_script00/
│   ├── multi_sample/
│   ├── per_sample_all_genes_all_res/
│   └── per_sample_canonical_markers_all_res/
├── ptmt_clean_script00/
│   └── multi_sample/
└── ptmt_pc55_clean_script00/
    ├── multi_sample/
    ├── per_sample_panel_all_res/run_1/
    └── per_sample_all_genes_all_res/
```

Use `04_plot_multi_sample_dotplot_from_config_local.py --config <config.json>` for these configs.
Script 03 and 04 are described in `README_00_04_XENIUM_ANALYSIS_WORKFLOW.md`;
the options below are newer than that manual.

## Script 03 Options (Clean-Object Exporter)

`03_export_dotplot_data_from_clean_script00_config.py` summarises `.raw` when
present, otherwise `.X`, per `.obs` grouping column. It attaches script 01
`new_labels` to `group_label` by default when `project` is set and
`config/01_QC/<project>/<sample>.json` exists. Turn that off with
`"use_qc_annotations": false`, or point elsewhere with `qc_config_dir` or a
per-sample/per-grouping `qc_config`. Disable it for recluster branches so
branch-specific numeric cluster IDs are kept.

Summary CSV columns used by script 04: `sample`, `resolution`, `cluster_id`,
`group_id`, `group_label`, `groupby`, `gene`, `marker_group`,
`mean_expression`, `percent_expressing`, `n_cells`, `expression_source`,
`adata_path`.

## Script 04 Options

- **Rows.** `sample_filters` selects samples and resolutions; a filter's
  `"resolution": "all"`, `"include_all_resolutions": true` or a `"resolutions"`
  list keeps several. Resolutions compare numerically (`0.7` == `0.70`).
  `sample_order` stacks pseudo-samples (e.g. QC branches);
  `y_label_template` formats row labels.
- **Genes.** `gene_file` (+ `gene_column`, `gene_group_column`) sets the genes
  and their groups; `gene_review_file` (+ `gene_review_keep_column`) drops
  genes. `gene_order_column` orders genes (e.g. `plot_order`);
  `gene_label_columns` + `gene_label_separator` build display labels such as
  `ANPEP - LE-KLK3; Mono`. `highlight_gene_file` + `highlight_gene_column`
  colour genes, styled by `figure.highlight_label_color`/`highlight_label_weight`.
  These travel together in the `GeneMetadata` dataclass.
- **Colour.** Mean expression is z-scored per gene across the plotted rows,
  clipped at `z_score_clip` (colour limits default to ±2.5). Dot size is
  `percent_expressing`.
- **Annotation sets.** An `annotation` block
  (`{"set": <set dir>, "sample": <sample>, "label_level": "tier_1",
  "source_rows": {"sample": "min15_area95"}}`) keeps that sample's rows at the
  cohort's clustering column and selected resolution, adds the tier label from
  the set's `mapping.csv`, and fails if any cluster is unmapped or any mapped
  cluster is missing. It groups rows by label, turns on `rows_top_to_bottom`,
  and derives the output path and title from the set. `--annotation <set dir>`
  and `--label-level <tier>` replot the same config for another version or tier.

## Annotation Review Metadata

Script 04 can join one or more cluster-level metadata tables without modifying
the script 03 expression exports. Each entry must have one row per configured
join key:

```json
"row_metadata_files": [
  {
    "file": "path/to/cluster_review.csv",
    "join_columns": ["sample", "cluster_id"],
    "columns": ["proposed_label", "confidence"]
  }
],
"y_label_template": "{sample} | c{cluster_id} | {proposed_label} | {confidence}"
```

`row_filters` accepts allowed values by metadata column. Set
`row_filter_mode` to `"any"` to retain rows matching at least one filter;
the default `"all"` requires every filter. Use `row_order_columns`,
`row_order_values`, and `row_separator_column` to group review rows by
annotation rather than the default sample/resolution/cluster order.

For marker tables with a rank that restarts within each cell type, set
`group_gene_order` to `true` and provide `gene_group_order`. Script 04
then keeps marker groups together and applies `gene_order_column` within each
group.

The PTMT PC55 run 1 broad-analysis overview and priority-review examples are:

```text
active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_run1_r1p5_broad_marker_review_local.json
active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_run1_r1p5_broad_marker_priority_review_local.json
```

The corresponding minor-class examples are:

```text
active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_run1_r1p5_minor_class_marker_review_local.json
active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_run1_r1p5_minor_class_marker_priority_review_local.json
```

All four automated broad/minor marker-review configs write their current PNG
outputs to `figures/dotplots/automated_marker_selection/`.

For PTMT PC55 run 1 per-sample panel dotplots across all available resolutions, run:

```bash
for cfg in config/04_plot_dotplot/active/ptmt_pc55_clean_script00/per_sample_panel_all_res/run_1/*.json; do
  conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
    --config "$cfg"
done
```
