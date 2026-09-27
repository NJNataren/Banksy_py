# Script 04 Dotplot Configs

`active/` contains current runnable plot configs for the script 00 clean-object workflow.

`archive/` contains older dated runs, archive-label configs, and old HPC/legacy layouts retained for provenance.

## Active Layout

```text
active/
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
